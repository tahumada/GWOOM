#!/usr/bin/env python
"""
boom_filter_cone_search.py
==========================

Run a BOOM filter (a MongoDB aggregation pipeline such as ztf_basic.json) on
the alerts inside ONE cone -- e.g. one tile from tile_skymap.py -- in a time
window after a GW event.

Uses BOOM's own filter-test endpoint, POST /filters/test, so the pipeline runs
exactly as BOOM's filter workers would run it: BOOM adds the auxiliary-data
lookups the filter needs (here `cross_matches.milliquas_v8`) and the ZTF
programid permissions. The cone is sent as a MOC (`moc_ascii`), which BOOM
matches on its `coordinates.hpx` index together with the `candidate.jd`
window. The MOC slightly overfills the cone at its edge, so results are cut to
the exact radius afterwards.

Server-side limits: each request may cover at most 7 days of JD, so longer
windows are split into 7-day chunks automatically.

Credentials (same as blastwave):
    export BOOM_API_USER="..."
    export BOOM_API_PASSWORD="..."      # or put both in a .env file

Usage
-----
    python boom_filter_cone_search.py --filter ztf_basic.json \\
        --tiles outputs/tiles.csv --tile-id 0 --t0 2026-10-12T03:14:15 --dt-days 7

    python boom_filter_cone_search.py --filter ztf_basic.json \\
        --ra 160.1 --dec 14.26 --radius 2 --t0 2461325.63 --programids 1 2 3
"""

from __future__ import annotations

import argparse
import copy
import csv
import json
import os
import sys

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.time import Time
from mocpy import MOC

from blastwave import BOOMClient

MAX_JD_WINDOW = 7.0          # enforced by /filters/test
MAX_MOC_RANGES = 4500        # server limit is 5000 HEALPix index ranges
EXTRA_FIELDS = ["candidate.ra", "candidate.dec", "candidate.magpsf",
                "candidate.sigmapsf", "candidate.fid", "candidate.band"]
ZTF_BANDS = {1: "g", 2: "r", 3: "i"}


def load_credentials() -> None:
    """Load BOOM credentials from a .env file.

    blastwave calls load_dotenv() itself, but when run from a script that
    searches upward from blastwave's *installed* location, not your folder,
    so your .env is never found. Load it here first: from the current
    directory (or its parents), then next to this script.

    blastwave reads BOOM_API_USER / BOOM_API_PASSWORD. Its README mentions
    BOOM_USERNAME / BOOM_PASSWORD instead, so accept those as a fallback."""
    from dotenv import find_dotenv, load_dotenv

    here = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".env")
    for path in (find_dotenv(usecwd=True), here):
        if path and os.path.exists(path):
            load_dotenv(path, override=False)
    for new, old in (("BOOM_API_USER", "BOOM_USERNAME"), ("BOOM_API_PASSWORD", "BOOM_PASSWORD")):
        if not os.getenv(new) and os.getenv(old):
            os.environ[new] = os.environ[old]
    missing = [v for v in ("BOOM_API_USER", "BOOM_API_PASSWORD") if not os.getenv(v)]
    if missing:
        sys.exit(f"Missing {', '.join(missing)}: put them in a .env file in the folder you "
                 f"run from (or next to this script), or export them in your shell.")


# --------------------------------------------------------------------------- #
# Inputs
# --------------------------------------------------------------------------- #
def read_tile(path: str, tile_id: int) -> tuple[float, float, float]:
    """(ra, dec, radius) in degrees for one row of tile_skymap.py's CSV."""
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            if int(row["tile_id"]) == tile_id:
                return float(row["ra_deg"]), float(row["dec_deg"]), float(row["radius_deg"])
    raise ValueError(f"tile_id {tile_id} not found in {path}")


def to_jd(t: str) -> float:
    """ISO string or JD number -> JD."""
    try:
        return Time(float(t), format="jd").jd
    except ValueError:
        return Time(t).jd


def load_pipeline(path: str) -> list[dict]:
    """Load a BOOM filter (a JSON list of aggregation stages) from a file."""
    with open(path) as f:
        pipeline = json.load(f)
    if isinstance(pipeline, dict) and isinstance(pipeline.get("pipeline"), list):
        pipeline = pipeline["pipeline"]          # tolerate {"pipeline": [...]} wrappers
    if not isinstance(pipeline, list) or not all(isinstance(s, dict) for s in pipeline):
        raise ValueError(f"{path} is not a BOOM filter: expected a JSON list of stages")
    return prepare_pipeline(pipeline)


def prepare_pipeline(pipeline: list[dict], extra_fields: list[str] = EXTRA_FIELDS) -> list[dict]:
    """Copy of `pipeline` whose $project stages keep the fields we need to
    report and to apply the exact cone cut (most filters keep neither RA nor Dec).

    In MongoDB a $project is an *exclusion* only if every field (except _id)
    is 0/false; those already keep RA/Dec, so only inclusion projections
    (including ones with computed fields) are extended."""
    pipeline = copy.deepcopy(pipeline)
    for stage in pipeline:
        proj = stage.get("$project")
        if not proj:
            continue
        exclusion = all(v in (0, False) for k, v in proj.items() if k != "_id")
        if exclusion or "candidate" in proj:     # "candidate" kept whole already
            continue
        for fld in extra_fields:
            proj.setdefault(fld, 1)
    return pipeline


def cone_moc(ra: float, dec: float, radius_deg: float) -> tuple[str, int]:
    """IVOA ASCII MOC covering the cone, at the finest depth (<=12, as BOOM uses
    for its own cone searches) that stays under the server's range limit."""
    for depth in range(12, 2, -1):
        moc = MOC.from_cone(lon=ra * u.deg, lat=dec * u.deg, radius=radius_deg * u.deg,
                            max_depth=depth, delta_depth=2)
        if len(moc.to_depth29_ranges) <= MAX_MOC_RANGES:
            return moc.to_string(format="ascii"), depth
    raise ValueError("Cone too large to express as a MOC within the server limit")


# --------------------------------------------------------------------------- #
# Query
# --------------------------------------------------------------------------- #
def run_filter_in_cone(client: BOOMClient, pipeline: list[dict], ra: float, dec: float,
                       radius_deg: float, jd_start: float, jd_end: float,
                       programids: list[int] | None, survey: str = "ZTF",
                       limit: int = 10_000, verbose: bool = True) -> list[dict]:
    """Run `pipeline` on `survey` alerts inside the cone and JD window.
    Returns the raw alert documents, cut to the exact cone radius.
    programids=None sends no permissions (surveys without programids, e.g. LSST)."""
    moc_ascii, depth = cone_moc(ra, dec, radius_deg)
    results = []
    edges = np.arange(jd_start, jd_end, MAX_JD_WINDOW).tolist() + [jd_end]
    for j0, j1 in zip(edges[:-1], edges[1:]):
        body = {
            "pipeline": pipeline,
            "survey": survey,
            "permissions": {survey: programids} if programids else {},
            "moc_ascii": moc_ascii,
            "start_jd": j0,
            "end_jd": j1,
        }
        count = client.api("post", "filters/test/count", data=body)
        count.raise_for_status()
        n = count.json()["data"]["count"]

        body["limit"] = limit
        res = client.api("post", "filters/test", data=body)
        if not res.ok:
            raise RuntimeError(f"BOOM returned {res.status_code}: {res.text[:500]}")
        docs = res.json()["data"]["results"]
        if verbose:
            print(f"  JD {j0:.3f}-{j1:.3f}: {n} pass the filter in the MOC "
                  f"(depth {depth}), {len(docs)} returned", file=sys.stderr)
        if n > len(docs):
            print(f"  WARNING: truncated at --limit {limit}; raise it or shorten the window",
                  file=sys.stderr)
        results.extend(docs)

    # chunk edges are inclusive on both sides: drop repeated candids
    results = list({d.get("_id", id(d)): d for d in results}.values())

    # exact cone cut (the MOC overfills by up to one cell at the edge)
    centre = SkyCoord(ra * u.deg, dec * u.deg)
    kept = []
    for d in results:
        c = d.get("candidate", {})
        if c.get("ra") is None:
            kept.append(d)            # no coordinates projected: keep, cannot test
            continue
        sep = centre.separation(SkyCoord(c["ra"] * u.deg, c["dec"] * u.deg)).deg
        if sep <= radius_deg:
            d["sep_deg"] = float(sep)
            kept.append(d)
    return kept


def to_rows(docs: list[dict], jd0: float | None) -> list[dict]:
    rows = []
    for d in docs:
        c, p = d.get("candidate", {}), d.get("properties", {})
        rows.append({
            "objectId": d.get("objectId"),
            "candid": d.get("_id"),
            "ra": c.get("ra"),
            "dec": c.get("dec"),
            "sep_deg": round(d["sep_deg"], 5) if "sep_deg" in d else None,
            "jd": c.get("jd"),
            "dt_days": round(c["jd"] - jd0, 4) if jd0 and c.get("jd") else None,
            "band": ZTF_BANDS.get(c.get("fid"), c.get("fid")),
            "magpsf": c.get("magpsf"),
            "sigmapsf": c.get("sigmapsf"),
            "rb": c.get("rb"),
            "drb": c.get("drb"),
            "ndethist": c.get("ndethist"),
            "age_days": round(c["jd"] - c["jdstarthist"], 3)
                        if c.get("jd") and c.get("jdstarthist") else None,
            "star": p.get("star"),
            "near_brightstar": p.get("near_brightstar"),
            "Is_Quasar": d.get("Is_Quasar"),
        })
    return sorted(rows, key=lambda r: (r["objectId"] or "", r["jd"] or 0))


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("Usage")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--filter", required=True, help="filter pipeline JSON (e.g. ztf_basic.json)")
    g = p.add_argument_group("cone centre (a tile from tiles.csv, or explicit)")
    g.add_argument("--tiles", help="CSV written by tile_skymap.py")
    g.add_argument("--tile-id", type=int, default=0, help="default 0 = most probable tile")
    g.add_argument("--ra", type=float)
    g.add_argument("--dec", type=float)
    g.add_argument("--radius", type=float, help="deg; overrides the tile radius")
    p.add_argument("--t0", required=True, help="GW event time: ISO or JD")
    p.add_argument("--dt-days", type=float, default=7.0, help="days after t0 (default 7)")
    p.add_argument("--before-days", type=float, default=0.0, help="days before t0 (default 0)")
    p.add_argument("--survey", default="ZTF", help="BOOM survey (this filter is ZTF-specific)")
    p.add_argument("--programids", type=int, nargs="+", default=[1,2,3],
                   help="ZTF programids you have access to: 1=public 2=partnership 3=Caltech")
    p.add_argument("--limit", type=int, default=10_000, help="max alerts per 7-day request")
    p.add_argument("-o", "--out", default="filter_cone_alerts.csv")
    p.add_argument("--save-json", help="also save the full returned documents")
    p.add_argument("--outdir", default="outputs", help="directory for outputs (default outputs/)")
    a = p.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    a.out = os.path.join(a.outdir, a.out)
    if a.save_json:
        a.save_json = os.path.join(a.outdir, a.save_json)

    if a.tiles:
        ra, dec, radius = read_tile(a.tiles, a.tile_id)
    elif a.ra is not None and a.dec is not None:
        ra, dec, radius = a.ra, a.dec, None
    else:
        p.error("give --tiles, or --ra and --dec")
    radius = a.radius if a.radius is not None else radius
    if radius is None:
        p.error("--radius is required without --tiles")

    jd0 = to_jd(a.t0)
    jd_start, jd_end = jd0 - a.before_days, jd0 + a.dt_days
    pipeline = load_pipeline(a.filter)

    print(f"Filter {a.filter} on {a.survey} (programids {a.programids})\n"
          f"Cone RA={ra:.4f} Dec={dec:.4f} r={radius} deg, JD {jd_start:.4f}-{jd_end:.4f}",
          file=sys.stderr)

    load_credentials()
    client = BOOMClient()
    docs = run_filter_in_cone(client, pipeline, ra % 360.0, dec, radius,
                              jd_start, jd_end, a.programids, a.survey, a.limit)
    rows = to_rows(docs, jd0)

    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]) if rows else ["objectId"])
        w.writeheader()
        w.writerows(rows)
    if a.save_json:
        with open(a.save_json, "w") as f:
            json.dump(docs, f, indent=1, default=str)

    n_obj = len({r["objectId"] for r in rows})
    print(f"{len(rows)} alerts from {n_obj} objects inside the cone -> {a.out}",
          file=sys.stderr)


if __name__ == "__main__":
    main()