"""
Cone-search ZTF / LSST (DECam later) over the tiles of a skymap, keep sources
first detected no earlier than 1 hr before the trigger, merge, write a text file.

    python gwoom/multisurvey_query.py skymaps/GRB261013a.fits -r 5 --level 0.9
"""

import argparse
import copy
import os
import sys

import astropy.units as u
import numpy as np
import requests
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from dotenv import load_dotenv
from mocpy import MOC

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from filters.filters import all_filters
from gwoom.tile_skymap import read_skymap, tile_skymap

# project root = the folder above this script, so outputs land in the same
# place whatever directory the script is run from
DEFAULT_OUTDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "outputs")

BOOM_URL = "https://api.kaboom.caltech.edu"
SURVEYS = ("ztf", "lsst", "decam")
PRE_TRIGGER_DAYS = 1 / 24     # first detection may be at most 1 hr before the trigger
MAX_WINDOW_DAYS = 7           # BOOM's /filters/test limit per request
MATCH_ARCSEC = 1.5            # ZTF/LSST sources closer than this are the same transient
# fields every filter's final $project must keep so we can report the source
KEEP_FIELDS = ["candidate.ra", "candidate.dec", "candidate.magpsf", "candidate.band"]


def auth():
    # login to boom and get the token
    load_dotenv()
    user, password = os.getenv("BOOM_USERNAME"), os.getenv("BOOM_PASSWORD")
    if not user or not password:
        sys.exit("Set BOOM_USERNAME and BOOM_PASSWORD (environment or .env)")
    r = requests.post(f"{BOOM_URL}/auth", data={"username": user, "password": password})
    r.raise_for_status()
    return r.json()["access_token"]


def get_trigger_time(skymap_path):
    # trigger time from the skymap header (LVK and GBM maps both have DATE-OBS)
    with fits.open(skymap_path) as hdul:
        for hdu in hdul:
            if hdu.header.get("DATE-OBS"):
                return Time(hdu.header["DATE-OBS"]).jd
    sys.exit(f"No DATE-OBS in {skymap_path}: pass --tstart")


def load_cones(skymap_path, radius, level):
    # tile the skymap directly with tile_skymap; returns [(ra, dec, radius), ...]
    tiles = tile_skymap(read_skymap(skymap_path, radius_deg=radius), radius, level)
    return [(ra, dec, radius) for ra, dec in zip(tiles.ra, tiles.dec)]


def pick_filters(surveys=SURVEYS):
    # pick one filter for each survey from the list of filters available. ask as input so its forced step.
    names = list(all_filters)
    print("Available filters:")
    for i, name in enumerate(names):
        print(f"  [{i}] {name}")
    picked = {survey: None for survey in SURVEYS}
    for survey in surveys:
        choice = input(f"Filter number for {survey.upper()} (blank to skip): ").strip()
        picked[survey] = all_filters[names[int(choice)]] if choice else None
    return picked


def add_tstart_to_filters(filters, tstart):
    # for each filter in the picked filters, add the tstart cut to the filter. and call them as ztf_filters, lsst_filters, decam_filters. return them as a dict.
    cut = tstart - PRE_TRIGGER_DAYS
    out = {}
    for survey, pipeline in filters.items():
        if pipeline is None:
            out[f"{survey}_filters"] = None
            continue
        pipeline = copy.deepcopy(pipeline)
        for field in KEEP_FIELDS:
            pipeline[-1]["$project"].setdefault(field, 1)
        if survey == "ztf":
            # jdstarthist = first detection at this position
            pipeline.insert(0, {"$match": {"candidate.jdstarthist": {"$gte": cut}}})
        elif survey == "lsst":
            # jd_min_prv (computed by the LSST filters) = first previous detection;
            # no previous detections -> this alert is the first
            pipeline.append({"$match": {"$expr": {"$gte": [
                {"$ifNull": ["$jd_min_prv", "$candidate.jd"]}, cut]}}})
        else:
            print(f"No tstart cut defined for {survey}, using filter as is")
        out[f"{survey}_filters"] = pipeline
    return out


def _query_boom(survey, cones, pipeline, tstart, tend, token, permissions):
    # run the filter in every cone, in <=7 day chunks; one row per alert
    headers = {"Authorization": f"Bearer {token}"}
    alerts = {}
    for ra, dec, radius in cones:
        moc = MOC.from_cone(lon=ra * u.deg, lat=dec * u.deg, radius=radius * u.deg, max_depth=10)
        for jd0 in np.arange(tstart - PRE_TRIGGER_DAYS, tend, MAX_WINDOW_DAYS):
            body = {
                "survey": survey,
                "pipeline": pipeline,
                "permissions": permissions,
                "moc_ascii": moc.to_string(format="ascii"),
                "start_jd": jd0,
                "end_jd": min(jd0 + MAX_WINDOW_DAYS, tend),
            }
            r = requests.post(f"{BOOM_URL}/filters/test", json=body, headers=headers)
            if not r.ok:
                print(f"{survey} cone ({ra:.2f}, {dec:.2f}) failed: {r.status_code} {r.text[:200]}")
                continue
            for a in r.json()["data"]["results"]:
                alerts[a["_id"]] = a          # tiles overlap: same alert only once
        print(f"{survey} cone ({ra:.2f}, {dec:.2f}): {len(alerts)} alerts so far")
    return list(alerts.values())


def query_ztf(cones, ztf_filters, tstart, tend, token):
    # use the cones, and the filters to query the ztf database.
    return _query_boom("ZTF", cones, ztf_filters, tstart, tend, token, {"ZTF": [1, 2, 3]})


def query_lsst(cones, lsst_filters, tstart, tend, token):
    # use the cones, and the filters to query the lsst database.
    return _query_boom("LSST", cones, lsst_filters, tstart, tend, token, {})


def query_decam(cones, decam_filters, tstart, tend, token):
    # use the cones, and the filters to query the decam database.
    print("DECam not available yet, skipping")
    return []


def merge_candidates(ztf_candidates, lsst_candidates, decam_candidates):
    # merge the candidates from the three surveys into a single list of candidates.
    # FIXME: need checks
    objects = {}
    for survey, alerts in (("ZTF", ztf_candidates), ("LSST", lsst_candidates),
                           ("DECAM", decam_candidates)):
        for a in alerts:
            key = (survey, str(a["objectId"]))
            if key not in objects:
                objects[key] = {"survey": survey, "objectId": str(a["objectId"]),
                                "ra": a["candidate"]["ra"], "dec": a["candidate"]["dec"]}

    # sources from different surveys within MATCH_ARCSEC are the same transient
    merged = []
    for obj in objects.values():
        pos = SkyCoord(obj["ra"] * u.deg, obj["dec"] * u.deg)
        for t in merged:
            if (obj["survey"] not in t["surveys"]
                    and pos.separation(SkyCoord(t["ra"] * u.deg, t["dec"] * u.deg)).arcsec < MATCH_ARCSEC):
                t["surveys"][obj["survey"]] = obj["objectId"]
                break
        else:
            merged.append({"ra": obj["ra"], "dec": obj["dec"],
                           "surveys": {obj["survey"]: obj["objectId"]}})
    return merged


def write_candidates(candidates, path):
    with open(path, "w") as f:
        f.write("# ra_deg dec_deg ids\n")
        for c in candidates:
            ids = " ".join(f"{s}:{i}" for s, i in c["surveys"].items())
            f.write(f"{c['ra']:.6f} {c['dec']:+.6f} {ids}\n")
    print(f"{len(candidates)} candidates -> {path}")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("skymap")
    p.add_argument("-r", "--radius", type=float, required=True, help="cone radius, deg")
    p.add_argument("--level", type=float, default=0.9, help="credible level to tile")
    p.add_argument("--tstart", help="trigger time (ISO); default: DATE-OBS of the skymap")
    p.add_argument("--dt-days", type=float, default=7, help="search up to this long after trigger")
    p.add_argument("--surveys", nargs="+", choices=SURVEYS, default=list(SURVEYS),
                   help="surveys to query, e.g. --surveys ztf lsst (default: all)")
    p.add_argument("-o", "--out", default="candidates.txt")
    p.add_argument("--outdir", default=DEFAULT_OUTDIR,
                   help="output directory, created if missing (default: <project root>/outputs)")
    args = p.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    args.out = os.path.join(args.outdir, args.out)

    tstart = Time(args.tstart).jd if args.tstart else get_trigger_time(args.skymap)
    tend = min(tstart + args.dt_days, Time.now().jd)

    cones = load_cones(args.skymap, args.radius, args.level)
    filters = add_tstart_to_filters(pick_filters(args.surveys), tstart)
    token = auth()

    ztf = query_ztf(cones, filters["ztf_filters"], tstart, tend, token) if filters["ztf_filters"] else []
    lsst = query_lsst(cones, filters["lsst_filters"], tstart, tend, token) if filters["lsst_filters"] else []
    decam = query_decam(cones, filters["decam_filters"], tstart, tend, token) if filters["decam_filters"] else []

    write_candidates(merge_candidates(ztf, lsst, decam), args.out)
