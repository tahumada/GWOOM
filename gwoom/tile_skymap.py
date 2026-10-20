#!/usr/bin/env python
"""
tile_skymap.py
==============

Cover the credible region of a HEALPix skymap (e.g. an LVK GW localization)
with circles ("tiles") of fixed angular radius, so that

  * every pixel of the region is covered (hard guarantee, verified at the end)
  * the number of circles, and therefore the overlap, is as small as we can
    reasonably get.

Circles cannot tile a surface without overlapping, so "minimal overlap" means
the thinnest *covering*. In the plane the optimum is a hexagonal lattice with
centre spacing sqrt(3)*r, which covers ~1.21x the area it needs to. This code:

  0. splits the region into connected islands (GW maps are often multimodal)
     and tiles each island separately:
  1. rotates the sky so the island's probability-weighted centroid sits on
     the equator (keeps the grid as undistorted as possible where it matters),
  2. lays down a hexagonal grid of centres in that frame,
  3. adds a circle on any pixel still uncovered (spherical-distortion patching),
  4. greedily removes circles whose pixels are all covered by other circles,
  5. repeats 1-4 for several random grid phases/rotations (and slightly
     tightened grids, which need less patching on the sphere) and keeps the
     solution with the fewest circles, then the least overlap,
  6. finally drops any circle made redundant by a neighbouring island.

"Covered" means the pixel centre lies within a circle. Use --strict to require
that the *whole* pixel lies inside a circle (radius is shrunk by the maximum
pixel radius when placing and checking circles).

Usage
-----
    python tile_skymap.py path/to/skymap.fits --radius 1.5
    python tile_skymap.py bayestar.multiorder.fits -r 2 --level 0.95 --plot
    python tile_skymap.py skymap.fits -r 1 --level 1.0 --out tiles.csv   # whole map, explicit name

As a module
-----------
    from tile_skymap import read_skymap, tile_skymap
    prob = read_skymap("skymap.fits")
    result = tile_skymap(prob, radius_deg=1.5, level=0.9)
    result.ra, result.dec          # circle centres, degrees
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from dataclasses import dataclass, field
from astropy.io import fits
import healpy as hp
import numpy as np
import logging
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

# project root = the folder above this script, so outputs land in the same
# place whatever directory the script is run from
DEFAULT_OUTDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "outputs")

SQRT3 = np.sqrt(3.0)

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="[%(levelname)s at %(asctime)s]: %(message)s",
                    datefmt="%H:%M:%S")


# --------------------------------------------------------------------------- #
# Reading
# --------------------------------------------------------------------------- #
def read_skymap(path: str, nside: int | None = None, radius_deg: float | None = None,
                verbose: bool = True) -> np.ndarray:
    """Read a HEALPix probability map and return it in RING ordering, normalised.

    Handles flat ("PROB" column) and LVK multi-order (NUNIQ) files when
    ligo.skymap is installed; falls back to healpy for flat maps otherwise.

    nside      : resample to this NSIDE (probability-conserving).
    radius_deg : if nside is None, pick the coarsest NSIDE whose pixels are
                 <= radius/4 across (never finer than the native map). This
                 keeps the tiling fast without changing the answer noticeably.
    """
    try:
        from ligo.skymap.io import read_sky_map
        logger.info(f"Reading skymap file {path} with ligo.skymap")
        prob, _meta = read_sky_map(path, moc=False, nest=False)
        logger.info(f"Done reading skymap file {path} with ligo.skymap")
    except ImportError:
        if verbose:
            logger.warning("ligo.skymap not found, reading with healpy "
                  "(multi-order/NUNIQ files will not work)", file=sys.stderr)
        prob = hp.read_map(path, field=0, nest=False, dtype=np.float64)

    prob = np.asarray(prob, dtype=np.float64)
    prob[~np.isfinite(prob) | (prob == hp.UNSEEN) | (prob < 0)] = 0.0
    native = hp.npix2nside(prob.size)
    logger.info(f"native NSIDE={native}")
    if nside is None and radius_deg is not None:
        nside = native
        while nside > 1 and np.degrees(hp.nside2resol(nside // 2)) <= radius_deg / 4:
            nside //= 2
    if nside is not None and nside != native:
        if nside > native:
            raise ValueError(f"Requested nside={nside} is finer than the map ({native}).")
        prob = hp.ud_grade(prob, nside, order_in="RING", order_out="RING", power=-2)
        if verbose:
            print(f"Resampled map from NSIDE={native} to NSIDE={nside}", file=sys.stderr)

    total = prob.sum()
    logger.info(f"total probability={total:.6g}")
    if total <= 0:
        raise ValueError("Map has no positive probability.")
    return prob / total


def event_name(path: str) -> str:
    """Event name from the OBJECT header"""
    # FIXME: ligo skymaps have the event name written in the format we want, whereas fermi skymaps have the object name in a different format. need to think about standardizing these.
    try:
        with fits.open(path) as hdul:
            for hdu in hdul:
                obj = str(hdu.header.get("OBJECT", "")).strip()
                if obj and obj.upper() not in ("FULLSKY", "PARTIAL"):   # healpy defaults
                    return obj.rstrip("/").split("/")[-1]      # LVK may give a URL
    except OSError:
        pass
    name = os.path.basename(path)
    for ext in (".gz", ".fits", ".fit", ".multiorder"):
        name = name.removesuffix(ext)
    return name


def credible_region(prob: np.ndarray, level: float) -> np.ndarray:
    """Boolean mask of the smallest set of pixels holding >= `level` probability.
    level >= 1 returns every pixel with non-zero probability."""
    if level >= 1.0:
        return prob > 0
    order = np.argsort(prob)[::-1]
    n = int(np.searchsorted(np.cumsum(prob[order]), level)) + 1
    mask = np.zeros(prob.size, dtype=bool)
    mask[order[:min(n, prob.size)]] = True
    return mask


def region_components(nside: int, region_pix: np.ndarray) -> list[np.ndarray]:
    """Split region pixels (RING) into connected islands. GW maps are often
    multimodal; each island gets its own grid frame, which keeps the hex grid
    undistorted on every island instead of only near the global centroid."""

    local = -np.ones(hp.nside2npix(nside), dtype=np.int64)
    local[region_pix] = np.arange(region_pix.size)
    nb = hp.get_all_neighbours(nside, region_pix)          # (8, n)
    rows = np.repeat(np.arange(region_pix.size)[None, :], 8, 0).ravel()
    cols = np.where(nb.ravel() >= 0, local[np.clip(nb.ravel(), 0, None)], -1)
    ok = cols >= 0
    g = coo_matrix((np.ones(ok.sum()), (rows[ok], cols[ok])),
                   shape=(region_pix.size, region_pix.size))
    n, labels = connected_components(g, directed=False)
    return [region_pix[labels == k] for k in range(n)]


# --------------------------------------------------------------------------- #
# Geometry helpers
# --------------------------------------------------------------------------- #
def _frame_matrix(centre_vec: np.ndarray, roll: float) -> np.ndarray:
    """Rotation R such that R @ centre_vec = (1, 0, 0), followed by a roll of
    `roll` radians about that axis. Rows are the new basis vectors."""
    e1 = centre_vec / np.linalg.norm(centre_vec)
    ref = np.array([0.0, 0.0, 1.0]) if abs(e1[2]) < 0.99 else np.array([1.0, 0.0, 0.0])
    e2 = np.cross(ref, e1)
    e2 /= np.linalg.norm(e2)
    e3 = np.cross(e1, e2)
    R0 = np.vstack([e1, e2, e3])
    c, s = np.cos(roll), np.sin(roll)
    Rx = np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
    return Rx @ R0


def _lonlat(vec: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    lon = np.degrees(np.arctan2(vec[..., 1], vec[..., 0]))
    lat = np.degrees(np.arcsin(np.clip(vec[..., 2], -1, 1)))
    return lon, lat


def _vec(lon_deg, lat_deg) -> np.ndarray:
    lon, lat = np.radians(lon_deg), np.radians(lat_deg)
    return np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)], axis=-1)


def _hex_centres(lon, lat, r, phase_lat, phase_lon) -> np.ndarray:
    """Hexagonal grid of centres (rotated-frame lon/lat, degrees) covering the
    bounding box of the region pixels, with a margin of one radius."""
    lat_min, lat_max = max(lat.min() - r, -90.0), min(lat.max() + r, 90.0)
    full_ring = (lon.max() - lon.min()) > 270.0
    lon_min, lon_max = lon.min(), lon.max()

    dlat = 1.5 * r
    rows = np.arange(lat_min - dlat * (1 - phase_lat), lat_max + dlat, dlat)
    out = []
    for i, phi in enumerate(rows):
        phi = float(np.clip(phi, -90.0, 90.0))
        # Use the latitude in this row's band closest to the equator, where the
        # true spacing between neighbours is largest -> no gaps between rows.
        phi_eff = max(abs(phi) - r, 0.0)
        if phi_eff >= 89.999:
            out.append((0.0, phi))
            continue
        dlon = SQRT3 * r / np.cos(np.radians(phi_eff))
        # Every row is anchored on the region's central meridian (lon=0 in this
        # frame) so neighbouring rows stay properly staggered where it matters.
        offset = 0.5 * (i % 2) + phase_lon
        margin = r / max(np.cos(np.radians(min(abs(phi) + r, 89.9))), 1e-3)
        if full_ring or dlon >= 360.0 or (lon_max - lon_min + 2 * margin) >= 360.0:
            n = max(1, int(np.ceil(360.0 / dlon)))
            lons = (np.arange(n) + offset) * 360.0 / n
        else:
            k0 = int(np.floor((lon_min - margin) / dlon - offset)) - 1
            k1 = int(np.ceil((lon_max + margin) / dlon - offset)) + 1
            lons = (np.arange(k0, k1 + 1) + offset) * dlon
        out.extend((lo, phi) for lo in lons)
    out = np.asarray(out)
    out[:, 0] = (out[:, 0] + 180.0) % 360.0 - 180.0
    return out


# --------------------------------------------------------------------------- #
# Core
# --------------------------------------------------------------------------- #
@dataclass
class Tiling:
    ra: np.ndarray                 # circle centres, deg
    dec: np.ndarray
    radius_deg: float              # requested radius
    level: float
    nside: int
    region_mask: np.ndarray        # credible-region pixels (RING)
    coverage: np.ndarray           # number of circles covering each pixel (RING)
    prob_in_circle: np.ndarray     # map probability inside each circle (overlaps double-count)
    prob_assigned: np.ndarray      # region probability assigned to nearest circle (sums w/o double counting)
    summary: dict = field(default_factory=dict)


def _discs(nside, vecs, r_rad, region_idx_of):
    """For each centre, the region-local indices of pixels whose centres are inside."""
    discs = []
    for v in vecs:
        pix = hp.query_disc(nside, v, r_rad, inclusive=False, nest=False)
        loc = region_idx_of[pix]
        discs.append(loc[loc >= 0])
    return discs


def _solve(nside, region_pix, region_vec, weights, r_deg, rng, trial, shrink=1.0):
    centre = (region_vec * weights[:, None]).sum(0)
    if np.linalg.norm(centre) < 1e-6:          # (near) all-sky: any frame is fine
        centre = np.array([1.0, 0.0, 0.0])
    if trial == 0:
        roll, ph_lat, ph_lon = 0.0, 0.0, 0.0
    else:
        roll, ph_lat, ph_lon = rng.uniform(0, np.pi / 3), rng.uniform(), rng.uniform()
    R = _frame_matrix(centre, roll)

    lon, lat = _lonlat(region_vec @ R.T)
    grid = _hex_centres(lon, lat, r_deg * shrink, ph_lat, ph_lon)
    vecs = list(_vec(grid[:, 0], grid[:, 1]) @ R)   # back to equatorial frame

    npix = hp.nside2npix(nside)
    region_idx_of = -np.ones(npix, dtype=np.int64)
    region_idx_of[region_pix] = np.arange(region_pix.size)
    r_rad = np.radians(r_deg)

    discs = _discs(nside, vecs, r_rad, region_idx_of)
    keep = [i for i, d in enumerate(discs) if d.size]
    vecs = [vecs[i] for i in keep]
    discs = [discs[i] for i in keep]

    counts = np.zeros(region_pix.size, dtype=np.int64)
    for d in discs:
        np.add.at(counts, d, 1)

    # Patch any holes left by spherical distortion of the grid.
    logger.debug(f"{np.count_nonzero(counts == 0)} holes after grid placement")
    while True:
        holes = np.flatnonzero(counts == 0)
        if holes.size == 0:
            break
        # Centre the patch on the local cluster of holes around the first hole
        # (one circle then usually fixes the whole cluster).
        hv = region_vec[holes]
        near = hv[hv @ hv[0] >= np.cos(r_rad)]
        v = near.mean(0)
        v /= np.linalg.norm(v)
        if np.min(near @ v) < np.cos(r_rad):
            v = hv[0]
        d = _discs(nside, [v], r_rad, region_idx_of)[0]
        vecs.append(v)
        discs.append(d)
        np.add.at(counts, d, 1)

    # Greedy pruning: drop circles that add nothing, edge-clippers first.
    alive = np.ones(len(discs), dtype=bool)
    for _ in range(2):
        unique = np.array([np.count_nonzero(counts[d] == 1) for d in discs])
        size = np.array([d.size for d in discs])
        for i in np.lexsort((size, unique)):
            if alive[i] and np.all(counts[discs[i]] > 1):
                counts[discs[i]] -= 1
                alive[i] = False

    vecs = np.asarray(vecs)[alive]
    overlap = int(counts.sum() - counts.size)
    return vecs, counts, overlap


def tile_skymap(prob: np.ndarray, radius_deg: float, level: float = 0.9,
                strict: bool = False, n_trials: int = 10, seed: int = 0,
                verbose: bool = True) -> Tiling:
    """Cover the `level` credible region of `prob` (RING) with circles of
    radius `radius_deg`. See module docstring for the method."""
    nside = hp.npix2nside(prob.size)
    mask = credible_region(prob, level)
    region_pix = np.flatnonzero(mask)
    region_vec = np.array(hp.pix2vec(nside, region_pix)).T

    r_place = radius_deg
    if strict:
        r_place = radius_deg - np.degrees(hp.max_pixrad(nside))
        if r_place <= 0:
            raise ValueError("Radius smaller than a pixel: use a finer --nside or drop --strict.")

    rng = np.random.default_rng(seed)
    islands = region_components(nside, region_pix)
    if verbose:
        logger.info(f"Region: {region_pix.size} pixels in {len(islands)} island(s)")

    # Spherical distortion breaks the perfect hex lattice away from the frame
    # equator. Slightly tighter grids leave fewer holes to patch; which one wins
    # depends on the region, so try a few and keep the best -- per island.
    all_vecs = []
    logger.info(f"tiling with {n_trials} random grid phases and shrink factors")
    for k, isl in enumerate(islands):
        isl_vec = np.array(hp.pix2vec(nside, isl)).T
        best = None
        for t in range(max(1, n_trials)):
            state = rng.bit_generator.state
            for shrink in (1.0, 0.97, 0.94, 0.9):
                rng.bit_generator.state = state      # same phase for every shrink
                vecs, counts, overlap = _solve(nside, isl, isl_vec, prob[isl],
                                               r_place, rng, t, shrink)
                key = (len(vecs), overlap)
                if best is None or key < best[0]:
                    best = (key, vecs)
            rng.random()                             # advance to a new phase
        if verbose:
            print(f"  island {k}: {isl.size} pix -> {best[0][0]} circles", file=sys.stderr)
        all_vecs.append(best[1])

    # Global clean-up: islands closer than 2r can share circles.
    vecs = np.vstack(all_vecs)
    npix = hp.nside2npix(nside)
    region_idx_of = -np.ones(npix, dtype=np.int64)
    region_idx_of[region_pix] = np.arange(region_pix.size)
    discs = _discs(nside, list(vecs), np.radians(r_place), region_idx_of)
    counts = np.zeros(region_pix.size, dtype=np.int64)
    for d in discs:
        np.add.at(counts, d, 1)
    alive = np.ones(len(discs), dtype=bool)
    unique = np.array([np.count_nonzero(counts[d] == 1) for d in discs])
    for i in np.lexsort((np.array([d.size for d in discs]), unique)):
        if np.all(counts[discs[i]] > 1):
            counts[discs[i]] -= 1
            alive[i] = False
    vecs = vecs[alive]

    assert np.all(counts >= 1), "Coverage check failed"  # should never trigger

    ra, dec = hp.vec2ang(vecs, lonlat=True)
    ra = ra % 360.0

    coverage = np.zeros(prob.size, dtype=np.int64)
    coverage[region_pix] = counts

    # Probability per circle (whole map, overlaps counted in every circle) ...
    r_rad = np.radians(r_place)
    prob_in = np.array([prob[hp.query_disc(nside, v, r_rad, inclusive=False)].sum()
                        for v in vecs])
    # ... and a true partition of the region: each pixel goes to its nearest centre.
    nearest = np.empty(region_pix.size, dtype=np.int64)
    for s in range(0, region_pix.size, 200_000):
        nearest[s:s + 200_000] = np.argmax(region_vec[s:s + 200_000] @ vecs.T, axis=1)
    prob_assigned = np.bincount(nearest, weights=prob[region_pix], minlength=len(vecs))

    pix_area = hp.nside2pixarea(nside, degrees=True)
    cap_area = 2 * np.pi * (1 - np.cos(np.radians(radius_deg))) * (180 / np.pi) ** 2
    region_area = region_pix.size * pix_area
    summary = dict(
        n_tiles=len(vecs),
        radius_deg=radius_deg,
        level=level,
        nside=nside,
        strict=strict,
        region_area_deg2=region_area,
        region_prob=float(prob[region_pix].sum()),
        tiles_area_deg2=len(vecs) * cap_area,
        covering_density=len(vecs) * cap_area / region_area,   # 1.21 = ideal hex lattice
        frac_region_pix_multiply_covered=float(np.mean(counts > 1)),
        mean_circles_per_pixel=float(counts.mean()),
        uncovered_pixels=int(np.count_nonzero(counts == 0)),
    )
    return Tiling(ra, dec, radius_deg, level, nside, mask, coverage,
                  prob_in, prob_assigned, summary)


# --------------------------------------------------------------------------- #
# Output
# --------------------------------------------------------------------------- #
def write_csv(t: Tiling, path: str) -> None:
    order = np.argsort(-t.prob_assigned)          # most probable tiles first
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["tile_id", "ra_deg", "dec_deg", "radius_deg",
                    "prob_in_circle", "prob_assigned"])
        for k, i in enumerate(order):
            w.writerow([k, f"{t.ra[i]:.6f}", f"{t.dec[i]:.6f}", t.radius_deg,
                        f"{t.prob_in_circle[i]:.6e}", f"{t.prob_assigned[i]:.6e}"])


def plot(t: Tiling, path: str) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    m = np.where(t.region_mask, t.coverage, hp.UNSEEN).astype(float)
    vecs = np.array(hp.ang2vec(t.ra, t.dec, lonlat=True))
    c = vecs.mean(0)
    c = c / np.linalg.norm(c) if np.linalg.norm(c) > 1e-6 else np.array([1.0, 0, 0])
    region_vec = np.array(hp.pix2vec(t.nside, np.flatnonzero(t.region_mask))).T
    extent = np.degrees(np.arccos(np.clip(region_vec @ c, -1, 1)).max()) + t.radius_deg

    title = (f"{t.summary['n_tiles']} circles, r={t.radius_deg}°, "
             f"{int(t.level * 100)}% region, density {t.summary['covering_density']:.2f}")
    vmax = max(2, int(t.coverage.max()))
    if extent < 60:
        rot_ra, rot_dec = hp.vec2ang(c, lonlat=True)
        hp.gnomview(m, rot=(float(rot_ra[0]), float(rot_dec[0])), xsize=900,
                    reso=2.2 * extent * 60 / 900, title=title, unit="circles per pixel",
                    min=0, max=vmax, cmap="viridis", notext=True)
    else:
        hp.mollview(m, title=title, unit="circles per pixel", min=0, max=vmax, cmap="viridis")
    hp.graticule(dpar=15, dmer=15, alpha=0.3)

    if len(vecs) <= 3000:                      # circle outlines
        th = np.linspace(0, 2 * np.pi, 72)
        r = np.radians(t.radius_deg)
        for v in vecs:
            R = _frame_matrix(v, 0.0)
            ring = np.stack([np.full_like(th, np.cos(r)), np.sin(r) * np.cos(th),
                             np.sin(r) * np.sin(th)], -1) @ R
            lo, la = hp.vec2ang(ring, lonlat=True)
            hp.projplot(lo, la, lonlat=True, color="w", lw=0.5, alpha=0.8)
    hp.projscatter(t.ra, t.dec, lonlat=True, s=4, c="r")
    plt.savefig(path, dpi=150, bbox_inches="tight")
    logger.info(f"Wrote tiling plot to {path}")
    plt.close("all")


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("Usage")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("skymap", help="HEALPix .fits / .fits.gz (flat or LVK multi-order)")
    p.add_argument("-r", "--radius", type=float, required=True, help="circle radius, deg")
    p.add_argument("-l", "--level", type=float, default=0.9,
                   help="credible level to cover (default 0.9; 1.0 = all non-zero pixels)")
    p.add_argument("--nside", type=int, default=None,
                   help="work at this NSIDE (default: auto, pixels <= radius/4)")
    p.add_argument("--native", action="store_true", help="use the map's native NSIDE")
    p.add_argument("--strict", action="store_true",
                   help="require whole pixels (not just centres) inside a circle")
    p.add_argument("--trials", type=int, default=10, help="grid phases to try (default 10)")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("-o", "--out", default=None,
                   help="output CSV (default <event>_tiles.csv, event from the OBJECT header)")
    p.add_argument("--plot", nargs="?", const="", default=None,
                   help="save a PNG of the tiling (default name <event>_tiles.png)")
    p.add_argument("--outdir", default=DEFAULT_OUTDIR,
                   help="output directory, created if missing (default: <project root>/outputs)")
    a = p.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    event = event_name(a.skymap)
    a.out = os.path.join(a.outdir, a.out or f"{event}_tiles.csv")
    if a.plot is not None:
        a.plot = os.path.join(a.outdir, a.plot or f"{event}_tiles.png")

    prob = read_skymap(a.skymap, nside=a.nside,
                       radius_deg=None if a.native else a.radius)
    t = tile_skymap(prob, a.radius, a.level, strict=a.strict, n_trials=a.trials, seed=a.seed)
    write_csv(t, a.out)
    if a.plot:
        plot(t, a.plot)

    rows = [(k, f"{v:.4g}" if isinstance(v, float) else str(v)) for k, v in t.summary.items()]
    wk = max(len(k) for k, _ in rows + [("quantity", "")])
    wv = max(len(v) for _, v in rows + [("", "value")])
    print(f"\n┌─{'─' * wk}─┬─{'─' * wv}─┐")
    print(f"│ {'quantity':<{wk}} │ {'value':>{wv}} │")
    print(f"├─{'─' * wk}─┼─{'─' * wv}─┤")
    for k, v in rows:
        print(f"│ {k:<{wk}} │ {v:>{wv}} │")
    print(f"└─{'─' * wk}─┴─{'─' * wv}─┘")
    logger.info(f"Wrote {a.out}" + (f" and {a.plot}" if a.plot else ""))
    logger.info(f"{t.summary['n_tiles']} circles, radius {t.radius_deg} deg")

if __name__ == "__main__":
    main()
