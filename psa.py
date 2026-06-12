#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "astrometry>=4.1.2,<5",
#     "astropy>=5.3",
#     "numpy>=1.24",
#     "sep>=1.2.1",
#     "pillow>=10.1",
# ]
# ///
"""psa.py -- portable plate-solve + annotate.

Self-contained replacement for the classic local pipeline of a shell
wrapper around solve-field + plot-constellations. Solving uses the
astrometry.net engine compiled into the `astrometry` wheel; index files
and annotation catalogs are downloaded into a local cache on first use
and everything runs offline afterwards.

Usage:
    psa.py image.fit                  solve + annotate
    psa.py image.jpg --hd             also label Henry Draper stars
    psa.py --prefetch --hd            warm the cache for offline use
    psa.py --check                    verify cached catalogs

Outputs (in "<name> Solved/"):
    <name>.wcs        FITS WCS header of the solution
    annotations.png   annotated image (--transparent for overlay only)
    solution.json     machine-readable solve summary
"""

from __future__ import annotations

import argparse
import json
import math
import multiprocessing
import re
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np

try:
    import astrometry
except ImportError:
    astrometry = None

if astrometry is not None and not hasattr(astrometry, "Solver"):
    sys.exit(
        "error: the 'astrometry' module that was imported is the astrometry.net\n"
        "utility package (system site-packages), not the PyPI solver. Run this\n"
        "script in an isolated environment:  uv run psa.py ...  (or ./psa.py ...)"
    )

if astrometry is None:
    sys.exit(
        "error: dependencies missing. Run via uv (recommended):\n"
        "    uv run psa.py ...     or just  ./psa.py ...\n"
        "or install manually:\n"
        "    pip install 'astrometry>=4.1.2,<5' astropy numpy sep pillow"
    )

import sep
from astropy.io import fits
from astropy.wcs import WCS
from PIL import Image, ImageDraw, ImageFont

Image.MAX_IMAGE_PIXELS = None  # astro mosaics exceed PIL's decompression-bomb default

VERSION = "1.1.0"

DEFAULT_CACHE = Path.home() / ".cache" / "psa"

URLS = {
    "constellation-lines.json": "https://raw.githubusercontent.com/ofrohn/d3-celestial/master/data/constellations.lines.json",
    "iau-csn.txt": "https://www.pas.rochester.edu/~emamajek/WGSN/IAU-CSN.txt",
    "openngc.csv": "https://raw.githubusercontent.com/mattiaverga/OpenNGC/master/database_files/NGC.csv",
    "hd.fits": "http://data.astrometry.net/hd.fits",
    "stars.6.json": "https://raw.githubusercontent.com/ofrohn/d3-celestial/master/data/stars.6.json",
    "starnames.json": "https://raw.githubusercontent.com/ofrohn/d3-celestial/master/data/starnames.json",
}

# Index series shipped by data.astrometry.net via the astrometry package.
# Values are the valid scale range. Default mirrors a typical wide-field rig
# (quad sizes ~1.4 - 33 deg; suits fields from a few degrees to all-sky).
SERIES_SCALES = {
    "4100": (7, 19),
    "4200": (0, 19),
    "5000": (0, 7),
    "5200": (0, 6),
    "5200_heavy": (0, 6),
}

CONSTELLATION_NAMES = {
    "And": "Andromeda", "Ant": "Antlia", "Aps": "Apus", "Aqr": "Aquarius",
    "Aql": "Aquila", "Ara": "Ara", "Ari": "Aries", "Aur": "Auriga",
    "Boo": "Bootes", "Cae": "Caelum", "Cam": "Camelopardalis", "Cnc": "Cancer",
    "CVn": "Canes Venatici", "CMa": "Canis Major", "CMi": "Canis Minor",
    "Cap": "Capricornus", "Car": "Carina", "Cas": "Cassiopeia",
    "Cen": "Centaurus", "Cep": "Cepheus", "Cet": "Cetus", "Cha": "Chamaeleon",
    "Cir": "Circinus", "Col": "Columba", "Com": "Coma Berenices",
    "CrA": "Corona Australis", "CrB": "Corona Borealis", "Crv": "Corvus",
    "Crt": "Crater", "Cru": "Crux", "Cyg": "Cygnus", "Del": "Delphinus",
    "Dor": "Dorado", "Dra": "Draco", "Equ": "Equuleus", "Eri": "Eridanus",
    "For": "Fornax", "Gem": "Gemini", "Gru": "Grus", "Her": "Hercules",
    "Hor": "Horologium", "Hya": "Hydra", "Hyi": "Hydrus", "Ind": "Indus",
    "Lac": "Lacerta", "Leo": "Leo", "LMi": "Leo Minor", "Lep": "Lepus",
    "Lib": "Libra", "Lup": "Lupus", "Lyn": "Lynx", "Lyr": "Lyra",
    "Men": "Mensa", "Mic": "Microscopium", "Mon": "Monoceros", "Mus": "Musca",
    "Nor": "Norma", "Oct": "Octans", "Oph": "Ophiuchus", "Ori": "Orion",
    "Pav": "Pavo", "Peg": "Pegasus", "Per": "Perseus", "Phe": "Phoenix",
    "Pic": "Pictor", "Psc": "Pisces", "PsA": "Piscis Austrinus",
    "Pup": "Puppis", "Pyx": "Pyxis", "Ret": "Reticulum", "Sge": "Sagitta",
    "Sgr": "Sagittarius", "Sco": "Scorpius", "Scl": "Sculptor",
    "Sct": "Scutum", "Ser": "Serpens", "Ser1": "Serpens Caput",
    "Ser2": "Serpens Cauda", "Sex": "Sextans", "Tau": "Taurus",
    "Tel": "Telescopium", "Tri": "Triangulum", "TrA": "Triangulum Australe",
    "Tuc": "Tucana", "UMa": "Ursa Major", "UMi": "Ursa Minor", "Vel": "Vela",
    "Vir": "Virgo", "Vol": "Volans", "Vul": "Vulpecula",
}

# OpenNGC object types worth annotating (skips stars, duplicates, non-existent).
NGC_KEEP_TYPES = {
    "G", "GPair", "GTrpl", "GGroup", "GCl", "OCl", "Cl+N", "PN", "SNR",
    "Neb", "EmN", "RfN", "HII", "DrkN", "Nova",
}

COLORS = {
    "lines": (255, 210, 130),
    "conname": (255, 210, 130),
    "star": (255, 255, 255),
    "ngc": (140, 255, 140),
    "hd": (170, 200, 255),
}


def log(msg: str) -> None:
    print(msg, flush=True)


# ----------------------------------------------------------------------------
# Cache + downloads
# ----------------------------------------------------------------------------

def fetch(url: str, dest: Path, label: str) -> bool:
    """Download url -> dest unless dest already exists. Returns availability."""
    if dest.exists() and dest.stat().st_size > 0:
        return True
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    log(f"  downloading {label} ...")
    try:
        req = urllib.request.Request(url, headers={"User-Agent": f"psa.py/{VERSION}"})
        with urllib.request.urlopen(req, timeout=120) as r, open(tmp, "wb") as f:
            while chunk := r.read(1 << 18):
                f.write(chunk)
        tmp.rename(dest)
        return True
    except Exception as e:
        tmp.unlink(missing_ok=True)
        log(f"  warning: could not download {label}: {e}")
        return False


def parse_scales(spec: str) -> set[int]:
    if "-" in spec:
        lo, hi = spec.split("-", 1)
        return set(range(int(lo), int(hi) + 1))
    return {int(s) for s in spec.split(",")}


def get_index_files(cache: Path, series_csv: str, scales: set[int]) -> list[Path]:
    """Download (if needed) and return index file paths for the chosen series."""
    index_dir = cache / "index"
    index_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for name in [s.strip() for s in series_csv.split(",") if s.strip()]:
        series = getattr(astrometry, f"series_{name}", None)
        if series is None:
            sys.exit(f"error: unknown index series '{name}' "
                     f"(known: {', '.join(SERIES_SCALES)})")
        lo, hi = SERIES_SCALES.get(name, (min(scales), max(scales)))
        use = {s for s in scales if lo <= s <= hi}
        if not use:
            log(f"  note: series {name} has no scales in requested range, skipping")
            continue
        paths += series.index_files(cache_directory=index_dir, scales=use)
    if not paths:
        sys.exit("error: no index files selected")
    return paths


# ----------------------------------------------------------------------------
# Catalog loading (raw download -> compiled .npz, offline afterwards)
# ----------------------------------------------------------------------------

def _catdir(cache: Path) -> Path:
    d = cache / "catalogs"
    d.mkdir(parents=True, exist_ok=True)
    return d


def load_constellation_lines(cache: Path):
    """Return (segments[N,4] ra1,dec1,ra2,dec2 deg, seg_con[N] idx, con_ids)."""
    cat = _catdir(cache)
    npz = cat / "conlines.npz"
    if not npz.exists():
        raw = cat / "constellation-lines.json"
        if not fetch(URLS["constellation-lines.json"], raw, "constellation lines"):
            return None
        data = json.loads(raw.read_text())
        segs, seg_con, ids = [], [], []
        for feat in data["features"]:
            cid = feat["id"]
            if cid not in ids:
                ids.append(cid)
            ci = ids.index(cid)
            for line in feat["geometry"]["coordinates"]:
                for (lo1, la1), (lo2, la2) in zip(line, line[1:]):
                    segs.append((lo1 % 360.0, la1, lo2 % 360.0, la2))
                    seg_con.append(ci)
        np.savez_compressed(npz, segs=np.array(segs, "f8"),
                            seg_con=np.array(seg_con, "i2"),
                            ids=np.array(ids, "U8"))
    d = np.load(npz)
    return d["segs"], d["seg_con"], list(d["ids"])


def load_bright_stars(cache: Path):
    """IAU named stars -> (name, bayer-con designation, ra, dec, hip)."""
    cat = _catdir(cache)
    npz = cat / "bright-stars-v2.npz"
    if not npz.exists():
        raw = cat / "iau-csn.txt"
        if not fetch(URLS["iau-csn.txt"], raw, "IAU star names"):
            return None
        names, desigs, ras, decs, hips = [], [], [], [], []
        date_re = re.compile(r"^\d{4}-\d{2}-\d{2}$")
        for line in raw.read_text(encoding="utf-8", errors="replace").splitlines():
            if not line.strip() or line.lstrip()[0] in "#$":
                continue
            toks = line.split()
            di = next((i for i, t in enumerate(toks) if date_re.match(t)), None)
            if di is None or di < 10:
                continue
            try:
                ra, dec = float(toks[di - 2]), float(toks[di - 1])
            except ValueError:
                continue
            name = line[:18].strip()
            greek, con = toks[di - 10], toks[di - 9]
            desig = f"{greek} {con}" if greek != "_" and con != "_" else ""
            hip = int(toks[di - 4]) if toks[di - 4].isdigit() else 0
            names.append(name); desigs.append(desig)
            ras.append(ra); decs.append(dec); hips.append(hip)
        if len(names) < 100:
            log(f"  warning: IAU-CSN parse looks wrong ({len(names)} stars)")
        np.savez_compressed(npz, name=np.array(names, "U24"),
                            desig=np.array(desigs, "U16"),
                            ra=np.array(ras, "f8"), dec=np.array(decs, "f8"),
                            hip=np.array(hips, "u4"))
    d = np.load(npz)
    return d["name"], d["desig"], d["ra"], d["dec"], d["hip"]


def load_bayer_stars(cache: Path):
    """Bright stars with Bayer/Flamsteed designations but no IAU proper name.

    Built from d3-celestial stars.6.json (positions, magnitudes, HIP ids)
    joined with starnames.json (designations). Restores plot-constellations'
    behavior of labeling stars like "γ Cas" that the IAU list omits.
    """
    cat = _catdir(cache)
    npz = cat / "bayer-stars.npz"
    if not npz.exists():
        raw_pos = cat / "stars.6.json"
        raw_nam = cat / "starnames.json"
        if not (fetch(URLS["stars.6.json"], raw_pos, "bright star positions")
                and fetch(URLS["starnames.json"], raw_nam, "star designations")):
            return None
        namedb = json.loads(raw_nam.read_text(encoding="utf-8"))
        labels, ras, decs, mags, hips = [], [], [], [], []
        for feat in json.loads(raw_pos.read_text(encoding="utf-8"))["features"]:
            hip = feat.get("id")
            mag = feat.get("properties", {}).get("mag")
            if not isinstance(hip, int) or mag is None:
                continue
            info = namedb.get(str(hip))
            if not info:
                continue
            con = info.get("c", "")
            desig = info.get("bayer") or info.get("flam")
            if not desig or not con:
                continue
            lon, lat = feat["geometry"]["coordinates"]
            labels.append(f"{desig} {con}")
            ras.append(lon % 360.0); decs.append(lat)
            mags.append(float(mag)); hips.append(hip)
        np.savez_compressed(npz, label=np.array(labels, "U12"),
                            ra=np.array(ras, "f8"), dec=np.array(decs, "f8"),
                            mag=np.array(mags, "f4"), hip=np.array(hips, "u4"))
    d = np.load(npz)
    return d["label"], d["ra"], d["dec"], d["mag"], d["hip"]


def _sex2deg(s: str, hours: bool) -> float:
    sign = -1.0 if s.strip().startswith("-") else 1.0
    p = [float(x) for x in s.strip().lstrip("+-").split(":")]
    val = p[0] + p[1] / 60.0 + (p[2] if len(p) > 2 else 0.0) / 3600.0
    return sign * val * (15.0 if hours else 1.0)


def load_ngc(cache: Path):
    """OpenNGC -> (label, common, ra, dec, majax_arcmin, mag)."""
    cat = _catdir(cache)
    npz = cat / "ngc-v2.npz"
    if not npz.exists():
        raw = cat / "openngc.csv"
        if not fetch(URLS["openngc.csv"], raw, "OpenNGC catalog"):
            return None
        labels, commons, ras, decs, majs, mags = [], [], [], [], [], []
        lines = raw.read_text(encoding="utf-8", errors="replace").splitlines()
        cols = {c: i for i, c in enumerate(lines[0].split(";"))}

        def field(row, name):
            i = cols.get(name)
            return row[i].strip() if i is not None and i < len(row) else ""

        def num(row, name, default):
            try:
                return float(field(row, name) or default)
            except ValueError:
                return default

        name_re = re.compile(r"^(NGC|IC)(\d{4})(.*)$")
        for line in lines[1:]:
            row = line.split(";")
            if field(row, "Type") not in NGC_KEEP_TYPES:
                continue
            ra_s, dec_s = field(row, "RA"), field(row, "Dec")
            if not ra_s or not dec_s:
                continue
            m = name_re.match(field(row, "Name"))
            if m:
                label = f"{m.group(1)} {int(m.group(2))}{m.group(3).strip()}"
            else:
                label = field(row, "Name")
            messier = field(row, "M")
            if messier:
                label = f"M {int(messier)}"
            try:
                ra, dec = _sex2deg(ra_s, True), _sex2deg(dec_s, False)
            except ValueError:
                continue
            labels.append(label); commons.append(field(row, "Common names"))
            ras.append(ra); decs.append(dec)
            majs.append(num(row, "MajAx", 0.0))
            mags.append(min(num(row, "V-Mag", 99.0), num(row, "B-Mag", 99.0)))
        np.savez_compressed(npz, label=np.array(labels, "U16"),
                            common=np.array(commons, "U48"),
                            ra=np.array(ras, "f8"), dec=np.array(decs, "f8"),
                            majax=np.array(majs, "f4"),
                            mag=np.array(mags, "f4"))
    d = np.load(npz)
    return d["label"], d["common"], d["ra"], d["dec"], d["majax"], d["mag"]


def load_hd(cache: Path):
    """Henry Draper positions from astrometry.net's hd.fits (libkd kd-tree).

    The tree stores unit-sphere xyz as scaled uint32 with ranges in the
    kdtree_range HDU: x = minval[d] + raw / scale. HD number of tree point i
    is perm[i] + 1. Compiled to hd.npz (ra, dec, hd) on first use.
    """
    cat = _catdir(cache)
    npz = cat / "hd.npz"
    if not npz.exists():
        raw = cat / "hd.fits"
        if not fetch(URLS["hd.fits"], raw, "Henry Draper catalog (hd.fits)"):
            return None
        hdus = {}
        endian = ">"  # libkd writes native order and records it in a card
        with fits.open(raw) as hdul:
            for h in hdul[1:]:
                ttype = h.header.get("TTYPE1", "")
                if ttype == "kdtree_header" and h.header.get("ENDIAN") == "04:03:02:01":
                    endian = "<"
                hdus[ttype] = (np.frombuffer(h.data.tobytes(), dtype=np.uint8),
                               h.header.get("NAXIS1", 0), h.header.get("NAXIS2", 0))
        if "kdtree_data" not in hdus or "kdtree_perm" not in hdus:
            log("  warning: hd.fits has unexpected structure, skipping HD")
            return None
        dbytes, w, n = hdus["kdtree_data"]
        if w == 12:  # 3 x uint32, scaled by kdtree_range
            pts = dbytes.view(endian + "u4").reshape(n, 3).astype("f8")
            rb, _, _ = hdus["kdtree_range"]
            rng = rb.view(endian + "f8")  # minval[3], maxval[3], scale
            pts = rng[:3] + pts / rng[6]
        elif w == 24:
            pts = dbytes.view(endian + "f8").reshape(n, 3)
        else:
            log(f"  warning: hd.fits data width {w} not handled, skipping HD")
            return None
        perm = hdus["kdtree_perm"][0].view(endian + "u4")[:n]
        norm = np.linalg.norm(pts, axis=1)
        if not (0.9 < float(np.median(norm)) < 1.1):
            log("  warning: hd.fits decode failed sanity check, skipping HD")
            return None
        pts = pts / norm[:, None]
        ra = np.degrees(np.arctan2(pts[:, 1], pts[:, 0])) % 360.0
        dec = np.degrees(np.arcsin(np.clip(pts[:, 2], -1, 1)))
        np.savez_compressed(npz, ra=ra, dec=dec,
                            hd=(perm.astype("u4") + 1))
    d = np.load(npz)
    return d["ra"], d["dec"], d["hd"]


# ----------------------------------------------------------------------------
# Image loading + star extraction
# ----------------------------------------------------------------------------

FITS_SUFFIXES = {".fit", ".fits", ".fts", ".fz"}


def load_image(path: Path):
    """Return (gray 2-D float32 for solving, PIL RGB for display, scale hint).

    The scale hint (lo, hi arcsec/pixel) comes from image metadata when
    available -- EXIF 35 mm-equivalent focal length for camera files,
    FOCALLEN + XPIXSZ for FITS -- and greatly narrows the blind search.
    """
    hint = None
    if path.suffix.lower() in FITS_SUFFIXES:
        with fits.open(path) as hdul:
            hdu = next((h for h in hdul if h.data is not None
                        and getattr(h.data, "ndim", 0) >= 2), None)
            if hdu is None:
                sys.exit(f"error: no image data in {path}")
            data = np.asarray(hdu.data, dtype=np.float32)
            focal = hdu.header.get("FOCALLEN") or hdul[0].header.get("FOCALLEN")
            pixsz = hdu.header.get("XPIXSZ") or hdul[0].header.get("XPIXSZ")
        if focal and pixsz:
            app = 206.265 * float(pixsz) / float(focal)
            hint = (app * 0.7, app * 1.4)
        if data.ndim == 3:  # Siril writes (3, H, W); some tools (H, W, 3)
            ax = int(np.argmin(data.shape))
            data = data.mean(axis=ax)
        data = np.nan_to_num(data, nan=float(np.nanmin(data)))
        lo, hi = np.percentile(data, (0.5, 99.8))
        disp = np.clip((data - lo) / max(hi - lo, 1e-9), 0, 1) ** 0.5
        rgb = Image.fromarray((disp * 255).astype(np.uint8), "L").convert("RGB")
        return data, rgb, hint

    img = Image.open(path)  # EXIF orientation deliberately NOT applied:
    img.load()              # pixel coords must match the raw image grid
    try:
        fl35 = img.getexif().get_ifd(0x8769).get(41989)  # FocalLengthIn35mmFilm
        if fl35:
            fov_deg = math.degrees(2 * math.atan(36.0 / 2.0 / float(fl35)))
            app = fov_deg * 3600.0 / img.size[0]
            hint = (app * 0.6, app * 1.6)
    except Exception:
        pass
    if img.mode in ("I;16", "I;16B", "I;16L", "I", "F"):
        data = np.asarray(img, dtype=np.float32)
        rgb = Image.fromarray(
            (np.clip(data / max(data.max(), 1e-9), 0, 1) * 255).astype(np.uint8),
            "L").convert("RGB")
    else:
        rgb = img.convert("RGB")
        data = np.asarray(rgb, dtype=np.float32).mean(axis=2)
    return data, rgb, hint


def _gaussian_smooth(a: np.ndarray, sigma: float = 1.0) -> np.ndarray:
    """Separable Gaussian blur (matched filter for star detection)."""
    r = max(2, int(3.0 * sigma))
    x = np.arange(-r, r + 1, dtype=np.float32)
    k = np.exp(-0.5 * (x / sigma) ** 2)
    k /= k.sum()
    out = a
    for axis in (0, 1):
        pad = [(0, 0), (0, 0)]
        pad[axis] = (r, r)
        ap = np.pad(out, pad, mode="reflect")
        acc = np.zeros_like(a)
        sl = [slice(None), slice(None)]
        for j, wgt in enumerate(k):
            sl[axis] = slice(j, j + a.shape[axis])
            acc += wgt * ap[tuple(sl)]
        out = acc
    return out


def extract_stars(gray: np.ndarray, downsample: int, max_objs: int,
                  threshold: float):
    """Detect stars -> (x, y) in full-resolution pixel coords.

    Mirrors solve-field's simplexy strategy: detect on a PSF-smoothed image
    (matched filter -- suppresses both single-pixel noise spikes and extended
    foreground blobs), rank by smoothed peak, and spatially uniformize the
    selection so every image region contributes candidates.
    """
    work = gray
    if downsample > 1:
        h, w = gray.shape
        h2, w2 = h - h % downsample, w - w % downsample
        work = gray[:h2, :w2].reshape(
            h2 // downsample, downsample, w2 // downsample, downsample
        ).mean(axis=(1, 3))
    work = np.ascontiguousarray(work, dtype=np.float32)

    sep.set_extract_pixstack(3_000_000)
    sep.set_sub_object_limit(8192)
    bkg = sep.Background(work)
    smooth = _gaussian_smooth(work - bkg.back(), sigma=1.0)
    err = float(sep.Background(smooth).globalrms) or 1.0
    try:
        objs = sep.extract(smooth, thresh=threshold, err=err, minarea=5,
                           filter_kernel=None)
    except Exception:
        # dense fields can overflow the deblender; merged neighbors are
        # perfectly fine centroids for plate solving
        objs = sep.extract(smooth, thresh=threshold, err=err, minarea=5,
                           filter_kernel=None, deblend_cont=1.0)
    if len(objs) == 0:
        return np.empty((0, 2))
    round_enough = objs["a"] / np.maximum(objs["b"], 1e-6) < 3.0
    objs = objs[round_enough] if round_enough.any() else objs

    # Uniformize: round-robin the brightest star of every grid cell, then the
    # second-brightest, ... so foreground or edge artifacts in one region
    # cannot crowd real stars out of the candidate list.
    G = 10
    gx = np.clip((objs["x"] / work.shape[1] * G).astype(int), 0, G - 1)
    gy = np.clip((objs["y"] / work.shape[0] * G).astype(int), 0, G - 1)
    cell = gy * G + gx
    rank_in_cell = np.empty(len(objs), dtype=int)
    seen: dict[int, int] = {}
    for pos in np.argsort(objs["peak"])[::-1]:
        c = int(cell[pos])
        rank_in_cell[pos] = seen.get(c, 0)
        seen[c] = rank_in_cell[pos] + 1
    order = np.lexsort((-objs["peak"], rank_in_cell))[:max_objs]
    objs = objs[order]
    # map downsampled pixel centers back to full-res coordinates
    xy = np.column_stack([
        objs["x"] * downsample + (downsample - 1) / 2.0,
        objs["y"] * downsample + (downsample - 1) / 2.0,
    ])
    return xy


# ----------------------------------------------------------------------------
# Solving
# ----------------------------------------------------------------------------

def solve(xy: np.ndarray, index_files: list[Path], args, auto_hint=None):
    size_hint = position_hint = None
    if args.scale_low or args.scale_high:
        size_hint = astrometry.SizeHint(
            lower_arcsec_per_pixel=args.scale_low or 0.05,
            upper_arcsec_per_pixel=args.scale_high or 1000.0)
    elif auto_hint and not args.no_auto_hint:
        size_hint = astrometry.SizeHint(
            lower_arcsec_per_pixel=auto_hint[0],
            upper_arcsec_per_pixel=auto_hint[1])
        log(f"  scale hint from image metadata: "
            f"{auto_hint[0]:.1f}-{auto_hint[1]:.1f} arcsec/px")
    if args.ra is not None and args.dec is not None:
        position_hint = astrometry.PositionHint(
            ra_deg=args.ra, dec_deg=args.dec, radius_deg=args.radius)

    base_params = dict(
        sip_order=args.sip_order,
        output_logodds_threshold=math.log(1e9),
    )
    if args.sip_order == 0:
        base_params["tune_up_logodds_threshold"] = None

    if args.timeout > 0:
        # The engine call is not interruptible in-process, so a hard
        # deadline requires running it in a worker we can kill.
        recv, send = multiprocessing.Pipe(duplex=False)
        proc = multiprocessing.Process(
            target=_solve_worker,
            args=(send, index_files, xy.tolist(), size_hint, position_hint,
                  base_params))
        proc.start()
        send.close()
        if recv.poll(args.timeout):
            payload = recv.recv()
            proc.join()
        else:
            proc.terminate()
            proc.join(2)
            if proc.is_alive():
                proc.kill()
                proc.join()
            sys.exit(f"error: no solution within --timeout {args.timeout:g}s "
                     "(hints help: --scale-low/--scale-high, --ra/--dec)")
        if isinstance(payload, dict) and "error" in payload:
            sys.exit(f"error: solver worker failed: {payload['error']}")
        return astrometry.Solution.from_json(payload)

    params = dict(base_params,
                  logodds_callback=lambda logodds: astrometry.Action.STOP)
    solver = astrometry.Solver(index_files)
    return solver.solve(
        stars=xy.tolist(),
        size_hint=size_hint,
        position_hint=position_hint,
        solution_parameters=astrometry.SolutionParameters(**params),
    )


def _solve_worker(conn, index_files, stars, size_hint, position_hint,
                  base_params):
    """Engine run in a child process so --timeout can hard-kill it."""
    try:
        params = dict(base_params,
                      logodds_callback=lambda logodds: astrometry.Action.STOP)
        solver = astrometry.Solver(index_files)
        solution = solver.solve(
            stars=stars,
            size_hint=size_hint,
            position_hint=position_hint,
            solution_parameters=astrometry.SolutionParameters(**params),
        )
        conn.send(solution.to_json())
    except Exception as e:
        conn.send({"error": f"{type(e).__name__}: {e}"})


def wcs_orientation(w: WCS, cx: float, cy: float):
    """(pixscale arcsec/px, east-of-north rotation deg, parity) at image center."""
    sky = w.all_pix2world([[cx, cy]], 0)[0]
    ra, dec = float(sky[0]), float(sky[1])
    step = 0.02
    pn = w.all_world2pix([[ra, min(dec + step, 89.9)]], 0, quiet=True)[0]
    dxn, dyn = pn[0] - cx, pn[1] - cy
    orient = math.degrees(math.atan2(-dxn, dyn))  # wcsinfo convention, E of N
    cd = w.pixel_scale_matrix
    pixscale = 3600.0 * math.sqrt(abs(np.linalg.det(cd)))
    parity = "neg" if np.linalg.det(cd) > 0 else "pos"  # astronomical convention
    return ra, dec, pixscale, orient, parity


def write_new_fits(input_path: Path, gray: np.ndarray, wcs_header,
                   dest: Path) -> None:
    """solve-field's .new equivalent: the input image with the WCS embedded.

    FITS input keeps its original data and header (WCS cards merged in);
    raster input is written as a float32 grayscale FITS.
    """
    if input_path.suffix.lower() in FITS_SUFFIXES:
        with fits.open(input_path) as hdul:
            hdu = next(h for h in hdul if h.data is not None
                       and getattr(h.data, "ndim", 0) >= 2)
            data, header = hdu.data.copy(), hdu.header.copy()
    else:
        data, header = gray.astype(np.float32), fits.Header()
    header.update(wcs_header)
    fits.PrimaryHDU(data=data, header=header).writeto(dest, overwrite=True)


def hms(ra: float) -> str:
    h = ra / 15.0
    m = (h % 1) * 60
    return f"{int(h):02d}:{int(m):02d}:{(m % 1) * 60:06.3f}"


def dms(dec: float) -> str:
    s = "+" if dec >= 0 else "-"
    d = abs(dec)
    m = (d % 1) * 60
    return f"{s}{int(d):02d}:{int(m):02d}:{(m % 1) * 60:06.3f}"


# ----------------------------------------------------------------------------
# Annotation
# ----------------------------------------------------------------------------

def radec_to_xyz(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec)
    return np.stack([np.cos(dec) * np.cos(ra),
                     np.cos(dec) * np.sin(ra),
                     np.sin(dec)], axis=-1)


def in_field_mask(ra, dec, ctr_ra, ctr_dec, radius_deg):
    v = radec_to_xyz(np.asarray(ra), np.asarray(dec))
    c = radec_to_xyz(ctr_ra, ctr_dec)
    return v @ c >= math.cos(math.radians(radius_deg))


def world2pix(w: WCS, ra, dec):
    pts = np.column_stack([np.atleast_1d(ra), np.atleast_1d(dec)])
    try:
        out = w.all_world2pix(pts, 0, quiet=True, maxiter=30)
    except Exception:
        out = w.wcs_world2pix(pts, 0)
    return np.asarray(out)


class Annotator:
    def __init__(self, base: Image.Image, w: WCS, transparent: bool,
                 font_size: int, line_width: int):
        self.W, self.H = base.size
        self.wcs = w
        if transparent:
            self.img = Image.new("RGBA", (self.W, self.H), (0, 0, 0, 0))
        else:
            self.img = base.convert("RGB")
        self.draw = ImageDraw.Draw(self.img)
        diag = math.hypot(self.W, self.H)
        self.lw = line_width or max(2, round(diag / 1800))
        fs = font_size or max(14, round(diag / 110))
        self.font = self._font(fs)
        self.font_small = self._font(max(11, int(fs * 0.72)))
        self.label_boxes: list[tuple[float, float, float, float]] = []

    @staticmethod
    def _font(size: int):
        for cand in ("DejaVuSans.ttf", "Arial.ttf", "FreeSans.ttf",
                     "/usr/share/fonts/TTF/DejaVuSans.ttf",
                     "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"):
            try:
                return ImageFont.truetype(cand, size)
            except OSError:
                continue
        try:
            return ImageFont.load_default(size=size)
        except TypeError:
            return ImageFont.load_default()

    def on_canvas(self, x, y, margin=0.0):
        return -margin <= x < self.W + margin and -margin <= y < self.H + margin

    def line(self, p1, p2, color):
        self.draw.line([p1, p2], fill=(0, 0, 0), width=self.lw + 2)
        self.draw.line([p1, p2], fill=color, width=self.lw)

    def circle(self, xy, r, color):
        bb = [xy[0] - r, xy[1] - r, xy[0] + r, xy[1] + r]
        self.draw.ellipse(bb, outline=(0, 0, 0), width=self.lw + 2)
        self.draw.ellipse(bb, outline=color, width=self.lw)

    def text(self, xy, txt, color, small=False, force=False):
        font = self.font_small if small else self.font
        bbox = self.draw.textbbox(xy, txt, font=font)
        if not force:
            for b in self.label_boxes:
                if not (bbox[2] < b[0] or bbox[0] > b[2]
                        or bbox[3] < b[1] or bbox[1] > b[3]):
                    return False
        self.label_boxes.append(bbox)
        o = max(1, self.lw // 2)
        for dx, dy in ((-o, 0), (o, 0), (0, -o), (0, o),
                       (-o, -o), (o, o), (-o, o), (o, -o)):
            self.draw.text((xy[0] + dx, xy[1] + dy), txt, font=font, fill=(0, 0, 0))
        self.draw.text(xy, txt, font=font, fill=color)
        return True


def subdivide_segment(ra1, dec1, ra2, dec2, step_deg=2.0):
    """Great-circle interpolation so long segments follow projection curvature."""
    v1, v2 = radec_to_xyz(ra1, dec1), radec_to_xyz(ra2, dec2)
    ang = math.acos(float(np.clip(np.dot(v1, v2), -1, 1)))
    n = max(2, int(math.degrees(ang) / step_deg) + 1)
    ts = np.linspace(0, 1, n)
    if ang < 1e-9:
        pts = np.outer(1 - ts, v1) + np.outer(ts, v2)
    else:
        pts = (np.outer(np.sin((1 - ts) * ang), v1)
               + np.outer(np.sin(ts * ang), v2)) / math.sin(ang)
    pts /= np.linalg.norm(pts, axis=1)[:, None]
    ra = np.degrees(np.arctan2(pts[:, 1], pts[:, 0])) % 360
    dec = np.degrees(np.arcsin(np.clip(pts[:, 2], -1, 1)))
    return ra, dec


def annotate(ann: Annotator, cache: Path, ctr_ra: float, ctr_dec: float,
             pixscale: float, args, detected_xy: np.ndarray, counts: dict):
    radius = math.hypot(ann.W, ann.H) * pixscale / 3600.0 / 2.0 * 1.15

    # --- constellation lines + names ---------------------------------------
    if not args.no_constellations:
        res = load_constellation_lines(cache)
        if res:
            segs, seg_con, ids = res
            m = (in_field_mask(segs[:, 0], segs[:, 1], ctr_ra, ctr_dec, radius)
                 | in_field_mask(segs[:, 2], segs[:, 3], ctr_ra, ctr_dec, radius))
            verts: dict[int, list] = {}
            for (ra1, dec1, ra2, dec2), ci in zip(segs[m], seg_con[m]):
                rr, dd = subdivide_segment(ra1, dec1, ra2, dec2)
                px = world2pix(ann.wcs, rr, dd)
                drew = False
                for p1, p2 in zip(px, px[1:]):
                    if (ann.on_canvas(*p1, margin=ann.W) and
                            ann.on_canvas(*p2, margin=ann.W)):
                        ann.line(tuple(p1), tuple(p2), COLORS["lines"])
                        drew = True
                if drew:
                    counts["constellation_segments"] += 1
                    for p in (px[0], px[-1]):
                        if ann.on_canvas(*p):
                            verts.setdefault(int(ci), []).append(p)
            for ci, pts in verts.items():
                if len(pts) < 2:
                    continue
                c = np.mean(pts, axis=0)
                name = CONSTELLATION_NAMES.get(ids[ci], ids[ci])
                tw = ann.draw.textlength(name, font=ann.font)
                if ann.text((c[0] - tw / 2, c[1]), name, COLORS["conname"]):
                    counts["constellations"] += 1

    # --- named bright stars --------------------------------------------------
    iau_hips: set[int] = set()
    if not args.no_bright:
        res = load_bright_stars(cache)
        if res:
            names, desigs, ras, decs, hips = res
            iau_hips = {int(h) for h in hips if h}
            m = in_field_mask(ras, decs, ctr_ra, ctr_dec, radius)
            px = world2pix(ann.wcs, ras[m], decs[m])
            r = max(6, ann.lw * 4)
            for (x, y), nm, dg in zip(px, names[m], desigs[m]):
                if not ann.on_canvas(x, y):
                    continue
                ann.circle((x, y), r, COLORS["star"])
                label = f"{nm} / {dg}" if dg else str(nm)
                if ann.text((x + r * 1.3, y - r), label, COLORS["star"]):
                    counts["bright_stars"] += 1

    # --- Bayer/Flamsteed-only bright stars (no IAU proper name) --------------
    if not args.no_bright and args.bright_mag > 0:
        res = load_bayer_stars(cache)
        if res:
            labels, ras, decs, mags, hips = res
            keep = (mags <= args.bright_mag) & ~np.isin(hips, list(iau_hips))
            m = in_field_mask(ras[keep], decs[keep], ctr_ra, ctr_dec, radius)
            px = world2pix(ann.wcs, ras[keep][m], decs[keep][m])
            r = max(5, ann.lw * 3)
            for (x, y), lb in zip(px, labels[keep][m]):
                if not ann.on_canvas(x, y):
                    continue
                ann.circle((x, y), r, COLORS["star"])
                if ann.text((x + r * 1.3, y - r), str(lb), COLORS["star"],
                            small=True):
                    counts["bayer_stars"] += 1

    # --- NGC / IC / Messier ---------------------------------------------------
    if not args.no_ngc:
        res = load_ngc(cache)
        if res:
            labels, commons, ras, decs, majs, mags = res
            # keep objects that are plausibly visible/interesting: bright
            # enough, named, Messier, or with substantial angular extent
            keep = ((mags <= args.ngc_mag) | (commons != "") | (majs >= 5.0)
                    | np.char.startswith(labels, "M "))
            labels, commons, ras, decs, majs = (
                labels[keep], commons[keep], ras[keep], decs[keep], majs[keep])
            m = in_field_mask(ras, decs, ctr_ra, ctr_dec, radius)
            px = world2pix(ann.wcs, ras[m], decs[m])
            for (x, y), lb, cm, mj in zip(px, labels[m], commons[m], majs[m]):
                if not ann.on_canvas(x, y):
                    continue
                r = max(8.0, float(mj) * 60.0 / 2.0 / pixscale)
                if r > 0.9 * max(ann.W, ann.H):  # absurdly large overlay
                    continue
                ann.circle((x, y), r, COLORS["ngc"])
                label = f"{lb} ({cm.split(',')[0]})" if cm else str(lb)
                if ann.text((x + r * 0.74 + 4, y + r * 0.74 + 4), label,
                            COLORS["ngc"], small=True):
                    counts["ngc"] += 1

    # --- Henry Draper labels (--hd) -------------------------------------------
    if args.hd:
        res = load_hd(cache)
        if res:
            ras, decs, hd_num = res
            m = in_field_mask(ras, decs, ctr_ra, ctr_dec, radius)
            ras, decs, hd_num = ras[m], decs[m], hd_num[m]
            px = world2pix(ann.wcs, ras, decs)
            on = np.array([ann.on_canvas(x, y) for x, y in px])
            px, hd_num = px[on], hd_num[on]
            if len(px) > args.hd_max:
                # prefer HD stars nearest to actual detections (likely visible)
                if len(detected_xy):
                    d2 = ((px[:, None, :] - detected_xy[None, :500, :]) ** 2
                          ).sum(axis=2).min(axis=1)
                    keep = np.argsort(d2)[:args.hd_max]
                else:
                    keep = np.arange(args.hd_max)
                px, hd_num = px[keep], hd_num[keep]
                log(f"  note: {int(on.sum())} HD stars in field, "
                    f"labeling {args.hd_max} (raise with --hd-max)")
            r = max(3, ann.lw * 2)
            for (x, y), n in zip(px, hd_num):
                ann.circle((x, y), r, COLORS["hd"])
                if ann.text((x + r + 3, y + r + 3), f"HD {int(n)}",
                            COLORS["hd"], small=True):
                    counts["hd"] += 1


# ----------------------------------------------------------------------------
# Diagnostics (--check)
# ----------------------------------------------------------------------------

def check_catalogs(cache: Path, with_hd: bool) -> None:
    ok = True
    res = load_constellation_lines(cache)
    if res:
        segs, seg_con, ids = res
        lyr = "Lyr" in ids
        log(f"  constellation lines: {len(segs)} segments, "
            f"{len(ids)} constellations, Lyra present: {lyr}")
        ok &= lyr
    else:
        ok = False
    res = load_bright_stars(cache)
    if res:
        names, desigs, ras, decs, hips = res
        i = list(names).index("Vega") if "Vega" in names else -1
        log(f"  bright stars: {len(names)} named")
        if i >= 0:
            d = math.hypot((ras[i] - 279.2347) * math.cos(math.radians(38.78)),
                           decs[i] - 38.7837) * 3600
            log(f"    Vega: {desigs[i]} at RA {ras[i]:.4f} Dec {decs[i]:+.4f} "
                f"(offset {d:.1f} arcsec, expect < 5)")
            ok &= d < 5
        else:
            ok = False
    else:
        ok = False
    res = load_bayer_stars(cache)
    if res:
        labels, ras, decs, mags, hips = res
        i = np.where(hips == 4427)[0]  # gamma Cas: bright, no IAU name
        log(f"  bayer/flamsteed stars: {len(labels)}")
        if len(i):
            i = int(i[0])
            log(f"    HIP 4427 -> '{labels[i]}' mag {mags[i]:.2f} "
                f"(expect 'γ Cas' ~2.2)")
            ok &= str(labels[i]) == "γ Cas" and mags[i] < 3.0
        else:
            ok = False
    else:
        ok = False
    res = load_ngc(cache)
    if res:
        labels, commons, ras, decs, majs, mags = res
        idx = [i for i, l in enumerate(labels) if l == "M 57"]
        log(f"  NGC/IC objects: {len(labels)}")
        if idx:
            i = idx[0]
            d = math.hypot((ras[i] - 283.396) * math.cos(math.radians(33.03)),
                           decs[i] - 33.029) * 3600
            log(f"    M 57 at RA {ras[i]:.4f} Dec {decs[i]:+.4f} "
                f"(offset {d:.1f} arcsec, expect < 60)")
            ok &= d < 60
        else:
            ok = False
    else:
        ok = False
    if with_hd:
        res = load_hd(cache)
        if res:
            ras, decs, hd_num = res
            log(f"  HD stars: {len(ras)}")
            for name, hra, hdec, hd_expect in (
                    ("Vega", 279.2347, 38.7837, 172167),
                    ("Sirius", 101.2872, -16.7161, 48915)):
                v = radec_to_xyz(ras, decs) @ radec_to_xyz(hra, hdec)
                i = int(np.argmax(v))
                sep_as = math.degrees(math.acos(float(np.clip(v[i], -1, 1)))) * 3600
                log(f"    {name}: nearest HD {int(hd_num[i])} at {sep_as:.1f} arcsec "
                    f"(expect HD {hd_expect})")
                ok &= int(hd_num[i]) == hd_expect and sep_as < 30
        else:
            ok = False
    log("catalog check: " + ("OK" if ok else "FAILED"))
    sys.exit(0 if ok else 1)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="psa.py",
        description="Portable plate-solve + annotate (astrometry.net engine, "
                    "self-bootstrapping index files and catalogs).",
        epilog="examples:\n"
               "  psa.py result.fit                 solve + annotate\n"
               "  psa.py vega.jpg --hd              with Henry Draper labels\n"
               "  psa.py --prefetch --hd            warm cache for offline use\n"
               "  psa.py img.jpg --ra 279 --dec 38 --radius 10   hinted (faster)\n",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("image", nargs="?", help="FITS/JPEG/PNG/TIFF image to solve")
    p.add_argument("--hd", action="store_true",
                   help="annotate Henry Draper catalog stars")
    p.add_argument("--out", metavar="DIR",
                   help="output directory (default: '<name> Solved' in CWD)")
    p.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE,
                   help=f"index/catalog cache (default {DEFAULT_CACHE})")
    p.add_argument("--downsample", type=int, default=2, metavar="N",
                   help="downsample factor before extraction (default 2)")
    p.add_argument("--objs", type=int, default=1000, metavar="N",
                   help="max detected stars passed to solver (default 1000)")
    p.add_argument("--threshold", type=float, default=5.0, metavar="SIGMA",
                   help="extraction threshold in background sigma (default 5)")
    p.add_argument("--series", default="4100,4200",
                   help="index series (default 4100,4200)")
    p.add_argument("--scales", default="11-19",
                   help="index scales, e.g. 11-19 or 8,9,10 (default 11-19)")
    p.add_argument("--sip-order", type=int, default=3,
                   help="SIP distortion order, 0 disables (default 3)")
    p.add_argument("--timeout", type=float, default=0, metavar="SEC",
                   help="hard solve deadline; the engine runs in a worker "
                        "process killed at the deadline (default 0 = none)")
    p.add_argument("--write-new", action="store_true",
                   help="also write <name>.new.fits -- the input image with "
                        "the WCS embedded (like solve-field's .new)")
    p.add_argument("--ra", type=float, help="position hint RA (deg)")
    p.add_argument("--dec", type=float, help="position hint Dec (deg)")
    p.add_argument("--radius", type=float, default=15.0,
                   help="position hint radius (deg, default 15)")
    p.add_argument("--scale-low", type=float, metavar="APP",
                   help="scale hint lower bound (arcsec/pixel)")
    p.add_argument("--scale-high", type=float, metavar="APP",
                   help="scale hint upper bound (arcsec/pixel)")
    p.add_argument("--no-auto-hint", action="store_true",
                   help="ignore scale hints derived from EXIF/FITS metadata")
    p.add_argument("--transparent", action="store_true",
                   help="annotations on transparent background "
                        "(plot-constellations style) instead of over the image")
    p.add_argument("--font-size", type=int, default=0,
                   help="label font px (default: scaled to image)")
    p.add_argument("--line-width", type=int, default=0,
                   help="line width px (default: scaled to image)")
    p.add_argument("--hd-max", type=int, default=1000,
                   help="max HD labels (default 1000)")
    p.add_argument("--ngc-mag", type=float, default=12.0,
                   help="NGC/IC magnitude cutoff; named, Messier, and large "
                        "objects are always kept (default 12)")
    p.add_argument("--bright-mag", type=float, default=4.0,
                   help="label Bayer/Flamsteed stars (no IAU name) down to "
                        "this magnitude; 0 disables (default 4.0)")
    p.add_argument("--no-constellations", action="store_true")
    p.add_argument("--no-bright", action="store_true")
    p.add_argument("--no-ngc", action="store_true")
    p.add_argument("--no-annotate", action="store_true",
                   help="plate solve only, skip annotation image")
    p.add_argument("--prefetch", action="store_true",
                   help="download index files + catalogs, then exit")
    p.add_argument("--check", action="store_true",
                   help="verify cached catalogs against known objects, then exit")
    p.add_argument("--version", action="version", version=f"psa.py {VERSION}")
    return p


def main() -> None:
    args = build_parser().parse_args()
    cache: Path = args.cache_dir
    cache.mkdir(parents=True, exist_ok=True)

    if args.check:
        check_catalogs(cache, with_hd=True)

    if args.prefetch:
        log("Prefetching index files ...")
        files = get_index_files(cache, args.series, parse_scales(args.scales))
        log(f"  {len(files)} index files ready")
        log("Prefetching annotation catalogs ...")
        load_constellation_lines(cache)
        load_bright_stars(cache)
        load_bayer_stars(cache)
        load_ngc(cache)
        if args.hd:
            load_hd(cache)
        log(f"Cache ready at {cache} -- psa.py now runs offline.")
        return

    if not args.image:
        build_parser().print_help()
        sys.exit(1)

    input_path = Path(args.image)
    if not input_path.exists():
        sys.exit(f"error: {input_path} not found")
    out_dir = Path(args.out) if args.out else Path(f"{input_path.stem} Solved")
    out_dir.mkdir(parents=True, exist_ok=True)
    wcs_file = out_dir / f"{input_path.stem}.wcs"
    ann_file = out_dir / "annotations.png"

    log(f"Running plate solving on: {input_path}")
    log(f"Output directory: {out_dir}")

    log("Step 1: solving ...")
    t0 = time.time()
    gray, rgb, auto_hint = load_image(input_path)
    h, w = gray.shape
    xy = extract_stars(gray, max(1, args.downsample), args.objs, args.threshold)
    log(f"  extracted {len(xy)} sources "
        f"({w}x{h}, downsample {args.downsample})")
    if len(xy) < 8:
        sys.exit("error: too few stars detected -- adjust --threshold/--downsample?")

    index_files = get_index_files(cache, args.series, parse_scales(args.scales))
    solution = solve(xy, index_files, args, auto_hint)
    if not solution.has_match():
        sys.exit("error: no solution found (try different --scales/--series, "
                 "more --objs, or position/scale hints)")
    best = solution.best_match()
    sol_wcs = best.astropy_wcs()
    ra_c, dec_c, pixscale, orient, parity = wcs_orientation(
        sol_wcs, (w - 1) / 2.0, (h - 1) / 2.0)
    t_solve = time.time() - t0

    log(f"  solved in {t_solve:.1f}s using {Path(best.index_path).name}")
    log(f"  center: RA {hms(ra_c)} ({ra_c:.6f}), Dec {dms(dec_c)} ({dec_c:+.6f})")
    log(f"  pixel scale: {pixscale:.4f} arcsec/px, "
        f"field {w * pixscale / 3600:.3g} x {h * pixscale / 3600:.3g} deg")
    log(f"  rotation: up is {orient:.2f} deg E of N, parity {parity}")

    hdr = sol_wcs.to_header(relax=True)
    hdr["IMAGEW"] = (w, "image width in pixels")
    hdr["IMAGEH"] = (h, "image height in pixels")
    hdr["COMMENT"] = f"solved by psa.py {VERSION} (astrometry.net engine)"
    fits.PrimaryHDU(header=hdr).writeto(wcs_file, overwrite=True)
    log(f"  wrote {wcs_file}")
    if args.write_new:
        new_file = out_dir / f"{input_path.stem}.new.fits"
        write_new_fits(input_path, gray, hdr, new_file)
        log(f"  wrote {new_file}")

    counts = {"constellations": 0, "constellation_segments": 0,
              "bright_stars": 0, "bayer_stars": 0, "ngc": 0, "hd": 0}
    if not args.no_annotate:
        log("Step 2: annotating ...")
        ann = Annotator(rgb, sol_wcs, args.transparent,
                        args.font_size, args.line_width)
        annotate(ann, cache, ra_c, dec_c, pixscale, args, xy, counts)
        ann.img.save(ann_file)
        log(f"  drew {counts['constellations']} constellations, "
            f"{counts['bright_stars']} named + {counts['bayer_stars']} bayer "
            f"stars, {counts['ngc']} NGC/IC"
            + (f", {counts['hd']} HD" if args.hd else ""))
        log(f"  wrote {ann_file}")

    summary = {
        "input": str(input_path), "image_size": [w, h],
        "stars_extracted": int(len(xy)), "solve_seconds": round(t_solve, 2),
        "ra_center_deg": ra_c, "dec_center_deg": dec_c,
        "ra_center_hms": hms(ra_c), "dec_center_dms": dms(dec_c),
        "pixscale_arcsec": pixscale,
        "field_deg": [w * pixscale / 3600.0, h * pixscale / 3600.0],
        "rotation_deg_E_of_N": orient, "parity": parity,
        "index": Path(best.index_path).name,
        "logodds": float(getattr(best, "logodds", float("nan"))),
        "annotations": counts,
    }
    (out_dir / "solution.json").write_text(json.dumps(summary, indent=2) + "\n")

    log("Success! Files created:")
    log(f"  - WCS solution:  {wcs_file}")
    if not args.no_annotate:
        log(f"  - Annotations:   {ann_file}")
    log(f"  - Summary:       {out_dir / 'solution.json'}")


if __name__ == "__main__":
    main()
