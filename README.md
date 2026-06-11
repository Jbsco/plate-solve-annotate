# psa-port

Self-contained, portable Python port of `plate-solve-annotate.sh` (the `psa`
alias): blind plate solving + sky annotation in a single file, no system
astrometry.net install required.

```
./psa.py image.fit            # solve + annotate
./psa.py image.jpg --hd       # also label Henry Draper stars
./psa.py --prefetch --hd      # warm the cache for offline use
./psa.py --check              # verify cached catalogs
```

The shebang runs the script through [uv](https://docs.astral.sh/uv/), which
resolves the inline dependencies (PEP 723) into a cached environment
automatically — nothing to install beyond uv itself. `uv run psa.py ...`
works too, as does a plain venv with
`pip install 'astrometry>=4.1.2,<5' astropy numpy sep pillow`.

## How it replaces the system packages

| original (`psa` alias) | this script |
|---|---|
| `solve-field` + `astrometry-engine` (C binaries, AUR package) | [`astrometry`](https://github.com/neuromorphicsystems/astrometry) wheel — same astrometry.net engine compiled as a Python extension |
| `image2pnm` + netpbm conversion (Python 3.13-locked, **currently broken on this system under Python 3.14**) | `astropy.io.fits` / Pillow |
| simplexy source extraction | `sep` (Source-Extractor), matched-filter detection + spatial uniformization |
| index files hand-managed in `/usr/share/astrometry/data` (`wget.sh`) | downloaded programmatically into the cache on first run |
| `plot-constellations` (built-in catalogs + `hd.fits`) | Pillow drawing through the solved WCS; catalogs auto-downloaded: d3-celestial constellation lines, IAU-CSN star names, OpenNGC, astrometry.net `hd.fits` |

First run downloads ~45 MB into `~/.cache/psa/` (override with
`--cache-dir` or pre-warm with `--prefetch`); afterwards everything runs
fully offline (verified with all HTTP traffic blocked).

Default index scales mirror the wide-field rig this replaces:
series 4100 + 4200, scales 11–19 (quad sizes ≈ 1.4°–33°, suiting fields
from a few degrees up to all-sky). For narrower fields download deeper
scales, e.g. `--series 4200 --scales 7-12`, or Gaia-based
`--series 5200 --scales 0-6` for telescope-scale FOVs.

## Outputs

Written to `"<name> Solved/"` (same convention as the shell script):

- `<name>.wcs` — FITS WCS header; **readable by astrometry.net tools**
  (`wcsinfo`, `plot-constellations -w`, etc.)
- `annotations.png` — annotated image. By default annotations are composited
  over the (stretched) source image; `--transparent` reproduces the original
  overlay-only behavior.
- `solution.json` — machine-readable summary (center, scale, rotation,
  field size, matched index, annotation counts).

## Annotation layers

- constellation stick figures + names (`--no-constellations` to disable)
- IAU named bright stars, labeled `Name / α Con` (`--no-bright`)
- NGC / IC / Messier objects with extent circles (`--no-ngc`); faint
  anonymous entries are suppressed — kept if mag ≤ `--ngc-mag` (default 12),
  named, Messier, or ≥ 5′ across
- `--hd`: Henry Draper numbers (capped by `--hd-max`, nearest-to-detected
  stars preferred)

Label/line sizes scale with image resolution; override with `--font-size` /
`--line-width`.

## Performance hints

Scale hints make wide-field solves dramatically faster. They are applied
automatically when the file carries metadata:

- JPEG/TIFF: EXIF 35 mm-equivalent focal length
- FITS: `FOCALLEN` + `XPIXSZ` headers (Siril/Ekos write these)

Without metadata (e.g. WhatsApp re-encodes strip EXIF), a 60°+ field can
take many minutes blind — same as the original solve-field. Supply hints:

```
./psa.py img.jpg --scale-low 40 --scale-high 90     # arcsec/pixel bounds
./psa.py img.jpg --ra 332 --dec 47 --radius 25      # position hint
```

`--no-auto-hint` disables metadata hints (true blind solve).

## Validation (2026-06-11)

Tested against reference solve-field/plot-constellations outputs in
`/srv/samba/public/astrophotography`:

| image | field | solve time | vs reference |
|---|---|---|---|
| cassiopeia.jpg (12455×8250) | 10.4°×6.9° | 22 s blind | center Δ≈5″, scale Δ0.08%, corners ≤34″ |
| m57_8_30_stacked.jpg | 10.5°×7.0° | 13 s blind | center Δ≈5″ |
| PXL night shot (EXIF hint) | 54°×73° | **1.8 s** | center Δ≈5 px, scale Δ0.6% |
| Siril stack .fits (FITS hint) | 15.7°×10.5° | **1.6 s** | center Δ≈5″ vs same-night solve |
| WhatsApp shot, 40% trees, no EXIF | 73°×55° | 8.5 min with `--scale-low/high` | solved; original solve-field **failed** after 29 min CPU |

Catalog spot checks (`--check`): Vega → HD 172167 (3″), Sirius → HD 48915,
M 57 position 1.5″, 451 IAU names, 88 constellations.

## Platform notes

- Linux/macOS wheels for the `astrometry` engine; no native Windows (use
  WSL). On Raspberry Pi (aarch64) it builds from source — needs a C
  compiler.
- The original image orientation is used (EXIF rotation NOT applied), so
  pixel coordinates in the `.wcs` match the raw file, like solve-field.
- Rotation is reported in `wcsinfo`'s convention (east of north).

## Suggested alias

```bash
alias psa2='~/cs/psa-port/psa.py'
```
