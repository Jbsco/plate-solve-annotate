# plate-solve-annotate

Self-contained, portable plate solving + sky annotation for astrophotography
in a single Python file. Blind-solves an image against the astrometry.net
index files, writes a standard WCS solution, and renders constellation
figures, star names, NGC/IC/Messier objects, and optional Henry Draper
labels over the image — with no system astrometry.net installation.

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

## Background: the pipeline this replaces

Local plate solving is conventionally assembled from a distribution package
of [astrometry.net](http://astrometry.net) plus a small shell wrapper, along
the lines of:

```bash
#!/bin/bash
# typical local plate-solve-and-annotate wrapper
out="${1%.*} Solved"
solve-field "$1" --downsample 2 --objs 1000 --tag-all -D "$out"
plot-constellations -w "$out/$(basename "${1%.*}").wcs" -N -C -B \
    -D -d /usr/share/astrometry/data/hd.fits -o "$out/annotations.png"
```

That setup works well — until the system around it moves. Failure modes
observed in practice with exactly such a wrapper:

- **Hidden interpreter coupling.** `solve-field` shells out to `image2pnm`,
  a Python helper installed into a *versioned* site-packages directory. A
  routine system Python upgrade (3.13 → 3.14) orphaned that module, and
  every solve — including FITS input, which passes through the same
  file-type sniffing — began failing with `ModuleNotFoundError`. The
  compiled binaries were untouched and the package manager reported nothing
  wrong, but the pipeline was dead until rebuilt against the new Python.
- **Hand-maintained survey data.** Index files and `hd.fits` accumulate in
  a system directory, fetched by ad-hoc `wget` scripts, with mixed
  root/user file ownership and no record of which index scales were chosen
  or why. Reproducing the setup on a second machine is archaeology, and
  nothing updates or validates the data afterwards.
- **A wide native dependency surface.** cfitsio, wcslib, GSL, cairo,
  netpbm, libcurl — a chain of C libraries that must stay in step with the
  OS, and a rebuild burden whenever one of them bumps.

This script keeps the proven solver core and replaces everything around it:

| system-install pipeline | this script |
|---|---|
| `solve-field` + `astrometry-engine` (C binaries from a distro package) | [`astrometry`](https://github.com/neuromorphicsystems/astrometry) wheel — the same astrometry.net engine compiled as a Python extension |
| `image2pnm` + netpbm image conversion (interpreter-coupled) | `astropy.io.fits` / Pillow |
| simplexy source extraction | `sep` (Source-Extractor), matched-filter detection + spatial uniformization |
| index files hand-managed in a system directory via wget scripts | downloaded programmatically into a per-user cache on first run, scales recorded by construction |
| `plot-constellations` (compiled-in catalogs + `hd.fits`) | Pillow drawing through the solved WCS; catalogs auto-downloaded: d3-celestial constellation lines, IAU-CSN star names, OpenNGC, astrometry.net `hd.fits` |

First run downloads ~45 MB into `~/.cache/psa/` (override with
`--cache-dir` or pre-warm with `--prefetch`); afterwards everything runs
fully offline (verified with all HTTP traffic blocked).

Default index scales suit wide-field camera-lens imaging: series
4100 + 4200, scales 11–19 (quad sizes ≈ 1.4°–33°, fields from a few degrees
up to all-sky). For narrower fields download deeper scales, e.g.
`--series 4200 --scales 7-12`, or Gaia-based `--series 5200 --scales 0-6`
for telescope-scale FOVs.

## Outputs

Written to `"<name> Solved/"`:

- `<name>.wcs` — FITS WCS header; **readable by astrometry.net tools**
  (`wcsinfo`, `plot-constellations -w`, etc.)
- `annotations.png` — annotated image. By default annotations are composited
  over the (stretched) source image; `--transparent` produces an
  overlay-only RGBA image like plot-constellations.
- `solution.json` — machine-readable summary (center, scale, rotation,
  field size, matched index, annotation counts).

## Annotation layers

- constellation stick figures + names (`--no-constellations` to disable)
- IAU named bright stars, labeled `Name / α Con` (`--no-bright`)
- Bayer/Flamsteed-designated stars without IAU names (e.g. `γ Cas`) down to
  `--bright-mag` (default 4.0; 0 disables)
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

Without metadata (messaging apps typically strip EXIF on re-encode), a 60°+
field can take many minutes blind — same as solve-field. Supply hints:

```
./psa.py img.jpg --scale-low 40 --scale-high 90     # arcsec/pixel bounds
./psa.py img.jpg --ra 332 --dec 47 --radius 25      # position hint
```

`--no-auto-hint` disables metadata hints (true blind solve).

## Validation (2026-06-11)

Tested against archived `solve-field`/`plot-constellations` solutions of the
same images (DSLR stacks, Siril FITS output, and phone frames), comparing
with the original `wcsinfo`:

| image | field | solve time | vs reference solution |
|---|---|---|---|
| Cassiopeia field, DSLR stack JPEG (12455×8250) | 10.4°×6.9° | 22 s blind | center Δ≈5″, scale Δ0.08%, corner mapping ≤34″ |
| Lyra/M57 field, DSLR stack JPEG | 10.5°×7.0° | 13 s blind | center Δ≈5″ |
| Pixel Night Sight frame (EXIF intact) | 54°×73° | **1.8 s** | center Δ≈5 px, scale Δ0.6% |
| Siril-stacked FITS (`FOCALLEN`/`XPIXSZ` hint) | 15.7°×10.5° | **1.6 s** | center Δ≈5″ vs same-night solve |
| Phone frame, ~40% tree occlusion, EXIF stripped | 73°×55° | 8.5 min with `--scale-low/high` | solved; reference pipeline **failed** after 29 min CPU on the same pixels |

Catalog spot checks (`--check`): Vega → HD 172167 (3″), Sirius → HD 48915,
M 57 position 1.5″, 451 IAU names, 2954 Bayer/Flamsteed stars, 88
constellations.

**[docs/HOW-IT-WORKS.md](docs/HOW-IT-WORKS.md)** walks the full pipeline
process-by-process — algorithms, formulas, catalogs, and the
star-extraction failure modes — with annotated examples from
`docs/images/` and sample artifacts in `docs/samples/`.

## Platform notes

- Linux/macOS wheels for the `astrometry` engine; no native Windows (use
  WSL). On Raspberry Pi (aarch64) it builds from source — needs a C
  compiler.
- The original image orientation is used (EXIF rotation NOT applied), so
  pixel coordinates in the `.wcs` match the raw file, like solve-field.
- Rotation is reported in `wcsinfo`'s convention (east of north).

## Suggested alias

```bash
alias psa='/path/to/plate-solve-annotate/psa.py'
```

## License

[MIT](LICENSE.md). Runtime dependencies and runtime-downloaded data are not
vendored here and carry their own licenses — including the **GPL-3.0**
solver wheel, which matters if you redistribute bundled artifacts. See
[NOTICE.md](NOTICE.md) for the full chain and data attributions.
