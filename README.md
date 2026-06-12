# plate-solve-annotate

<img width="1023" height="184" alt="image" src="https://github.com/user-attachments/assets/78d612be-d423-42fa-862e-462050d1696e" />

Self-contained, portable plate solving + sky annotation for astrophotography
in a single Python file, built on the
[astrometry.net](http://astrometry.net) solving engine
([Lang et al. 2010](https://arxiv.org/abs/0910.2233)). Blind-solves an
image against the astrometry.net index files, writes a standard WCS
solution, and renders constellation figures, star names, NGC/IC/Messier
objects, and optional Henry Draper labels over the image — with no system
astrometry.net installation.

<img width="1023" height="184" alt="annotations" src="https://github.com/user-attachments/assets/9dfa910f-980d-4fb9-8337-452058da9f69" />

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

## Background: the setup this replaces

This project supersedes a long-serving local pipeline: a shell alias
wrapping a distribution package of
[astrometry.net](http://astrometry.net) 0.97. The wrapper, essentially
verbatim:

```bash
#!/bin/bash
# predecessor: shell wrapper around the system astrometry.net install
out="${1%.*} Solved"
solve-field "$1" --downsample 2 --objs 1000 --tag-all -D "$out"
plot-constellations -w "$out/$(basename "${1%.*}").wcs" -N -C -B \
    -D -d /usr/share/astrometry/data/hd.fits -o "$out/annotations.png"
#   (the HD catalog arguments were appended only when an --hd flag was set)
```

The surrounding installation, specifically:

- **Engine config** in `/etc/astrometry.cfg`: `cpulimit 600` (a CPU budget
  per field), `add_path` entries for `/usr/share/astrometry/data/41XX` and
  `42XX`, and `autoindex`.
- **Index data** hand-fetched by `wget` scripts into that system directory:
  series 4100 scales 4111–4119 (Tycho-2, ~21 MB) and series 4200 scales
  4211–4219 (2MASS, ~17 MB) — wide-field skymarks only, quads ≈ 85′–2000′ —
  plus the 4.5 MB `hd.fits`, with mixed root/user file ownership and the
  scale choices recorded nowhere but the scripts themselves.
- **Invocation always fully blind.** The wrapper hard-coded
  `--downsample 2 --objs 1000 --tag-all` and exposed no pathway to
  solve-field's own hinting options (`--scale-low/--scale-high`,
  `--ra/--dec/--radius`), so every solve searched the entire scale ladder
  and the whole sky regardless of what was known about the image.

What that specific arrangement meant in practice, measured on real frames:

- On well-sampled ~10° DSLR fields it worked reliably — the engine plus
  Tycho-2/2MASS indexes are excellent in that regime.
- On very wide (50°+) phone frames, always-blind invocation made success
  marginal: on one 73° frame with heavy foreground occlusion, the identical
  pixels ran ~29 minutes of CPU before stopping at
  `Total CPU time limit reached` without a solution. **That is not an
  engine deficiency** — given a two-flag scale hint, the same engine solves
  that frame (and with intact EXIF metadata, this project derives such
  hints automatically; the same frame class solves in seconds).
- A routine system Python upgrade (3.13 → 3.14) orphaned `image2pnm`, the
  interpreter-coupled helper `solve-field` invokes for input conversion and
  file-type sniffing. Every solve — *including FITS input* — began failing
  with `ModuleNotFoundError`, while the compiled binaries were untouched
  and the package manager reported nothing wrong.
- Reproducing the data directory on a second machine was archaeology, and
  nothing validated or updated it afterwards.

None of this reflects on astrometry.net itself: the solver, the index
files, and even the Henry Draper kd-tree carry over into this project
unchanged. What changed is the packaging and the invocation around them —
hint plumbing surfaced (and automated from image metadata), index scales
chosen and recorded programmatically, and the interpreter coupling and
native-library chain removed:

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

## Acknowledgements & references

This tool is an orchestration layer; the heavy lifting is done by
astrometry.net and the open catalogs it draws on.

- **[astrometry.net](http://astrometry.net)** — Dustin Lang, David W. Hogg,
  Keir Mierle, Michael Blanton & Sam Roweis. The blind solver at the heart
  of this script, the [index files](https://data.astrometry.net), and the
  Henry Draper kd-tree all come from this project, which also operates the
  free hosted solver at [nova.astrometry.net](https://nova.astrometry.net).
  If results from this tool contribute to academic work, please cite:

  > Lang, D., Hogg, D. W., Mierle, K., Blanton, M., & Roweis, S. 2010,
  > *Astrometry.net: Blind astrometric calibration of arbitrary
  > astronomical images*, AJ 139, 1782 —
  > [doi:10.1088/0004-6256/139/5/1782](https://doi.org/10.1088/0004-6256/139/5/1782)
  > · [arXiv:0910.2233](https://arxiv.org/abs/0910.2233)

  <details><summary>BibTeX</summary>

  ```bibtex
  @article{Lang2010Astrometry,
    author  = {Lang, Dustin and Hogg, David W. and Mierle, Keir and
               Blanton, Michael and Roweis, Sam},
    title   = {Astrometry.net: Blind Astrometric Calibration of
               Arbitrary Astronomical Images},
    journal = {The Astronomical Journal},
    volume  = {139},
    number  = {5},
    pages   = {1782--1800},
    year    = {2010},
    doi     = {10.1088/0004-6256/139/5/1782},
    eprint  = {0910.2233},
    archivePrefix = {arXiv}
  }
  ```
  </details>

- **[`astrometry` wheel](https://github.com/neuromorphicsystems/astrometry)**
  (International Centre for Neuromorphic Systems, Western Sydney
  University) — packages the astrometry.net engine for pip and provides
  programmatic index downloads, which is what makes a self-bootstrapping
  single-file tool possible.
- **[SEP](https://sep.readthedocs.io)** (Kyle Barbary,
  [JOSS 1(6), 58](https://doi.org/10.21105/joss.00058)) — Source Extractor
  as a library, after **SExtractor** (Bertin & Arnouts,
  [A&AS 117, 393](https://ui.adsabs.harvard.edu/abs/1996A%26AS..117..393B)).
- **[Astropy](https://www.astropy.org)** — FITS I/O and all WCS/SIP
  mathematics ([citation info](https://www.astropy.org/acknowledging.html)).
- **Catalog data** — [d3-celestial](https://github.com/ofrohn/d3-celestial)
  (Olaf Frohn): constellation figures and star designations;
  [IAU Working Group on Star Names](https://www.iau.org/science/scientific_bodies/working_groups/280/):
  official star names; [OpenNGC](https://github.com/mattiaverga/OpenNGC)
  (Mattia Verga): NGC/IC/Messier objects; index files built from
  **Tycho-2** ([Høg et al. 2000](https://ui.adsabs.harvard.edu/abs/2000A%26A...355L..27H))
  and **2MASS** ([Skrutskie et al. 2006](https://ui.adsabs.harvard.edu/abs/2006AJ....131.1163S)).

Per-stage algorithm references live in
[docs/HOW-IT-WORKS.md](docs/HOW-IT-WORKS.md#12-references).

## License

[MIT](LICENSE.md). Runtime dependencies and runtime-downloaded data are not
vendored here and carry their own licenses — including the **GPL-3.0**
solver wheel, which matters if you redistribute bundled artifacts. See
[NOTICE.md](NOTICE.md) for the full chain and data attributions.
