# Third-party licenses and data attribution

The code and documentation in this repository are MIT-licensed (see
[LICENSE.md](LICENSE.md)). This notice covers everything the script *uses*
but does not contain.

## Runtime dependencies

None of these are vendored in this repository — uv/pip installs them
separately on first run. Each carries its own license:

| component | role | license |
|---|---|---|
| [`astrometry`](https://github.com/neuromorphicsystems/astrometry) | solver engine (compiled astrometry.net C code) | **GPL-3.0** |
| [`sep`](https://sep.readthedocs.io) | source extraction (SExtractor-derived) | LGPL-3.0+ |
| [`astropy`](https://www.astropy.org) | FITS / WCS | BSD-3-Clause |
| `numpy` | arrays | BSD-3-Clause |
| [`Pillow`](https://pillow.readthedocs.io) | imaging / drawing | MIT-CMU |

Using or modifying this script locally imposes nothing beyond MIT. If you
**redistribute a combined artifact that bundles the `astrometry` wheel** —
a frozen executable, container image, or vendored environment — that
distribution as a whole must satisfy GPL-3.0. The MIT license on this
repository's code is GPL-compatible, so that is permitted, but GPL
compliance becomes the redistributor's obligation.

## Runtime-downloaded data

Fetched into the local cache (`~/.cache/psa/` by default) on first use;
not part of this repository:

- astrometry.net **index files** (series 4100/4200) and **`hd.fits`**
  (Henry Draper catalog kd-tree) — <https://data.astrometry.net>
- **d3-celestial** constellation lines and star designation data,
  © Olaf Frohn, BSD-3-Clause — <https://github.com/ofrohn/d3-celestial>
- **IAU-CSN** star name list, CC-BY-4.0 — attribution: IAU Working Group
  on Star Names (WGSN) —
  <https://www.pas.rochester.edu/~emamajek/WGSN/IAU-CSN.txt>
- **OpenNGC**, CC-BY-SA-4.0 — attribution: Mattia Verga —
  <https://github.com/mattiaverga/OpenNGC>

## Example images

The images under `docs/` are the author's own photographs, annotated with
data derived from the catalogs above.
