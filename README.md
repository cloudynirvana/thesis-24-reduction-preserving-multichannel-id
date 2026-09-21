# Reduction-preserving multi-channel identifiability: which Fisher ranks survive a stiff–sloppy reduction?

**Thesis #24.** Computational research, set out in Nile University B.Sc. chapter order for handoff.

**Depends on:** [Thesis #9](https://github.com/cloudynirvana/thesis-09-ccle-metabolic-ode-identifiability) (multi-channel metabolic ODE identifiability) and [Thesis #12](https://github.com/cloudynirvana/thesis-12-stiff-sloppy-cancer-ode-reduction) (stiff–sloppy / MBAM-style reduction).

**Author:** Kelechi Emeka Ogbonna  
**Email:** kelechiogbonna300@gmail.com  
**GitHub:** https://github.com/cloudynirvana  
**Date:** 21 September 2026

After a documented stiff–sloppy or MBAM-style reduction, which multi-channel Fisher ranks from the full metabolic ODE survive on the reduced system, and which ranks are artefacts of the unreduced coordinates?

The calculations use one eight-rate metabolic generator. Glucose, pyruvate, lactate and a glutamine-linked pool share a right-hand side with a fast modifier. The documented reduction slaves the modifier and keeps the product κ = h a / b. The same channel schedules are ranked on both parameter vectors. Practical rank uses a cut of 10<sup>−3</sup> times the leading Fisher eigenvalue. Draws are a synthetic surrogate, not a CCLE file.

On this surrogate the shared practical ranks survive. The steady lactate/glucose ratio has practical rank 1 of 8 and 1 of 6. The four-channel snapshot has practical rank 5 of 8 and 5 of 6, so the multi-channel advantage of 4 survives. A lineage loading, profiled out, leaves practical rank 3 on both vectors. Time courses stay at practical rank 5. The schedule that also records the modifier has practical rank 6 of 8 and has no reduced counterpart: that rank vanishes. Time-course numerical ranks that sit in the deleted (ρ, speed) plane also vanish; they never cleared the practical cut.

This is research only. It is not a medical device, not clinical decision support, not a dose, and not a cure. No document DOI is registered.

See [DISCLAIMER.md](DISCLAIMER.md). The manuscript is [THESIS.md](THESIS.md).

## Files

| Path | Role |
| --- | --- |
| `THESIS.md` | Manuscript (Chapters 1 to 5, Vancouver citations) |
| `THESIS.pdf` | PDF built from the Markdown |
| `build_pdf.py` | Regenerates `THESIS.pdf` |
| `CITATION.cff` | Citation metadata, no document DOI |
| `DISCLAIMER.md` | Research-only boundary |
| `sim/reduction_ranks.py` | Seeded full-versus-reduced Fisher ranks (seed 20260921) |
| `sim/results.json` | Numbers cited in Chapter Four |
| `sim/figures/` | Spectra, rank bars, axis overlaps, null walk, lactate |

## Reproduce

```bash
python3 -m pip install -r sim/requirements.txt
python3 sim/reduction_ranks.py
python3 build_pdf.py
```

NumPy, SciPy and Matplotlib are required for the ranks. The PDF step also needs the `markdown` and `weasyprint` packages. Regenerating the script rewrites `sim/results.json` and `sim/figures/`.

## Cite

Ogbonna KE. Reduction-preserving multi-channel identifiability: which Fisher ranks survive a stiff–sloppy reduction? [Internet]. Thesis #24 computational research thesis. 21 September 2026 [cited YYYY Mon DD]. Available from: https://github.com/cloudynirvana/thesis-24-reduction-preserving-multichannel-id

Machine-readable fields are in `CITATION.cff`. Add a document DOI there only after one exists.

Hub index, for cataloguing only: [research-theses-hub](https://github.com/cloudynirvana/research-theses-hub).

## Licence

Text and sketch code are MIT, with attribution. Computational research only.
