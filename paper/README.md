# Paper sources

| file | content |
|---|---|
| `jetv0pt_model.tex`, `refs.bib` | manuscript (revtex4-2), *Probing jet-induced medium response with jet--underlying-event transverse-momentum correlations*, and its references |
| `figs/` | the figures the manuscript includes, exactly as used |
| `eec_note/eec_v0_relation.tex` | note: the mathematical relation between energy–energy correlators and v0^{jet-UE} |
| `eec_note/make_figs.py`, `eec_note/figs/` | the schematic figures of the note (matplotlib; no data) |
| `eec_note/refs_eec_note.bib` | the note's references (checked against INSPIRE) |

## Building

```bash
cd paper
pdflatex jetv0pt_model && bibtex jetv0pt_model && pdflatex jetv0pt_model && pdflatex jetv0pt_model
cd eec_note
python3 make_figs.py                     # needs numpy + matplotlib (e.g. an LCG view)
pdflatex eec_v0_relation && bibtex eec_v0_relation && pdflatex eec_v0_relation && pdflatex eec_v0_relation
```

`refs.bib` holds the manuscript's references.

## Where each figure comes from

| figure file | macro (in `../plotting/`) | input |
|---|---|---|
| `figV3_main_probebins` (Fig. 2), `figV3_main_soft_classes` (Fig. 3), `figV3_app_hard_classes`, `figV3_baseline`, rows of Table I | `PaperFigsV3.C` | Phase2 results of the 2026-10-28 production (`results_v3/`) |
| `deltaV0_vsR_refnom_{call,c0-40,c40-80}_jt0_jet40_60` (appendix) | `PlotDeltaV0vsR.C` | same |
| `fig_pp`, `fig_methods` | `PaperFigs.C` | Phase2 results of the first production (`results_coarse/`) |
| `qa_reference`, `qa_refgap`, `qa_annulus`, `qa_fmean`, `qa_closure`, `qa_recoilsub` | `AppendixQA.C` | same |
| `qa_jets`, `qa_spectra` | `AppendixQA.C` | trees of the first production (since deleted; the files here are the last ones made) |
| `evd_antikt`, `evd_perpcone`, `evd_cs`, `evd_recoil` | `PlotEventDisplay.C` | `jetv0ana -display` output for one PYTHIA event (with thermal UE) and one JEWEL O+O recoil event; see `../plotting/examples/` |

The Results section uses the 2026-10-28 production (`results_v3/`). The pp comparison with MPI off, the subtraction and reference checks, and the appendix QA figures use the first production, as stated in their captions. Table II (binned vs single-class estimator) was computed from the Phase1 moments with the diagnostic described in the text.
