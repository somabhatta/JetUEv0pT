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
| `fig_pp`, `fig_dv0_pt`, `fig_dv0_R`, `fig_methods` | `PaperFigs.C` | Phase2 results of the first production (`results_coarse/`) |
| `qa_reference`, `qa_refgap`, `qa_annulus`, `qa_fmean`, `qa_closure`, `qa_recoilsub` | `AppendixQA.C` | same |
| `qa_jets`, `qa_spectra` | `AppendixQA.C` | trees of the first production (since deleted; the files here are the last ones made) |
| `evd_antikt`, `evd_perpcone`, `evd_cs`, `evd_recoil` | `PlotEventDisplay.C` | `jetv0ana -display` output for one PYTHIA event (with thermal UE) and one JEWEL O+O recoil event; see `../plotting/examples/` |

The numbers quoted in the text come from the first production. The 2026-10-28 production (`results_v3/`) adds centrality classes, more UE models and replica-wise errors. `PlotDeltaV0vsR.C` and `PlotBackgroundBaseline.C` read it; the manuscript does not use it yet.
