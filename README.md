# JetUEv0pT — jet–underlying-event v0(pT) in JEWEL and PYTHIA

Code for the model study of the jet–UE transverse-momentum correlation

    v0^{jet-UE}(b) = Cov(f_b, [pT]_ref) / ( <f_b> sqrt(C2,ref) ),

where `f_b = n_b / N` is the fraction of the leading jet's in-cone charged particles in probe-pT bin `b`, and `[pT]_ref` is the mean pT of the reference particles outside the jet: 0.2 < pT < 3 GeV/c, excluding the jet cone and the two perpendicular cones. The samples are pp and O+O at √s_NN = 5.36 TeV:

- **JEWEL 2.4.0:** vacuum, medium without recoils, medium with recoils;
- **PYTHIA 8:** vacuum;
- **added UE models,** embedded into each of these samples.

The chain has three steps. Each generation job runs the generator and the per-event analyzer together and keeps only ROOT trees. The Phase1/Phase2 framework turns the trees into v0(pT) with bootstrap errors. The plotting macros make the figures.

```
generation/   JEWEL centrality patch + build script, medium config, jetv0ana (per-event analyzer → trees)
analysis/     Phase1 (ini / mean / corr moments vs fine Nch) and Phase2 (v0(pT), bootstrap, differences)
plotting/     figures from the Phase2 results
pipeline/     HTCondor job scripts and DAG writers (generation and analysis)
validation/   bit-level comparison macros used to check that code changes leave outputs unchanged
paper/        manuscript source and its figures; note on the relation to energy-energy correlators
```

`generation/README_generation.md` and `analysis/README_analysis.md` give the details of each part. This file explains how to run the full chain.

---

## Requirements

All steps run on the Nikhef stbc cluster (el9, HTCondor) using LCG views from CVMFS.

| step | environment |
|---|---|
| JEWEL build and run | `source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_107 x86_64-el9-gcc14-opt` (gcc 14; **gcc 16 miscompiles JEWEL**) |
| jetv0ana build and run (incl. PYTHIA 8, FastJet, fjcontrib ConstituentSubtractor) | `source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_110 x86_64-el9-gcc16-opt` |
| analysis and plotting (ROOT, ACLiC) | `source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_107 x86_64-el9-gcc14-opt` |

The paths in the job scripts point to the production area:
- `B=/data/alice/sbhatta/JewelV0pT`: binaries, logs, results;
- `D=/dcache/alice/sbhatta/JewelV0pT`: trees and Phase1 outputs.

To run elsewhere, change them at the top of `pipeline/gen_job.sh` and `pipeline/ana_config.sh`.

---

## 1. Generation: events → trees

### 1a. Build
```bash
bash generation/jewel/build_jewel.sh               # patched JEWEL 2.4.0 → jewel-2.4.0-vac, jewel-2.4.0-simple
source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_110 x86_64-el9-gcc16-opt
make -C generation                                 # → generation/jetv0ana
```
- **Centrality patch:** JEWEL 2.4.0 writes zeros into the HepMC heavy-ion (`H`) record. The patch writes the centrality that JEWEL sampled for the medium of that event (fraction 0–1; −1 in vacuum). `jetv0ana` stores it in percent in the branch `centr`.
- **Tables:** before running in parallel, make JEWEL's tables once per PDF set (see `generation/README_generation.md`).

### 1b. What one generation job does
`pipeline/gen_job.sh <sample> <jobid> <nevents>` runs four steps:
1. Generate `<nevents>` events of one sample:

   | sample | generator |
   |---|---|
   | `pp_vac` | JEWEL pp vacuum |
   | `OO_vac` | JEWEL O+O vacuum |
   | `OO_NR` | JEWEL O+O medium, recoils not kept |
   | `OO_RC` | JEWEL O+O medium, recoils kept |
   | `py_pp` | PYTHIA 8 pp (Monash, MPI on) |

   All samples use a p̂T^4.5 bias above 5 GeV/c with compensating event weights.
2. Run `jetv0ana` on the events. It writes one tree `t_<variant>` per particle list, all from the same generated events:

   | variant | particle list |
   |---|---|
   | `nom` | generator particles only. For `OO_RC` this is `cs10`, `cs05` or `nosub`: the recoil thermal momenta removed with constituent subtraction (ΔR_max = 1.0 or 0.5) or not removed |
   | `embMB`, `embC` | + independent thermal pions, dN/dη = 42 (O+O MB) or 135 (O+O 0–5%), ⟨pT⟩ = 0.7 GeV/c |
   | `embNp` | + thermal pions with dN/dη scaled event by event with N_part(centrality) (optical Glauber ¹⁶O, normalized to the CMS 0–5% and MB values) |
   | `embFlow` | + thermal pions whose multiplicity **and** ⟨pT⟩ follow the centrality JEWEL used for that event's medium, plus a correlated per-event fluctuation |
   | `embMBmatch`, `embCmatch` | as embMB/embC, with the UE ⟨pT⟩ in the reference window matched to the jet-correlated reference particles (0.822 GeV/c) |
   | `embMBpt060`, `embMBpt080` | as embMB, with reference-window ⟨pT⟩ 0.60 or 0.80 GeV/c (scan) |
   | `embAng`, `embAngC`, `embHij` | + charged particles of real Angantyr or HIJING events in a given multiplicity window (30–60, 110–160, 30–60) |

   For `OO_RC`, every embedded variant is prefixed `cs10_`.
3. Copy the tree file to `$D/trees_v2/<sample>/` and check the copy.
4. Delete the HepMC file.

Each tree entry is one event and has, for every (acceptance, R) slot:
- the leading anti-kT jet: raw pT, ρA-corrected pT, constituent-subtracted pT, and the pT of the matching jet in the same event without added UE;
- the in-cone, perpendicular-cone and constituent-subtracted particle counts per probe-pT bin;
- the reference sums N, Σ pT, Σ pT² for six reference definitions;
- event weight, centrality and multiplicities.

The header of `generation/jetv0ana.cc` lists every option and branch.

### 1c. Production
```bash
python3 pipeline/make_dag_generation.py /data/alice/sbhatta/JewelV0pT/dag_gen       # 5 samples × 500 jobs × 20k events
condor_submit_dag /data/alice/sbhatta/JewelV0pT/dag_gen/jetv0_v2.dag
```
About 3,800 CPU-h and 0.75 TB of trees. A failed job is retried twice with a new seed, and JEWEL runs are killed after 3 h. Resubmitting the same `.dag` picks up the rescue file.

---

## 2. Analysis: trees → v0(pT)

### 2a. Build
```bash
source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_107 x86_64-el9-gcc14-opt
cd analysis
root -l -b -q -e 'gSystem->AddIncludePath("-DV0_NBOOT=40")' -e '.L Phase1.C+O' -e '.L Phase2.C+O'
```
`V0_NBOOT` sets the number of Poisson bootstrap replicas; production used 40.

### 2b. What is computed
- **Phase1, per sample × tree variant × acceptance × configuration.** Phase1 runs in three stages: `ini` (QA), `mean` (⟨f_b⟩, ⟨[pT]_ref⟩) and `corr` (⟨f_b [pT]_ref⟩, two-particle C2,ref). Every moment is accumulated in fine bins of mid-rapidity multiplicity, for nominal weights and for each bootstrap replica in the same pass. Replica weights are seeded by the global event index, so differences between tree variants of the same events have correlated errors.
- **Phase2, per comparison triple.** Phase2 combines the multiplicity bins into v0(pT) as an event-weighted average of ratios, with errors from the replica spread. It writes:
  - the differences between the samples of the triple;
  - per-replica values (`v0boot*`) for any further difference.

An **analysis configuration** `cfg = <reference>_<centrality>_jt<n>` selects how the trees are analysed:

| part | values |
|---|---|
| reference | `refnom` nominal; `refgap` |η − η_jet| > R + 0.5; `refaw` nominal minus the away side |Δφ − π| < 1; `refnear`, `reftrans`, `refaway` Δφ slices of the nominal reference |
| centrality | `call` no cut; `c<lo>-<hi>` cut on `centr`. Events without a centrality (vacuum without UE) pass every class |
| jet-pT definition of the jet classes | `jt0` ρA-corrected (nominal); `jt1` constituent-subtracted; `jt2` truth (no added UE) |

Jet classes: 10–20, 20–40, 40–60, 60–80, 80–120 GeV/c.

Probes, each with its own matched reference mean:
- `f` raw;
- `f_sub = (n_b − b_b)/(N − B)`, with the perpendicular-cone background subtracted;
- `f_cs`, from the constituent-subtracted jet.

### 2c. Production
```bash
python3 pipeline/make_dag_analysis.py /data/alice/sbhatta/JewelV0pT/dag_ana
condor_submit_dag /data/alice/sbhatta/JewelV0pT/dag_ana/jetue_analysis.dag
```
The DAG starts with a gate job that waits until every sample has its 500 tree files. Then each Phase1 chain runs: 10 mean slices → hadd → 10 corr slices → hadd. The partial files are removed after merging. Last, Phase2 runs once per triple. The matrix is set in `make_dag_analysis.py`. Production:
- `refnom_call_jt0`: all tree variants;
- `refnom_c0-40_jt0` and `refnom_c40-80_jt0`: `nom`, `embNp`, `embFlow`.

That is 4,819 jobs, about 500 CPU-h. Outputs:
- Phase1: `$D/phase1_v3/<cfg>/acc<i>/<sample>_<variant>_{ini,mean,corr}.root`;
- Phase2: `$B/results_v3/<cfg>/acc<i>/v0pt_<triple>.root`.

`acc0` is |η_trk| < 0.9, used for R ≤ 0.4; `acc1` is |η_trk| < 2.5, used for R = 0.8.

To run one slice by hand, see `analysis/README_analysis.md`.

---

## 3. Plotting: v0(pT) → figures

```bash
source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_107 x86_64-el9-gcc14-opt
cd plotting
root -l -b -q 'PlotDeltaV0vsR.C+("/data/alice/sbhatta/JewelV0pT/results_v3","plots")'
root -l -b -q 'PlotBackgroundBaseline.C+("/data/alice/sbhatta/JewelV0pT/results_v3","plots")'
```

| macro | output |
|---|---|
| `PlotDeltaV0vsR.C` | Δv0 vs R with the probe-pT bins side by side. Columns: no recoils − vacuum, recoils − vacuum. Rows: no added UE / thermal UE / flow-correlated UE. One page per event class (MB, 0–40%, 40–80%) and jet class |
| `PlotBackgroundBaseline.C` | soft (1–6 GeV/c) and hard (12–100 GeV/c) range summaries vs R of the medium modification S(U) = v0[medium+U] − v0[vacuum+U] and of the UE-only change F_g(U) = v0[g+U] − v0[g] for JEWEL vacuum and PYTHIA. Pages: thermal UE with matched-⟨pT⟩ closure, flow-correlated UE, Angantyr/HIJING UE. Also writes a text table of all numbers |
| `PaperFigs.C`, `AppendixQA.C` | the figures of the current paper draft. They read the Phase2 outputs of the first production (`results_coarse/<ref>/acc<i>/`) and, for `qa_jets` and `qa_spectra`, its trees (now deleted; these two plots are skipped when the trees are missing). Axes cover all probe bins (0.2–100 GeV/c) and y-ranges follow the data |
| `JetUEResults.h` | shared helpers: finding a sample's curve in the Phase2 files, replica-wise differences, range summaries |
| `PlotEventDisplay.C` | per-event displays of the jet finding, background, constituent subtraction and recoil subtraction (see [Event displays](#event-displays)) |

All differences use replica-wise bootstrap errors. These include the correlation between samples built from the same events.

---

## Event displays

To show what the analyzer does to individual events, `jetv0ana` can write the intermediate objects of its first events to a separate file. The option is off by default and does not change the main output (checked bit by bit, see below).

```bash
source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_110 x86_64-el9-gcc16-opt
jetv0ana -pythia 200 -seed S <production options> -out py.root -display py_display.root            # PYTHIA
jetv0ana -hepmc OO_RC.hepmc -recoil -seed S <options> -out rc.root -display rc_display.root \
         -displayEvents 20                                                                        # JEWEL recoils
```

| option | meaning |
|---|---|
| `-display FILE` | write the display file |
| `-displayEvents N` | number of events written (default 5) |
| `-displayVariants A,B,...` | variants written (default: the no-UE variant `nom` or `cs10`, and `embMB` or `cs10_embMB`) |

The file holds, for |η_trk| < 0.9 with R = 0.4 and |η_trk| < 2.5 with R = 0.8:
- the tracks entering the clustering, with their origin (generator, added UE, reduced by the recoil subtraction), their role (leading-jet constituent, perpendicular cone, reference, refaw, η-gap reference, CS-jet constituent) and their pT after the constituent subtraction;
- all anti-kT jets above 3 GeV/c with their active-area ghosts, and which jet was selected;
- ρ, the per-probe-bin n_b, b_b, n_b − b_b, f_b and f_sub;
- the CS ghost grid and the raw, ρA-corrected, CS and truth jet pT;
- for `-recoil`: the final-state particles and thermal scattering centres before and after the thermal-momentum subtraction (ΔR_max = 1.0 and 0.5).

The values are recorded inside the production functions and are not recomputed. The jet ghosts are read from FastJet's cluster sequence, so no extra random numbers are drawn. The header of `generation/jetv0ana.cc` lists every branch.

Plot one event:
```bash
source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_107 x86_64-el9-gcc14-opt
cd plotting
root -l -b -q 'PlotEventDisplay.C+("py_display.root","plots","PYTHIA 8 pp #sqrt{s} = 5.36 TeV","pythia")'
root -l -b -q 'PlotEventDisplay.C+("rc_display.root","plots","JEWEL O+O #sqrt{s_{NN}} = 5.36 TeV, medium with recoils","jewelOO_RC",-1,"cs10,cs10_embMB")'
```
The arguments are: display file, output directory, title label, file tag, event (−1: the displayed event with the highest R = 0.4 leading-jet pT) and variants. The output is one PDF per event and one PNG per page:

| page | content |
|---|---|
| (a) particles | η–φ map, marker area ∝ pT, coloured by origin |
| (b) jets, (b2) jets R = 0.8 | jets with constituents and active-area ghosts, the fiducial band \|η_jet\| < acc − R and the selected leading jet |
| (c) background | jet and perpendicular cones, track roles, η-gap and away-side regions; n_b, b_b, n_b − b_b and f_b, f_sub vs probe pT |
| (d) constituent subtraction | track pT after vs before CS, the CS ghost grid near the jet, and the jet-pT definitions |
| (e) recoils | final-state particles and thermal momenta before/after the subtraction, pT spectra and bookkeeping |

`plotting/examples/` holds one PYTHIA event (`nom`, `embMB`) and one JEWEL O+O recoil event (`cs10`, `cs10_embMB`).

**Check (2026-10-29):** without `-display`, the modified `jetv0ana` gives output bit-identical to the version before the change. With `-display`, the main output is bit-identical to the run without it. This was checked on PYTHIA (200 events, production options; also with the display on for all 200 events and 7 variants) and on JEWEL O+O with recoils (18 events).

---

## Checking that a code change leaves the results unchanged

```bash
root -l -b -q 'validation/compare_generation_outputs.C+("old_trees.root","new_trees.root")'
root -l -b -q 'validation/compare_analysis_outputs.C+("old_mean.root","new_mean.root")'
```
On 2026-10-29, both cleaned-up versions were compared with the binaries that made the production:
- **jetv0ana:** bit-identical over 3 test runs (PYTHIA, O+O recoils, pp vacuum) of about 200 events each;
- **Phase1/Phase2:** bit-identical over 30 output files.

## Known limitations and changes relative to the 2026-10-28 production
- **Forward multiplicities `nchFwdA` and `nchFwdC`:** they now have a 0.2 < pT < 10 GeV/c cut. The production trees have no pT cut on these branches. They are not used by the analysis.
- **η range of the thermal UE:** it is generated only in |η| < 3.5, so it does not reach the FT0A-like window (3.5 < η < 4.9). Centrality classes from forward multiplicity would need the thermal UE extended.
- **η-gap reference (`gap*`):** it removes the jet cone but not the perpendicular cones.
- **Model limits:** JEWEL recoils are free-streaming, not a hydrodynamic wake. The added UE models are independent particles, apart from `embFlow`'s centrality correlation, Angantyr and HIJING.

---

## Putting this code in a new GitHub repository (github.com/somabhatta)

Run these commands on your laptop, in this directory (`JetUEv0pT/`). They create the repository `somabhatta/JetUEv0pT`, make the first commit and push it. Pick `--private` or `--public` as you prefer. You can change the visibility later in the repository settings.

**Option A: GitHub CLI (`gh`)**
```bash
brew install gh            # once, if gh is not installed
gh auth login              # once: choose GitHub.com, HTTPS, log in with the browser
cd ~/Downloads/JewelV0pT/JetUEv0pT
git init -b main
git add .
git commit -m "Jet-UE v0(pT) model study: generation, Phase1/Phase2 analysis, plotting, pipeline"
gh repo create somabhatta/JetUEv0pT --private --source=. --remote=origin --push
```

**Option B: web page + git**
1. Go to https://github.com/new. Set Owner `somabhatta` and Repository name `JetUEv0pT`, choose Private or Public, and do **not** add a README, .gitignore or licence (this directory already has them). Click *Create repository*.
2. Then, in this directory:
```bash
cd ~/Downloads/JewelV0pT/JetUEv0pT
git init -b main
git add .
git commit -m "Jet-UE v0(pT) model study: generation, Phase1/Phase2 analysis, plotting, pipeline"
git remote add origin https://github.com/somabhatta/JetUEv0pT.git    # or git@github.com:somabhatta/JetUEv0pT.git with an SSH key
git push -u origin main
```
GitHub asks for a personal access token instead of a password over HTTPS: https://github.com/settings/tokens (classic token with the `repo` scope). Alternatively, add an SSH key at https://github.com/settings/keys and use the `git@github.com:` URL.

**Later changes**
```bash
git add -A
git commit -m "describe the change"
git push
```

**Before the first push**, check `git status`. Only source files, READMEs and the small example plots should be listed. `.gitignore` excludes ROOT files, compiled libraries, HepMC files and plot output directories. `generation/jewel/jewel-2.4.0-centrality.patch` keeps the original comment line with the production path, so that it applies to the deployed source exactly.

To share the repository with collaborators, use the repository's *Settings → Collaborators* page.
