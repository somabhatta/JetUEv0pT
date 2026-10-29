# Jet–UE v0(pT): Phase1 / Phase2 analysis

This directory holds the offline analysis. It turns the per-event jetv0ana trees into

    v0(pT_b) = Cov(f_b, [pT]_ref) / ( <f_b> sqrt(C2,ref) )

where `f_b = n_b / N` is the fraction of the leading jet's in-cone charged particles that fall in probe bin `b`, and `[pT]_ref` is the mean pT of the reference particles outside the jet.

| file | role |
|---|---|
| `Rule.h` | constants, tree layout, run options (env), bootstrap weights |
| `Tool.C` | TChain from a file list, progress counter, Nch-cut file IO |
| `Phase1.h`, `Phase1.C` | `Event` record and the staged `Phase1` driver |
| `Phase1_evt.C` | unpacks one tree entry: reference sums, jet-class and centrality selection, spectra |
| `Phase1_ini.C` | stage ini: QA distributions |
| `Phase1_mean.C` | stage mean: `<f_b>`, `<[pT]_ref>` vs fine Nch |
| `Phase1_corr.C` | stage corr: `<f_b [pT]_ref>`, two-particle C2 vs fine Nch |
| `Phase2.C` | combines the Nch bins into v0(pT), bootstrap errors, sample differences |

The job scripts and the DAG writer are in `../pipeline/` (`ana_*.sh`, `make_dag_analysis.py`).

## What is computed

### Phase1 (one sample × tree variant × acceptance × configuration, sliced)

Each stage runs Init → Run (event loop: `Fill_Evt` + `Fill_<stage>`) → Final.

* **Event unpacking** (`Phase1_evt.C`). For each R (0.1, 0.2, 0.4, 0.8) the code reads the reference sums `N_ref`, `S1 = Σ pT` and `S2 = Σ pT²` of the chosen reference. The leading jet goes into jet class `c` (`jet10_20 … jet80_120`) if two things hold:
  * its pT is in `[SELLO, SELHI)`, using the jet-pT definition set by `V0_JETPT`;
  * the event passes the centrality cut. Events with `centr < 0` always pass.

  For that class the code copies the in-cone spectrum `n_b`, the perpendicular-cone spectrum `p_b` and the spectrum of the constituent-subtracted jet `cs_b`.
* **Event gate** (mean and corr use the same gate): `[pT]_ref ≥ 0`, `N_ref ≥ 2` (the two-particle C2 needs it) and S2 present. Mean and corr therefore average over the same events, so `Σ_b Cov(f_b,[pT]_ref) = 0` holds exactly. The corr stage prints this check.
* **Three probes**, each with a reference mean taken over exactly its own events:

  | probe | definition | exists if | reference mean | C2 |
  |---|---|---|---|---|
  | f (`""`) | `n_b / Σn` | always | `pRef_Nch` | `pC2_Nch` |
  | f_sub (`_fsub`) | `(n_b − p_b) / Σ(n − p)` | `Σ(n − p) > 0` | `pRefFsub_Nch` | `pC2Sub_Nch` |
  | f_cs (`_fcs`) | `cs_b / Σcs` | `Σcs > 0` | `pRefFcs_Nch` | `pC2Sub_Nch` |

* **C2** (corr stage) is the two-particle dynamical variance with self-pairs removed: `C2 = [(S1 − N μ)² − (S2 − 2 μ S1 + N μ²)] / [N (N − 1)]`. Here `μ = μ(Nch)` comes from the merged nominal mean profile, at the event's Nch bin.
* **Fine Nch**: every moment is a profile vs `nchMid` (|η| < 0.8, 0.2 < pT < 5 GeV). The default binning is 100 bins on 0–500.
* **Poisson bootstrap**: mean and corr fill `NBOOT+1` copies of every profile in the same loop:
  * `k = 0` is the nominal set (event weight `w`);
  * `k ≥ 1` uses `w · Poisson(1)`, seeded by the global event index `entry + V0_EVTOFFSET` and by `k`.

  Because of this seeding, replica `k` is coherent between the mean and corr files. Re-running a job reproduces it, and `hadd` of the slices equals a single run. Tree variants of the same generated events (e.g. `nom` vs `embFlow`) get the same weights, so their differences have correlated errors.

### Phase2 (three merged samples)

For each sample, R, jet class and probe, and for each output pT bin Q:

1. `cov_j = Σ_{b∈Q} (<f_b ref>_j − <f_b>_j μ_j)` and `<f>_j = Σ_{b∈Q} <f_b>_j` in each fine Nch bin `j`. Summing over the fine bins of a coarse bin is exact.
2. **Average of ratios** over the Nch bins of the class, weighted by the event count `N_j`: `v0 = Σ_j N_j cov_j/(<f>_j σ_j) / Σ_j N_j` with `σ_j = sqrt(max(C2_j, 0))`. Also `v0refnorm = Σ_j N_j v0_j μ_j/σ_j / Σ_j N_j` (that is, v0 / v0_ref).
3. The same combination is repeated for every replica `k`, with all four inputs taken from replica `k`. A replica with any input missing is dropped (`nBoot_*` counts the ones used). The error is `CalErr` of the replica values:
   * drop non-finite values and exact zeros;
   * take the population SD;
   * remove values more than 3 SD from the mean and recompute the SD.

   The central value is always the nominal.
4. Differences to the first sample (`dv0*`): the error comes from the replica-wise differences, with the same trimmed SD but exact zeros kept.

The output pT grid is the probe grid. With `V0_PTGRID=coarse` it is 0.2, 0.5, 1, 2, 3, 6, 12, 30, 100 GeV/c (production default). Only the inclusive Nch class (`c1`) is formed. A wide class `c0` would need `INPUT/cutNch_<tag>.txt`, which the production does not have.

## Run options

| variable | meaning | default |
|---|---|---|
| `V0_TREE` | tree variant, tree `t_<V0_TREE>` | `nom` |
| `V0_IACC` | acceptance 0 = \|η\|<0.9, 1 = \|η\|<2.5 | 0 |
| `V0_IREF` | reference 0 refnom, 1 refgap, 2 refaw, 3 refnear, 4 reftrans, 5 refaway | 0 |
| `V0_CENTLO`, `V0_CENTHI` | keep `CENTLO ≤ centr < CENTHI`; `centr < 0` always passes | no cut |
| `V0_JETPT` | jet classes on 0 jetPtCorr, 1 csJetPt, 2 jetPtTruth | 0 |
| `V0_NCHBINS`, `V0_NCHMAX` | fine Nch axis | 100, 500 |
| `V0_EVTOFFSET` | global event offset of the slice (bootstrap seed) | 0 |
| `V0_MEANFILE` | merged mean read by corr | `OUTPUT/<ref>/<tag>_mean.root`, then `OUTPUT/<tag>_mean.root` |
| `V0_PTGRID` | `coarse` = coarse output grid in Phase2 | probe grid |
| `-DV0_NBOOT=n` | (compile time) number of replicas | 200 (production: 40) |

The references are refnom (outside the jet and both perpendicular cones), refgap (|Δη| > R + 0.5), refaw (refnom without |Δφ − π| < 1) and the Δφ slices of refnom: near (< π/3), trans and away (> 2π/3). All use 0.2–3 GeV/c.

## Output naming

`<tag> = <sample>_<tree>` (e.g. `OO_RC_cs10`). `<var>` is `""`, `_fsub` or `_fcs`. `<sel>` is the jet class (`jet20_40`, …). `R<i>` is the R index (0..3 = 0.1, 0.2, 0.4, 0.8). Replica `k ≥ 1` adds `_boot<k>` before `_R<i>`.

* ini `<tag>_ini_<slice>.root`: `hNch`, `hLeadHi` (centrality), `hNref_R<i>`, `hRef_R<i>`, `hNcone_R<i>` (leading-jet pT), `hNin_R<i>`, `hJetFlag_R<i>`, `hRefVsNch_R<i>`.
* mean `<tag>_mean[_<slice>].root`:
  * `pN_NchPt<var>[_boot<k>]_R<i>_<sel>` (TProfile2D)
  * `pRef_Nch…`, `pRefFsub_Nch…`, `pRefFcs_Nch…` (TProfile)
  * `REFWIN` (TNamed "ptlo pthi iref name")
* corr `<tag>_corr[_<slice>].root`:
  * `pNR_NchPt<var>[_boot<k>]_R<i>_<sel>`
  * `pC2_Nch…`, `pC2Sub_Nch…`
  * diagnostics `hCov_NchPt_*`, `Nf_NchPt_*`, `mu_Nch_*`
* Phase2 `v0pt_<name>.root`, with `<t> = <tag>_c1_R<i>_<sel>`:
  * `v0<var>_<t>`, `v0refnorm<var>_<t>`, `cov<var>_<t>`, `meanN<var>_<t>`, `v0errAnalytic<var>_<t>`
  * `v0boot<var>_<t>`, `v0refnormboot<var>_<t>` (TH2D, x = output pT bin, y bin k+1 = replica k)
  * `C2_<t>`, `refbar_<t>`, `nBoot_<t>`
  * `meanPtVsNch_<tag>_R<i>_<sel>`, `meanPtVsNchBoot_…`
  * `dv0<var>_<tag1>m<tag0>_c1_R<i>_<sel>`, `dv0refnorm<var>_…` (and the same for tag2)
  * canvases `c_jet`, `c_jet_refnorm` plus `v0pt_<name>_jet.pdf`, `v0pt_<name>_refnorm_jet.pdf`

## Build

```bash
source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_107 x86_64-el9-gcc14-opt
cd analysis
root -l -b -q -e 'gSystem->AddIncludePath("-DV0_NBOOT=40")' -e '.L Phase1.C+O' -e '.L Phase2.C+O'
```

This produces `Phase1_C.so` and `Phase2_C.so`. Leave out the `AddIncludePath` line for 200 replicas. Phase2 can be built with a larger NBOOT than the inputs: missing replicas are skipped, but `v0boot*` then gets NBOOT+1 y bins.

## Running one slice by hand

```bash
export V0_TREE=cs10 V0_IACC=0 V0_IREF=0 V0_JETPT=0      # refnom_call_jt0, |eta|<0.9
ls /dcache/alice/sbhatta/JewelV0pT/trees_v2/OO_RC/*.root | head -4 > files.list
mkdir -p out
root -l -b -q -e 'gSystem->Load("Phase1_C.so"); Phase1 p("files.list","OO_RC_cs10","out",0,-1); p.Exec_mean();'
mv out/OO_RC_cs10_mean_0.root out/OO_RC_cs10_mean.root   # with several slices: hadd them
export V0_MEANFILE=$PWD/out/OO_RC_cs10_mean.root
root -l -b -q -e 'gSystem->Load("Phase1_C.so"); Phase1 p("files.list","OO_RC_cs10","out",0,-1); p.Exec_corr();'
mv out/OO_RC_cs10_corr_0.root out/OO_RC_cs10_corr.root
```

The corr log ends with `corr CHECK: OK … max |Sum_pT cov| ≈ 1e-16`. Repeat this for the other two tags of a triple (e.g. `OO_vac_nom`, `OO_NR_nom`), then run:

```bash
V0_PTGRID=coarse root -l -b -q -e 'gSystem->Load("Phase2_C.so"); Phase2("out","OO_vac_nom,OO_NR_nom,OO_RC_cs10","v0pt_test.root");'
```

When you run several slices of one sample, give each slice its own `V0_EVTOFFSET` (the pipeline uses `slice × 10⁹`). Use the same offset for the slice's mean and corr job.

## Production (HTCondor DAG)

The scripts are in `pipeline/`:

* `ana_config.sh`: paths, the configuration-name parser `cfg_env`, `dcache_put` and the Phase1/Phase2 library locations. Override them with `PHASE1_LIB`, `PHASE2_LIB`, `TREES`, `PHASE1_OUT` and `RESROOT` (Phase2 results), e.g. for a test outside the production area.
* `ana_phase1_job.sh <stage> <sample> <tree> <iacc> <cfg> <islice> <nslice>`
* `ana_hadd_job.sh <stage> <sample> <tree> <iacc> <cfg> <nslice>`
* `ana_phase2_job.sh <iacc> <cfg> <tag0,tag1,tag2> <name>`
* `ana_gate.sh`: waits until every sample has 500 tree files.
* `make_dag_analysis.py <dagdir>`: writes `<dagdir>/jetue_analysis.dag`. It replaces `make_dag_ana_v3.py` and produces the same graph.

A configuration name looks like `<ref>_<cent>_jt<n>`, e.g. `refnom_call_jt0` or `refnom_c0-40_jt0`.

```bash
python3 pipeline/make_dag_analysis.py /data/alice/sbhatta/JewelV0pT/dag_ana
condor_submit_dag /data/alice/sbhatta/JewelV0pT/dag_ana/jetue_analysis.dag
```

The job scripts are taken from the directory of `make_dag_analysis.py`. The submit files pass that directory as `PIPELINE_DIR`.

Where the outputs land:

* Phase1 slices: `/dcache/alice/sbhatta/JewelV0pT/phase1_v3/<cfg>/acc<i>/parts/`. They are removed after a successful hadd.
* merged Phase1: `/dcache/alice/sbhatta/JewelV0pT/phase1_v3/<cfg>/acc<i>/<tag>_{ini,mean,corr}.root`
* Phase2: `/data/alice/sbhatta/JewelV0pT/results_v3/<cfg>/acc<i>/v0pt_<name>.root`, where `<name>` is `OO_main`, `OO_<emb>`, `pp_main`, `pp_<emb>` or `OO_recoilsub`.
* logs: `<dagdir>/log/<node>.{out,err}`
