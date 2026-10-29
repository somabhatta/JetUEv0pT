# Generation: JEWEL 2.4.0 / PYTHIA 8 events → jetv0ana trees

Jet–UE v0(pT) study, pp and O+O at √s_NN = 5.36 TeV. Each generation job makes
the events, runs `jetv0ana` on them in the same job, keeps only the ROOT trees and
deletes the HepMC file.

```
generation/
  jetv0ana.cc, Makefile          per-event analyzer (option/branch reference in the file header)
  jewel/
    jewel-2.4.0-centrality.patch writes the sampled centrality into the HepMC H record
    build_jewel.sh               builds the patched JEWEL (gcc 14, -O2)
    config/medium_OO_MB.dat      simple medium, O+O minimum bias (also passed to vacuum runs)
pipeline/
  gen_job.sh                     one generation job (sample, job id, events)
  make_dag_generation.py         HTCondor DAG of all generation jobs
```

All commands below are for stbc (el9). The production area is `/data/alice/sbhatta/JewelV0pT`.

## 1. Build JEWEL 2.4.0 with the centrality patch

```bash
bash generation/jewel/build_jewel.sh [jewel-2.4.0.tar.gz] [BUILD_DIR]
```

The script downloads the tarball from hepforge if you don't pass one, and checks it
against the production sha256. It then applies the patch and compiles with gfortran
14.2 from `LCG_107 x86_64-el9-gcc14-opt` at `-O2`, linking LHAPDF from the same view
with an rpath. It produces `jewel-2.4.0-vac` and `jewel-2.4.0-simple`. Rebuilt this
way, both binaries and all object files are byte-identical to the production
`jewel-2.4.0-gcc14-cent/`.

**Do not build JEWEL with gcc 16 (LCG_110).** At `-O2` it miscompiles JEWEL: medium
runs hang after a few events (NaN gluon momenta, endless loop in
MAKECASCADE/DOINSTATESCAT). Run JEWEL in the gcc 14 view as well.

The patch makes three changes. JEWEL 2.4.0 writes the `H` (heavy-ion) record with all
zeros; the patch writes `GETCENTRALITY()` into its impact-parameter slot (the 10th
number) instead. That value is the sampled centrality as a fraction (0–1) for the
simple medium and −1 in vacuum. `jetv0ana` converts it to percent (branch `centr`).
`GETCENTRALITY` is already defined in `medium-simple.f` and `medium-vac.f`, so only
`jewel-2.4.0.f` changes.

### Tables
JEWEL writes its tables on the first run when the files are missing. Production uses
`tables/g14_splitint.dat`, `g14_pdfs_<set>.dat` and `g14_xsecs_<set>.dat`
(76 MB / 0.2 MB / 11 MB), built by the gcc 14 binaries for the PDF sets 14600 (pp)
and 902750 (O+O). To make them, run one short job per PDF set *before* any parallel
running, with SPLITINTFILE/PDFFILE/XSECFILE pointing to non-existing files (see the
parameter block in `gen_job.sh`, e.g. NEVENT 50 with `jewel-2.4.0-vac`). Never reuse
tables across PDF sets or compilers.

The medium uses a simple Bjorken expansion with TAUI 0.4 fm/c and TI 0.43 GeV. TI is
the Pb default of 0.59 GeV scaled with CMS O+O 0–5% dN/dη = 135. The run is minimum
bias (CENTRMIN 0, CENTRMAX 100) with SIGMANN 6.8 fm².

## 2. Build jetv0ana

```bash
source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_110 x86_64-el9-gcc16-opt
cd generation && make
```

The analyzer is built with gcc 16 (LCG_110), the same view as the production binary.
Use that view and these flags to reproduce the existing trees bit by bit. The one exception is the forward multiplicities `nchFwdA` and `nchFwdC`. Since 2026-10-29 they count only particles with 0.2 < pT < 10 GeV/c, while the 2026-10-28 production trees have no pT cut on them.

## 3. Run one generation job by hand

```bash
bash pipeline/gen_job.sh <sample> <jobid> <nevents>     # e.g. OO_RC 7 200
```

The samples are `pp_vac`, `OO_vac`, `OO_NR`, `OO_RC` and `py_pp`. The job:

- works in `$_CONDOR_SCRATCH_DIR`, or in a temporary directory that it removes;
- sets the seed to 100000·cksum("<sample>_v2")[0:3] + jobid, plus 7919 per earlier
  attempt;
- runs JEWEL with a 3 h timeout, then `jetv0ana` with every UE option;
- copies the tree file to `/dcache/alice/sbhatta/JewelV0pT/trees_v2/<sample>/` and
  checks the copy with `cmp`;
- deletes the HepMC file.

**Warning:** a hand run writes to the production output and the attempt counters
under `prod_v2/`. For a private test, copy the script and change `B`, `OUT`,
`LOGDIR`, `ATTEMPTS` and `ANALYZER` at the top. You can also run the two steps
directly: write the `params.dat` block from the script, run
`jewel-2.4.0-simple params.dat`, then
`jetv0ana -hepmc F.hepmc [-recoil] -seed S <options> -out F.root`.

## 4. Write and submit the generation DAG

```bash
python3 pipeline/make_dag_generation.py <dagdir> [--test] [--exe /path/to/gen_job.sh]
condor_submit_dag <dagdir>/jetv0_v2.dag
```

The full DAG has 5 samples × 500 jobs × 20k events. `--test` gives 2 × 200 per sample.
Each node gets 1 CPU, 2 GB of memory and 4 GB of disk, runs on el9 in the long queue,
and has `RETRY ALL_NODES 2`. Monitor with `condor_q -dag` and
`tail <dagdir>/jetv0_v2.dag.dagman.out`. To resubmit failures, submit the same `.dag`
again; DAGMan picks up the rescue file.

## Checking a rebuilt analyzer

```bash
root -l -b -q 'validation/compare_generation_outputs.C+("old.root","new.root")'
```

The macro compares every tree entry and every histogram bin bitwise. It prints
`RESULT: bit-identical` only when nothing differs.

## Event displays

`jetv0ana ... -display display.root [-displayEvents N] [-displayVariants A,B]` writes
the intermediate objects of the first N events (default 5) to a separate file for
`plotting/PlotEventDisplay.C`. The main output does not change. See the
"Event displays" section of the top-level `README.md`.
