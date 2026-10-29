#!/usr/bin/env python3
"""Write the HTCondor DAG of the generation step (one node per gen_job.sh call).

usage: make_dag_generation.py <dagdir> [--test] [--exe PATH]

  <dagdir>   directory for gen.sub, jetv0_v2.dag and log/ (created)
  --test     2 jobs x 200 events per sample instead of the full production
  --exe      job script (default: gen_job.sh next to this file; the production
             DAG used /data/alice/sbhatta/JewelV0pT/pipeline/gen_job_v2.sh)

Production (2026-10-28): 5 samples x 500 jobs x 20k events = 2500 jobs, 10M
events per sample; O+O medium samples are minimum bias, centrality classes are
cut offline. Every node is retried twice; gen_job.sh reseeds each retry.
Submit with:  condor_submit_dag <dagdir>/jetv0_v2.dag
The analysis (Phase1/Phase2) DAGs are written separately once the trees exist.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

# sample -> (number of jobs, events per job); see gen_job.sh for the sample definitions
SAMPLES = {
    "pp_vac": (500, 20000),   # JEWEL vacuum pp
    "OO_vac": (500, 20000),   # JEWEL vacuum O+O (nPDF)
    "OO_NR":  (500, 20000),   # JEWEL medium O+O MB, no recoils
    "OO_RC":  (500, 20000),   # JEWEL medium O+O MB, recoils kept
    "py_pp":  (500, 20000),   # PYTHIA 8 pp
}
TEST_JOBS, TEST_EVENTS = 2, 200


def main(argv):
    args = [a for a in argv[1:]]
    test = "--test" in args
    exe = os.path.join(HERE, "gen_job.sh")
    if "--exe" in args:
        exe = os.path.abspath(args[args.index("--exe") + 1])
    positional = [a for i, a in enumerate(args) if not a.startswith("--") and (i == 0 or args[i - 1] != "--exe")]
    if len(positional) != 1:
        sys.exit(__doc__)
    dagdir = os.path.abspath(positional[0])
    os.makedirs(dagdir + "/log", exist_ok=True)
    samples = {s: (TEST_JOBS, TEST_EVENTS) for s in SAMPLES} if test else SAMPLES

    # submit description shared by all nodes: 1 core, 2 GB memory, 4 GB scratch
    # (HepMC of 20k JEWEL events), el9, long queue
    with open(f"{dagdir}/gen.sub", "w") as f:
        f.write(f"""executable     = {exe}
arguments      = $(args)
output         = {dagdir}/log/$(node).out
error          = {dagdir}/log/$(node).err
log            = {dagdir}/log/condor.log
request_cpus   = 1
request_memory = 2 GB
request_disk   = 4 GB
+UseOS         = "el9"
+JobCategory   = "long"
priority       = 50
getenv         = False
queue
""")

    # one node per (sample, job id); arguments: <sample> <jobid> <nevents>
    lines = []
    for sample, (njobs, nevents) in samples.items():
        for job in range(njobs):
            node = f"GEN_{sample}_{job}"
            lines += [f"JOB {node} {dagdir}/gen.sub", f'VARS {node} node="{node}" args="{sample} {job} {nevents}"']
    lines.append("RETRY ALL_NODES 2")
    with open(f"{dagdir}/jetv0_v2.dag", "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"wrote {dagdir}/jetv0_v2.dag: {sum(n for n, _ in samples.values())} jobs ({'TEST' if test else 'FULL'})")


if __name__ == "__main__":
    main(sys.argv)
