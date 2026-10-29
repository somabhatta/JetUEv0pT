#!/usr/bin/env python3
"""Write the HTCondor DAG of the Phase1/Phase2 analysis (does not submit it).

Graph:
  GATE (ana_gate.sh: waits until every sample has its 500 tree files)
    -> per (cfg, sample, tree variant, acceptance) one Phase1 chain:
         [ini slices -> hadd]                         (cfgs with with_ini only)
         mean slices -> hadd -> corr slices -> hadd   (parts removed after hadd)
    -> Phase2 per comparison triple and acceptance (needs the three merged corr)

cfg = <ref>_<cent>_jt<n> (see ana_config.sh). Matrix (author 2026-10-28: no
systematic variations for now; 40 bootstrap replicas everywhere):
  refnom_call_jt0                       nom + all embedded UE variants (+ ini QA)
  refnom_c0-40_jt0, refnom_c40-80_jt0   nom, embNp, embFlow

usage:  make_dag_analysis.py <dagdir>
submit: condor_submit_dag <dagdir>/jetue_analysis.dag
The job scripts are taken from the directory of this file (PIPELINE_DIR).
"""
import os
import sys

PIPELINE_DIR = os.path.dirname(os.path.abspath(__file__))
dagdir = os.path.abspath(sys.argv[1])
os.makedirs(dagdir + "/log", exist_ok=True)

SAMPLES = ["pp_vac", "OO_vac", "OO_NR", "OO_RC", "py_pp"]
EMBEDDINGS = ["embMB", "embC", "embNp", "embFlow", "embMBmatch", "embCmatch",
              "embAng", "embAngC", "embHij"]
NSLICE = 10          # Phase1 slices per chain (every NSLICE-th tree file)
ACCEPTANCES = (0, 1)  # V0_IACC: 0 |eta|<0.9, 1 |eta|<2.5


def tree(sample, base):
    """tree variant of a sample for base 'nom' or an embedding name
    (the recoil sample OO_RC carries the cs10 recoil subtraction)"""
    if sample == "OO_RC":
        return "cs10" if base == "nom" else "cs10_" + base
    return base


# cfg -> (tree bases, run the ini QA stage)
CFGS = {"refnom_call_jt0": (["nom"] + EMBEDDINGS, True)}
for cent in ["c0-40", "c40-80"]:
    CFGS[f"refnom_{cent}_jt0"] = (["nom", "embNp", "embFlow"], False)

MEM_PHASE1 = "1.5 GB"   # 40 replicas
MEM_HADD = "2 GB"
MEM_PHASE2 = "4 GB"


def write_sub(name, exe):
    with open(f"{dagdir}/{name}.sub", "w") as f:
        f.write(f"""executable     = {PIPELINE_DIR}/{exe}
arguments      = $(args)
environment    = "PIPELINE_DIR={PIPELINE_DIR}"
output         = {dagdir}/log/$(node).out
error          = {dagdir}/log/$(node).err
log            = {dagdir}/log/condor.log
request_cpus   = 1
request_memory = $(mem)
request_disk   = 4 GB
+UseOS         = "el9"
+JobCategory   = "long"
priority       = 50
getenv         = False
queue
""")


write_sub("gate", "ana_gate.sh")
write_sub("p1", "ana_phase1_job.sh")
write_sub("hadd", "ana_hadd_job.sh")
write_sub("p2", "ana_phase2_job.sh")

lines = []


def job(node, sub, args, mem):
    lines.append(f"JOB {node} {dagdir}/{sub}.sub")
    lines.append(f'VARS {node} node="{node}" args="{args}" mem="{mem}"')


job("GATE", "gate", "", "0.5 GB")

# ---- Phase1 chains ----------------------------------------------------------
corr_merge_node = {}   # (cfg, sample, tree, iacc) -> node name of the merged corr
for cfg, (bases, with_ini) in CFGS.items():
    for sample in SAMPLES:
        trees = [tree(sample, b) for b in bases]
        if cfg == "refnom_call_jt0" and sample == "OO_RC":
            trees += ["nosub", "cs05"]          # recoil-subtraction variants
        for tr in trees:
            for iacc in ACCEPTANCES:
                chain = f"{cfg}_{sample}_{tr}_a{iacc}".replace("-", "")
                stages = (["ini"] if with_ini else []) + ["mean", "corr"]
                slices, merge = {}, {}
                for stage in stages:
                    slices[stage] = []
                    for i in range(NSLICE):
                        node = f"P1{stage.upper()}_{chain}_{i}"
                        job(node, "p1", f"{stage} {sample} {tr} {iacc} {cfg} {i} {NSLICE}", MEM_PHASE1)
                        slices[stage].append(node)
                    merge[stage] = f"HADD{stage.upper()}_{chain}"
                    job(merge[stage], "hadd", f"{stage} {sample} {tr} {iacc} {cfg} {NSLICE}", MEM_HADD)
                    lines.append(f"PARENT {' '.join(slices[stage])} CHILD {merge[stage]}")
                first = slices["mean"] + (slices["ini"] if with_ini else [])
                lines.append(f"PARENT GATE CHILD {' '.join(first)}")
                lines.append(f"PARENT {merge['mean']} CHILD {' '.join(slices['corr'])}")   # corr needs the merged mean
                corr_merge_node[(cfg, sample, tr, iacc)] = merge["corr"]

# ---- Phase2 triples ---------------------------------------------------------
# first entry = reference of the stored differences; every tag's per-replica
# v0 (v0boot*) is stored, so other differences can be formed offline with
# correlated errors.
for cfg, (bases, with_ini) in CFGS.items():
    triples = {}
    for b in bases:
        label = "main" if b == "nom" else b
        triples[f"OO_{label}"] = [("OO_vac", tree("OO_vac", b)), ("OO_NR", tree("OO_NR", b)),
                                  ("OO_RC", tree("OO_RC", b))]
        triples[f"pp_{label}"] = [("pp_vac", b), ("py_pp", b), ("OO_vac", b)]
    if cfg == "refnom_call_jt0":
        triples["OO_recoilsub"] = [("OO_RC", "cs10"), ("OO_RC", "nosub"), ("OO_RC", "cs05")]
    for name, members in triples.items():
        for iacc in ACCEPTANCES:
            node = f"P2_{cfg}_{name}_a{iacc}".replace("-", "")
            tags = ",".join(f"{s}_{t}" for s, t in members)
            job(node, "p2", f"{iacc} {cfg} {tags} {name}", MEM_PHASE2)
            parents = " ".join(corr_merge_node[(cfg, s, t, iacc)] for s, t in members)
            lines.append(f"PARENT {parents} CHILD {node}")

lines.append("RETRY ALL_NODES 2")
dagfile = f"{dagdir}/jetue_analysis.dag"
with open(dagfile, "w") as f:
    f.write("\n".join(lines) + "\n")
njobs = sum(1 for l in lines if l.startswith("JOB "))
print(f"wrote {dagfile}: {njobs} jobs, {len(corr_merge_node)} Phase1 chains")
print(f"submit with: condor_submit_dag {dagfile}")
