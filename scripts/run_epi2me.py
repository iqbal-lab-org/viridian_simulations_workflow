#!/usr/bin/env python3

import argparse
import json
import os
import subprocess
import sys

import pyfastaq

def syscall(command, cwd=None, allow_fail=True):
    completed_process = subprocess.run(
        command,
        shell=True,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE,
        universal_newlines=True,
        cwd=cwd,
    )
    if completed_process.returncode != 0:
        print("Error running this command:", command, file=sys.stderr)
        print("Return code:", completed_process.returncode, file=sys.stderr)
        print(
            "Output from stdout:", completed_process.stdout, sep="\n", file=sys.stderr
        )
        print(
            "Output from stderr:", completed_process.stderr, sep="\n", file=sys.stderr
        )
        if not allow_fail:
            raise Exception("Error in system call. Cannot continue")
    return completed_process


def bash_out_to_time_and_memory(lines):
    """Parses output file `out` of command run like this:
          /usr/bin/time -v foo &> out
    returns a dictionary of time and memory used"""
    stats = {}

    in_time_lines = False
    for line in lines:
        if not in_time_lines:
            if line.startswith("\tCommand being timed:"):
                in_time_lines = True
        elif line.startswith("\tElapsed (wall clock) time (h:mm:ss or m:ss): "):
            time_fields = line.rstrip().split()[-1].split(":")
            if not 2 <= len(time_fields) <= 3:
                raise Exception(f"Error getting time from this line: {line}")
            time_in_seconds = float(time_fields[-1]) + 60 * float(time_fields[-2])
            if len(time_fields) == 3:
                time_in_seconds += 60 * 60 * float(time_fields[0])
            stats["wall_clock_time"] = time_in_seconds
        elif line.startswith("\tUser time (seconds): "):
            stats["user_time"] = float(line.rstrip().split()[-1])
        elif line.startswith("\tSystem time (seconds): "):
            stats["system_time"] = float(line.rstrip().split()[-1])
        elif line.startswith("\tMaximum resident set size (kbytes): "):
            stats["ram"] = float(line.rstrip().split()[-1])
        elif line.startswith("\tPercent of CPU this job got: "):
            stats["percent_CPU"] = float(line.rstrip().split()[-1].rstrip("%"))

    return stats



def run_epi2me(
  outdir,
  nf_script,
  scheme_version,
  reads,
  nxf_sing_cache,
  sample_name="sample",
  work_dir_root=None,
  debug=False,
  ):
    outdir = os.path.abspath(outdir)
    nf_script = os.path.abspath(nf_script)
    reads = os.path.abspath(reads)
    nxf_sing_cache = os.path.abspath(nxf_sing_cache)
    run_script = os.path.join(outdir, "run_nextflow.sh")
    log_info = {
        "inputs": {
            "outdir": outdir,
            "nf_script": nf_script,
            "reads": reads,
            "nxf_sing_cache": nxf_sing_cache,
            "nxf_sing_cache_files": os.listdir(nxf_sing_cache),
            "debug": debug,
            "sample_name": sample_name,
        }
    }
    log_file = os.path.join(outdir, "log.json")
    nf_out = "nextflow_out"
    os.mkdir(outdir)

    if work_dir_root is not None:
        assert os.path.exists(work_dir_root)
        work_dir = os.path.join(os.path.abspath(work_dir_root), f"{sample_name}.epi2me")
        subprocess.check_output(f"rm -rf {work_dir}", shell=True)
    else:
        work_dir = os.path.join(outdir, "tmp_nf_work")

    with open(run_script, "w") as f:
        print(f"export NXF_SINGULARITY_CACHEDIR={nxf_sing_cache}", file=f)
        print(
            "/hps/nobackup/iqbal/dander/viridian_simulations_workflow/nextflow run", os.path.abspath(nf_script),
            "-ansi-log false",
            "-work-dir", work_dir,
            "-profile singularity",
            "--out_dir", nf_out,
            "--scheme_version", scheme_version,
            "--fastq", reads,
            file=f
        )

    completed_process = syscall("bash run_nextflow.sh", cwd=outdir)
    log_info["nextflow"] = {
        "stdout": completed_process.stdout.rstrip().split("\n"),
        "stderr": completed_process.stderr.rstrip().split("\n"),
    }
    log_info["nextflow"].update(bash_out_to_time_and_memory(log_info["nextflow"]["stderr"]))
    with open(log_file, "w") as f:
        json.dump(log_info, f, indent=2)

    seqs = {}
    pyfastaq.tasks.file_to_dict(os.path.join(outdir, nf_out, "all_consensus.fasta"), seqs)
    assert len(seqs) == 1
    seq = list(seqs.values())[0]
    if sample_name is not None:
        seq.id = sample_name + " " + seq.id
    with open(os.path.join(outdir, "consensus.fa"), "w") as f:
        print(seq, file=f)

    if not debug:
        subprocess.check_output(f"rm -rf {work_dir} {nf_out}/*.bam {nf_out}/*.bam.bai", shell=True, cwd=outdir)
        subprocess.check_output(f"tar cvf {nf_out}.tar {nf_out}", shell=True, cwd=outdir)
        subprocess.check_output(f"gzip -9 {nf_out}.tar", shell=True, cwd=outdir)
        subprocess.check_output(f"rm -rf {nf_out}", shell=True, cwd=outdir)





parser = argparse.ArgumentParser(
    description = "Run epi2me/wf-artic on one sample",
    usage="%(prog)s [options]",
)
required_ops_group = parser.add_argument_group("Required options")
required_ops_group.add_argument("--outdir", required=True, help="Output directory")
required_ops_group.add_argument("--scheme_version", required=True, help="Primer scheme version. See --help from main nextflow script to get a list of possible values. Probably want one of: Midnight-ONT/V{1,2,3} ARTIC/V{1,2,3,4,4.1}", metavar="STR")
required_ops_group.add_argument("--reads", required=True, help="FASTQ file of ONT reads", metavar="FILENAME")

parser.add_argument("--main_nf", help="Path of main.nf file in the epi2me-labs/wf-artic repository [%(default)s]", metavar="FILENAME", default="/nfs/research/zi/mhunt/Covid/Epi2me/wf-artic/main.nf")
parser.add_argument("--nxf_sing_cache", help="Value of NXF_SINGULARITY_CACHEDIR [%(default)s]", metavar="DIRNAME", default="/nfs/research/zi/mhunt/Containers/nf_cache")
parser.add_argument("--debug", action="store_true", help="Keep temporary files")
parser.add_argument("--force", action="store_true", help="Overwrite outdir if it already exists")
parser.add_argument("--sample_name", help="Name of sample", metavar="STR")
parser.add_argument("--work_root_dir", help="Directory in which to put nextflow work directory", metavar="FILENAME")

options = parser.parse_args()

if options.force:
    subprocess.check_output(f"rm -rf {options.outdir}", shell=True)



run_epi2me(
        options.outdir,
        options.main_nf,
        options.scheme_version,
        options.reads,
        options.nxf_sing_cache,
        sample_name=options.sample_name,
        work_dir_root=options.work_root_dir,
        debug=options.debug,
)
