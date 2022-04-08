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



def run_ncov_2019_artic_nf(
  outdir,
  sif,
  nf_script,
  scheme_url,
  nextflow_path,
  schemeVersion,
  ilm1=None,
  ilm2=None,
  ont=None,
  debug=False,
  sample_name=None,
  work_dir_root=None,
  ):
    outdir = os.path.abspath(outdir)
    log_info = {
        "inputs": {
            "outdir": outdir,
            "sif": sif,
            "nf_script": nf_script,
            "ilm1": ilm1,
            "ilm2": ilm2,
            "ont": ont,
            "debug": debug,
            "sample_name": sample_name,
        }
    }
    sif = os.path.abspath(sif)
    nf_script = os.path.abspath(nf_script)

    os.mkdir(outdir)
    reads_dir = "Reads"
    reads_dir_abs = os.path.join(outdir, reads_dir)
    os.mkdir(reads_dir_abs)
    if ont is None:
        os.symlink(os.path.abspath(ilm1), os.path.join(reads_dir_abs, "reads_1.fq.gz"))
        os.symlink(os.path.abspath(ilm2), os.path.join(reads_dir_abs, "reads_2.fq.gz"))
        tech_opt = "--illumina"
        cpus = 2
        reads_opt = "--directory"
    else:
        # Make a fake guppy barcoded directory
        barcode_dir = os.path.join(reads_dir_abs, "barcode01")
        os.mkdir(barcode_dir)
        os.symlink(os.path.abspath(ont), os.path.join(os.path.abspath(barcode_dir), "sample_pass_barcode01_foo_0.fastq.gz"))
        tech_opt = "--medaka"
        cpus = 1
        reads_opt = "--basecalled_fastq"

    os.chdir(outdir)
    if work_dir_root is not None:
        assert os.path.exists(work_dir_root)
        work_dir = os.path.join(os.path.abspath(work_dir_root), f"{sample_name}.cog")
        subprocess.check_output(f"rm -rf {work_dir}", shell=True)
    else:
        work_dir = "tmp_nf_work"
        config_file = "nextflow.config"

    with open(config_file, "w") as f:
        print("executor.queueSize = 1",file=f)
        # For illumina, can't limit to 1 cpu, because the read trimming uses 2 cpus.
        # Looks like that's what trim galore uses (even though the nextflow
        # pipeline is using default of 1 thread). So even if we force
        # 1 cpu here, the actual job would use 2 anyway.
        print("executor.cpus =", cpus,file=f)
        print("process.cpus =", cpus, file=f)

    command = " ".join([
        "/usr/bin/time -v",
        nextflow_path + " run", os.path.abspath(nf_script),
        "-with-singularity", sif,
        "-work-dir", work_dir,
        "-ansi-log false",
        "-c", config_file,
        tech_opt,
        "--schemeRepoURL", scheme_url,
        "--schemeDir", "primer-schemes",
        "--schemeVersion", schemeVersion,
        "--scheme", "SARS-CoV-2",
        "--prefix out",
        reads_opt, reads_dir,
    ])
    completed_process = syscall(command)
    log_info["nextflow"] = {
        "stdout": completed_process.stdout.rstrip().split("\n"),
        "stderr": completed_process.stderr.rstrip().split("\n"),
        "command": command,
    }
    log_info["nextflow"].update(bash_out_to_time_and_memory(log_info["nextflow"]["stderr"]))
    with open("log.json", "w") as f:
        json.dump(log_info, f, indent=2)

    if tech_opt == "--illumina":
        os.rename("results/ncovIllumina_sequenceAnalysis_callVariants/reads.variants.tsv", "variants.tsv")
        old_fa = "results/ncovIllumina_sequenceAnalysis_makeConsensus/reads.primertrimmed.consensus.fa"
    elif tech_opt == "--medaka":
        medaka_dir = os.path.join("results", "articNcovNanopore_sequenceAnalysisMedaka_articMinIONMedaka")
        old_fa = os.path.join(medaka_dir, "out_barcode01.consensus.fasta")
        keep_dir = "vcfs_etc"
        os.mkdir(keep_dir)
        subprocess.check_output(f"mv {medaka_dir}/*vcf* {medaka_dir}/*preconsensus.fasta {keep_dir}/", shell=True)
        subprocess.check_output(f"tar cf {keep_dir}.tar {keep_dir}/", shell=True)
        subprocess.check_output(f"gzip -9 {keep_dir}.tar", shell=True)
        subprocess.check_output(f"rm -r {keep_dir}", shell=True)

    seqs = {}
    pyfastaq.tasks.file_to_dict(old_fa, seqs)
    assert len(seqs) == 1
    seq = list(seqs.values())[0]
    if sample_name is not None:
        seq.id = sample_name + " " + seq.id
    with open("consensus.fa", "w") as f:
        print(seq, file=f)

    if not debug:
        subprocess.check_output(f"rm -rf {work_dir} {config_file} .nextflow .nextflow.log {reads_dir} results", shell=True)

parser = argparse.ArgumentParser(
    description = "Run ncov2019-artic-nf on one sample",
    usage="%(prog)s [options]",
)
required_ops_group = parser.add_argument_group("Required options")
required_ops_group.add_argument("--sif", required=True, help="Name of artic-ncov2019.sif file", metavar="FILENAME")
required_ops_group.add_argument("--main_nf", required=True, help="Path of main.nf file in the ncov2019-artic-nf repository", metavar="FILENAME")
required_ops_group.add_argument("--outdir", required=True, help="Output directory")
required_ops_group.add_argument("--scheme_url", required=True, help="URL of repo defining the primer scheme")
required_ops_group.add_argument("--scheme_version", required=True, help="primer scheme version to use")
required_ops_group.add_argument("--nextflow_path", required=True, help="abs path for nextflow binary")

parser.add_argument("--ilm1", help="FASTQ file 1 of Illumina reads", metavar="FILENAME")
parser.add_argument("--ilm2", help="FASTQ file 2 of Illumina reads", metavar="FILENAME")
parser.add_argument("--ont", help="FASTQ file of ONT reads")
parser.add_argument("--force", action="store_true", help="Overwrite outdir if it already exists")
parser.add_argument("--sample_name", help="Name of sample", metavar="STRING")
parser.add_argument("--work_root_dir", help="Directory in which to put nextflow work directory", metavar="FILENAME")

options = parser.parse_args()

if options.ont is None:
    assert options.ilm1 is not None and options.ilm2 is not None
else:
    assert options.ilm1 is None and options.ilm2 is None


if options.force:
    subprocess.check_output(f"rm -rf {options.outdir}", shell=True)



run_ncov_2019_artic_nf(
        options.outdir,
        options.sif,
        options.main_nf,
        options.scheme_url,
        options.nextflow_path,
        options.scheme_version,
        ilm1=options.ilm1,
        ilm2=options.ilm2,
        ont=options.ont,
        sample_name=options.sample_name,
        work_dir_root=options.work_root_dir,
)



