import os
import subprocess

def run_cte(primer_scheme,
            assemblies,
            outdir,
            truth_vcf_dir,
            container_dir,
            method):
    """run covid-truth-eval"""
    # define the primer scheme
    if primer_scheme == "V3":
        scheme = "COVID-ARTIC-V3"
    if primer_scheme == "V4":
        scheme = "COVID-ARTIC-V4"
    if primer_scheme == "V4.1":
        scheme = "COVID-MIDNIGHT-1200"
    # iterate through assemblies
    manifest = ["name\ttruth_vcf\teval_fasta\tprimers"]
    count = 0
    for assem in assemblies:
        vcf_file = os.path.join(truth_vcf_dir, os.path.basename(assem), "04.truth_dropped.vcf")
        #vcf_file = os.path.join(truth_vcf_dir, os.path.basename(assem), "04.truth.vcf")
        if method == "artic":
            assembly_file = os.path.join(assem, "consensus_trimmed.fa")
        else:
            assembly_file = os.path.join(assem, "consensus.fa.gz")
            #assembly_file = os.path.join(assem, "masked.fasta")
        #for a in assemblies:
        #   manifest.append(os.path.basename(a) + "_" + str(count) + "\t" + vcf_file + "\t" + os.path.join(a, "consensus.fa") + "\t" + scheme)
        #  count += 1
        manifest.append(os.path.basename(assem) + "\t" + vcf_file + "\t" + assembly_file + "\t" + scheme)
    # save metadata as tsv file
    with open(os.path.join(outdir, "manifest.tsv"), "w") as manifestOut:
        manifestOut.write("\n".join(manifest))
    # run covid-truth-eval
    shell_command = "singularity run " + container_dir + "/covid-truth-eval/cte.img eval_runs \
        --outdir " + os.path.join(outdir, "OUT") + " " + os.path.join(outdir, "manifest.tsv")
    subprocess.run(shell_command, shell=True, check=True)

def run_varifier(assembly,
                truth_vcf_dir,
                reference_genome,
                simulated_genomes,
                output,
                container_dir,
                artic_vcf_dir=None):
    """Run varifier on viridian assembly to verify SNP calls"""
    if not artic_vcf_dir:
        vcf_file = os.path.join(assembly, "variants.vcf")
        output_dir = os.path.join(output, "viridian_assemblies", os.path.join(assembly.split("/")[-2], os.path.basename(assembly)))
    else:
        vcf_file = os.path.join(artic_vcf_dir, "/".join(assembly.split("/")[-2:]), "04.truth.vcf")
        output_dir = os.path.join(output, "artic_assemblies", os.path.join(assembly.split("/")[-2], os.path.basename(assembly)))
    truth_vcf = os.path.join(truth_vcf_dir, os.path.basename(assembly), "04.truth.vcf")
    truth_genome = os.path.join(simulated_genomes, os.path.basename(assembly) + ".fasta")
    varifier_command = "singularity run " + container_dir + "/varifier/varifier.img vcf_eval --truth_vcf " + truth_vcf + " "
    # varifier_command = 'singularity exec "docker://quay.io/iqballab/varifier" varifier vcf_eval --truth_vcf ' + truth_vcf + ' '
    varifier_command += truth_genome + " " + reference_genome + " " + vcf_file + " " + output_dir
    subprocess.run(varifier_command, shell=True, check=True)
