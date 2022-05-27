import os
import subprocess

import shutil 
import tempfile

def combine_vcfs_old(files,
                 output_dir,
                 file_no,
                 method,
                 vcf_dir=None):
    """ combine vcfs into a single vcf """
    vcf_rows = {}
    all_samples = []
    combined_vcf = []
    for f in files:
        if method == "viridian":
            vcf_file = os.path.join(f, "variants.vcf")
        if method == "simulated":
            vcf_file = os.path.join(vcf_dir, os.path.basename(f), "04.truth.vcf")
        if method == "artic":
            vcf_file = os.path.join(vcf_dir, os.path.basename(output_dir), os.path.basename(f), "04.truth.vcf")
        with open(vcf_file, "r") as inVCF:
            vcf_content = inVCF.read().split("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample")[1].splitlines()[1:]
        sample = os.path.basename(f)
        all_samples.append(sample)
        for v in vcf_content:
            if not v.replace("\t1/1", "") in vcf_rows:
                vcf_rows[v.replace("\t1/1", "")] = {sample: "1"}
            else:
                vcf_rows[v.replace("\t1/1", "")].update({sample: "1"})
    for row in vcf_rows:
        new_row = row
        for s in all_samples:
            if not s in vcf_rows[row]:
                new_row += "\t" + "0"
            if s in vcf_rows[row]:
                new_row += "\t" + vcf_rows[row][s]
        combined_vcf.append(new_row)
    sorted(combined_vcf, key=lambda x: x.split("	")[1])
    out_vcf = os.path.join(output_dir, str(file_no) + ".vcf")
    with open(out_vcf, "w") as outFile:
        outFile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(all_samples) + "\n" + "\n".join(combined_vcf))
    return out_vcf

def combine_vcfs(file_paths,
                 output_dir,
                 method,
                 assembly_dir,
                 vcf_dir=None):
    temp_dir = tempfile.mkdtemp(dir=output_dir)
    gzipped_files = []
    for f in file_paths:
        sample = os.path.basename(f)
        if method == "viridian":
            old_vcf_file = os.path.join(f, "variants.vcf")
            old_assembly = os.path.join(f, "consensus.fa")
        if method == "simulated":
            old_vcf_file = os.path.join(vcf_dir, os.path.basename(f), "04.truth.vcf")
            old_assembly = os.path.join(assembly_dir, os.path.basename(f) + ".fasta")
        if method == "artic":
            old_vcf_file = os.path.join(vcf_dir, os.path.basename(output_dir), os.path.basename(f), "04.truth.vcf")
            old_assembly = os.path.join(assembly_dir, os.path.basename(output_dir), os.path.basename(f), "consensus.fa")
        new_vcf_file = os.path.join(temp_dir, sample + ".vcf")
        new_assembly = os.path.join(temp_dir, sample + ".fa")
        # copy vcf file into temp dir
        subprocess.run("cp " + old_vcf_file + " " + new_vcf_file, shell=True, check=True)
        subprocess.run("cp " + old_assembly + " " + new_assembly, shell=True, check=True)
        # replace sample name in vcf
        with open(new_vcf_file, "r") as inVCF:
            content = inVCF.read().replace("sample", sample.replace(".fasta", ""))
        with open(new_vcf_file, "w") as outVCF:
            outVCF.write(content)
        # compress and index vcf file
        index_command = "singularity run singularity/images/bcftools.img convert -Oz -o " + new_vcf_file + ".gz" + " " + new_vcf_file + " " + " && " 
        index_command += "singularity run singularity/images/bcftools.img index " + new_vcf_file + ".gz"
        subprocess.run(index_command, shell=True, check=True)
        gzipped_files.append(new_vcf_file + '.gz')
    output_path = os.path.join(output_dir, os.path.basename(temp_dir) + ".vcf")
    if not len(gzipped_files) == 1:
        merge_command = "singularity run singularity/images/bcftools.img merge --no-version --output " + output_path + " " + " ".join(gzipped_files)
    else:
        merge_command = "cp " + new_vcf_file + " " + output_path
    subprocess.run(merge_command, shell=True, check=True)
    shutil.rmtree(temp_dir)
    return output_path

def initiate_phylogeny(method,
                       assemblies,
                       assembly_dir,
                       vcf_dir,
                       source,
                       output_dir,
                       reference_sequence):
    """ initiate the phylogeny with a single sample """
    with open(os.path.join(output_dir, method + "_assemblies", method + "_tree.nwk"), "w") as outTree:
            outTree.write("(" + os.path.basename(reference_sequence.replace(".fasta", "")) + ")")
    tree_dir = os.path.join(output_dir, method + "_assemblies")
    tree_out = os.path.join(tree_dir, method + "_phylo.pb")
    vcf_file = reference_sequence.replace(".fasta", ".vcf")
    MAT_command = "singularity exec singularity/usher/usher.sif usher --tree " + os.path.join(output_dir,  method + "_assemblies", method + "_tree.nwk") + \
    " --vcf " + vcf_file + " --save-mutation-annotated-tree " + \
    tree_out + " -d " + tree_dir
    subprocess.run(MAT_command, shell=True, check=True)    
    return tree_dir, tree_out


def add_samples(combined_vcf,
                tree_out,
                tree_dir,
                threads):
    """ add samples to phylogeny with Usher """
    usher_command = "singularity exec singularity/usher/usher.sif usher --sort-before-placement-3 --vcf " + combined_vcf + " --load-mutation-annotated-tree " + \
            tree_out + " -T " + str(threads) + " --save-mutation-annotated-tree " + tree_out +" --outdir " + tree_dir
    subprocess.run(usher_command, shell=True, check=True)

def optimise_phylogeny(tree_file,
                       threads):
    """ optimise the phylogeny using matoptimise """
    optimise_command = "singularity exec singularity/usher/usher.sif matOptimize -n -i " + tree_file + " -o " + tree_file + \
        "_optimised -T " + str(threads) + " -r 4 && rm -rf " + tree_file + " && mv " + tree_file + "_optimised " + tree_file 
    subprocess.run(optimise_command, shell=True, check=True)
