import os
import subprocess

def combine_vcfs(files,
                 output_dir,
                 file_no,
                 method,
                 vcf_dir=None):
    """ combine vcfs into a single vcf """
    vcf_rows = {}
    all_samples = []
    combined_vcf = []
    for f in files:
        if not method == "simulated":
            vcf_file = os.path.join(f, "variants.vcf")
        else:
            vcf_file = os.path.join(vcf_dir, os.path.basename(f), "04.truth.vcf")
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

def initiate_phylogeny(method,
                       assemblies,
                       assembly_dir,
                       vcf_dir,
                       source,
                       output_dir):
    """ initiate the phylogeny with a single sample """
    with open(os.path.join(output_dir, method + "_assemblies", method + "_tree.nwk"), "w") as outTree:
            outTree.write("(" + os.path.basename(assemblies[0]) + ")")
    tree_dir = os.path.join(output_dir, method + "_assemblies")
    tree_out = os.path.join(tree_dir, method + "_phylo.pb")
    if source == "viridian":
        vcf_file = os.path.join(vcf_dir, os.path.basename(assemblies[0]), "variants.vcf")
    if source == "simulated":
        vcf_file = os.path.join(vcf_dir, os.path.basename(assemblies[0]), "04.truth.vcf")
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
    usher_command = "singularity exec singularity/usher/usher.sif usher --vcf " + combined_vcf + " --load-mutation-annotated-tree " + \
            tree_out + " --write-uncondensed-final-tree -T " + str(threads) + " --save-mutation-annotated-tree " + tree_out +" --outdir " + tree_dir
    subprocess.run(usher_command, shell=True, check=True)

def optimise_phylogeny(tree_file,
                       threads):
    """ optimise the phylogeny using matoptimise """
    optimise_command = "singularity exec singularity/usher/usher.sif matOptimize -i " + tree_file + " -o " + tree_file + \
        "_optimised -T " + str(threads) + " -r 4 && rm -rf " + tree_file + " && mv " + tree_file + "_optimised " + tree_file 
    subprocess.run(optimise_command, shell=True, check=True)
