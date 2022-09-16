import os
import subprocess

# make images directory
if not os.path.exists("singularity/images"):
    os.mkdir("singularity/images")
#download nextflow binary
download_command = "wget https://github.com/nextflow-io/nextflow/releases/download/v21.04.3/nextflow-21.04.3-all"
download_command += " && mv nextflow-21.04.3-all nextflow && chmod 777 nextflow"
#subprocess.run(download_command, shell=True, check=True)
# iterate through list of recipes to build the image
recipes = ["VGsim",
        "phastSim",
        "ART",
        "Badread",
        "Map",
        "Mafft",
        "bcftools"]
# build image
for recipe in recipes:
    build_command = "singularity build --force --fakeroot singularity/images/"
    build_command += recipe + ".img singularity/recipes/"
    build_command += recipe + ".def"
    #subprocess.run(build_command, shell=True, check=True)
# build viridian image
viridian_command = "cd singularity && git clone https://github.com/iqbal-lab-org/viridian_workflow"
viridian_command += " && cd viridian_workflow && git checkout cb6c7a55cb145d74fc6352341a76a183b0c1b4e1 && "
viridian_command += "sudo singularity build --force viridian_workflow.img Singularity.def && cd .."
#subprocess.run(viridian_command, shell=True, check=True)
# build varifier image
varifier_command = "cd singularity && git clone https://github.com/iqbal-lab-org/varifier"
varifier_command += " && cd varifier && git checkout f0c0bf0ae3cf686a47bf1be1a37b3efab54ffcb3 && "
varifier_command += "sudo singularity build --force varifier.img Singularity.def && cd .."
#subprocess.run(varifier_command, shell=True, check=True)
# build cte image
cte_command = "cd singularity && git clone https://github.com/iqbal-lab-org/covid-truth-eval"
cte_command += " && cd covid-truth-eval && git checkout 9cd94b8ce2ef54d5f6168369567dd49e2db92d2e && "
cte_command += "sudo singularity build --force cte.img Singularity.def && cd .."
#subprocess.run(cte_command, shell=True, check=True)
# build usher image
usher_command = "cd singularity && git clone https://github.com/yatisht/usher.git"
usher_command += " && cd usher && git checkout 4f620fe3d8cd7e013675e15a0fb9129a0fc502c5 && "
usher_command += "docker build --no-cache -t usher_docker install/. && "
usher_command += "sudo singularity build --force usher.sif docker-daemon://usher_docker:latest && cd .."
#subprocess.run(usher_command, shell=True, check=True)
# fetch ncov2019-artic-nf
artic_nf_command = "cd singularity && git clone https://github.com/connor-lab/ncov2019-artic-nf && "
artic_nf_command += "cd ncov2019-artic-nf && git checkout 8af5152cf7107c3a369b996f5bad3473d329050c && "
artic_nf_command += "sudo singularity build --force artic_illumina.sif environments/illumina/Singularity && "
artic_nf_command += "sudo singularity build --force artic_nanopore.sif environments/nanopore/Singularity && cd ../.."
#subprocess.run(artic_nf_command, shell=True, check=True)
# build epi2me image
epi2me_command = "cd singularity && git clone https://github.com/epi2me-labs/wf-artic epi2me && "
epi2me_command += "cd epi2me && git checkout a60a1e1e73bde1971c680bd3d53076127dd63fc6 && cd ../.."
subprocess.run(epi2me_command, shell=True, check=True)
# build ushonium image
ushonium_command = "cd singularity && git clone https://github.com/martinghunt/ushonium && "
ushonium_command += "cd ushonium && git checkout 930a0d16b576ec815b71284c8387edebee6cd9dc && "
ushonium_command += "sudo singularity build --force ushonium.img Singularity.def && cd ../.."
#subprocess.run(ushonium_command, shell=True, check=True)
