from datetime import datetime
import os
import subprocess

# make images directory
if not os.path.exists("singularity/images"):
    os.mkdir("singularity/images")
#download nextflow binary
download_command = "wget https://github.com/nextflow-io/nextflow/releases/download/v21.04.3/nextflow-21.04.3-all"
download_command += " && mv nextflow-21.04.3-all nextflow && chmod 777 nextflow"
subprocess.run(download_command, shell=True, check=True)
end
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
    subprocess.run(build_command, shell=True, check=True)
# build viridian image
viridian_command = "cd singularity && git clone https://github.com/iqbal-lab-org/viridian_workflow"
viridian_command += " && cd viridian_workflow && "
viridian_command += "sudo singularity build --force viridian_workflow.img Singularity.def && cd .."
subprocess.run(viridian_command, shell=True, check=True)
# build varifier image
varifier_command = "cd singularity && git clone https://github.com/iqbal-lab-org/varifier"
varifier_command += " && cd varifier && "
varifier_command += "sudo singularity build --force varifier.img Singularity.def && cd .."
subprocess.run(varifier_command, shell=True, check=True)
# build cte image
cte_command = "cd singularity && git clone https://github.com/iqbal-lab-org/covid-truth-eval"
cte_command += " && cd covid-truth-eval && "
cte_command += "sudo singularity build --force cte.img Singularity.def && cd .."
subprocess.run(cte_command, shell=True, check=True)
# build usher image
usher_command = "cd singularity && git clone https://github.com/yatisht/usher.git"
usher_command += " && cd usher && "
usher_command += "docker build --no-cache -t usher_docker install/."
usher_command += " && sudo singularity build --force usher.sif docker-daemon://usher_docker:latest && cd .."
subprocess.run(usher_command, shell=True, check=True)
# fetch ncov2019-artic-nf
artic_nf_command = "cd singularity && git clone https://github.com/Danderson123/ncov2019-artic-nf && "
artic_nf_command += "cd ncov2019-artic-nf && sudo singularity build --force artic_illumina.sif environments/illumina/Singularity && "
artic_nf_command += "&& sudo singularity build --force artic_nanopore.sif environments/nanopore/Singularity && cd ../.."
subprocess.run(artic_nf_command, shell=True, check=True)
# write out file stating when containers were built for version control
now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
with open("singularity/build_date.txt", "w") as outDate:
    outDate.write(dt_string)
