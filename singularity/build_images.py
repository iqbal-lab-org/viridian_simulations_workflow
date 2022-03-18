from datetime import datetime
import os
import subprocess

# make images directory
if not os.path.exists("singularity/images"):
    os.mkdir("singularity/images")
# iterate through list of recipes to build the image
recipes = ["VGsim",
           "phastSim",
           "ART",
           "Badread",
           "Map"]
# build image
for recipe in recipes:
    build_command = "singularity build --force --fakeroot singularity/images/"
    build_command += recipe + ".img singularity/recipes/"
    build_command += recipe + ".def"
    subprocess.run(build_command, shell=True, check=True)
# build viridian image
viridian_command = "cd singularity && git clone https://github.com/iqbal-lab-org/viridian_workflow"
viridian_command += " && cd viridian_workflow && "
viridian_command += "sudo build --force viridian_workflow.img Singularity.def && cd .."
subprocess.run(viridian_command, shell=True, check=True)
# build varifier image
varifier_command = "git clone https://github.com/iqbal-lab-org/varifier"
varifier_command += " && cd varifier && "
varifier_command += "sudo build --force varifier.img Singularity.def && cd .."
subprocess.run(varifier_command, shell=True, check=True)
# build cte image
varifier_command = "git clone https://github.com/iqbal-lab-org/covid-truth-eval"
varifier_command += " && cd covid-truth-eval && "
varifier_command += "sudo singularity build --force cte.img Singularity.def && cd .."
subprocess.run(varifier_command, shell=True, check=True)
# build usher image
usher_command = "git clone https://github.com/yatisht/usher.git"
usher_command += " && cd usher && "
usher_command += "docker build --no-cache -t usher_docker install/."
usher_command += " && sudo singularity build --force usher.sif docker-daemon://usher_docker:latest && cd .."
subprocess.run(usher_command, shell=True, check=True)
# fetch ncov2019-artic-nf
artic_nf_command = "git clone https://github.com/connor-lab/ncov2019-artic-nf && "
artic_nf_command += "cd ncov2019-artic-nf && sudo singularity build --force environments/illumina/artic_illumina.sif environments/illumina/Singularity && "
artic_nf_command += "&& sudo singularity build --force environments/nanopore/artic_nanopore.sif environments/nanopore/Singularity && cd ../.."
subprocess.run(artic_nf_command, shell=True, check=True)
# write out file stating when containers were built for version control
now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
with open("singularity/build_date.txt", "w") as outDate:
    outDate.write(dt_string)
