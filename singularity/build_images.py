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
viridian_command += "singularity build --force --fakeroot viridian_workflow.img Singularity.def && "
viridian_command += "cd .."
subprocess.run(viridian_command, shell=True, check=True)
# build varifier image
varifier_command = "git clone https://github.com/iqbal-lab-org/varifier"
varifier_command += " && cd varifier && "
varifier_command += "singularity build --force --fakeroot varifier.img Singularity.def && cd .."
subprocess.run(varifier_command, shell=True, check=True)
# build usher image
usher_command = "git clone https://github.com/yatisht/usher.git"
usher_command += " && cd usher && "
usher_command += "docker build --no-cache -t usher_docker install/."
usher_command += " && sudo singularity build --force usher.sif docker-daemon://usher_docker:latest"
subprocess.run(usher_command, shell=True, check=True)