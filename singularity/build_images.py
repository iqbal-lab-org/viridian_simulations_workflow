import os
import subprocess

# make images directory
if not os.path.exists("singularity/images"):
    os.mkdir("signularity/images")
# iterate through list of recipes to build the image
recipes = ["VGsim", "phastSim", "ART", "Badread", "ErrorModes", "BWA", "Bedtools", "Samtools"]
for recipe in recipes:
    # build image
    build_command = "singularity build --force --fakeroot singularity/images/"
    build_command += recipe + ".img singularity/recipes/"
    build_command += recipe + ".def"
    subprocess.run(build_command, shell=True, check=True)
# build viridian image
viridian_command = "git clone https://github.com/iqbal-lab-org/viridian"
viridian_command += " && cd viridian && "
viridian_command += "singularity build --force --fakeroot viridian.img Singularity.def && "
viridian_command += "cd .."
subprocess.run(viridian_command, shell=True, check=True)
# build varifier image
varifier_command = "git clone https://github.com/iqbal-lab-org/varifier"
varifier_command += " && cd varifier && "
varifier_command += "singularity build --force --fakeroot varifier.img Singularity.def && "
varifier_command += "cd .."
subprocess.run(varifier_command, shell=True, check=True)