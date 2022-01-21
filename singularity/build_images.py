import subprocess

recipes = ["VGsim", "phastSim", "ART", "Badread"]
for recipe in recipes:
    # build image
    build_command = "singularity build --fakeroot singularity/images/"
    build_command += recipe + ".img singularity/recipes/"
    build_command += recipe + ".def"
    subprocess.run(build_command, shell=True, check=True)
# build viridian image
viridian_command = "git clone https://github.com/iqbal-lab-org/viridian"
viridian_command += " && cd viridian && "
viridian_command += "singularity build --fakeroot viridian.img Singularity.def && "
viridian_command += "cd .."
subprocess.run(viridian_command, shell=True, check=True)
# build varifier image
varifier_command = "git clone https://github.com/iqbal-lab-org/varifier"
varifier_command += " && cd varifier && "
varifier_command += "singularity build --fakeroot varifier.img Singularity.def && "
varifier_command += "cd .."
subprocess.run(varifier_command, shell=True, check=True)