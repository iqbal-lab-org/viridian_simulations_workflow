import os
import glob
import sys

from scripts.make_plots import run_cte

out_dir = "cte_viridian_output"
for o in ["cte_viridian_output", "cte_viridian_output/ART_assemblies", "cte_viridian_output/Badread_assemblies"]:
	if not os.path.exists(o):
		os.mkdir(o)
art_assemblies = glob.glob("viridian_ART_assemblies/*")
badread_assemblies = glob.glob("viridian_Badread_assemblies/*")
#run_cte("V4",
#	art_assemblies,
#	os.path.join(out_dir, "ART_assemblies"),
#	"truth_vcfs",
#	"singularity",
 #       "viridian")
run_cte("V4",
	badread_assemblies,
        os.path.join(out_dir, "Badread_assemblies"),
        "truth_vcfs",
        "singularity",
        "viridian")
sys.stderr.write("\nDone\n")
