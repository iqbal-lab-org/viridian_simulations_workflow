import os
import glob
import sys

from scripts.make_plots import run_cte

out_dir = "cte_artic_output"
for o in ["cte_artic_output", "cte_artic_output/ART_assemblies", "cte_artic_output/Badread_assemblies"]:
	if not os.path.exists(o):
		os.mkdir(o)
art_assemblies = glob.glob("artic_ART_assemblies/*")
badread_assemblies = glob.glob("artic_Badread_assemblies/*")
run_cte("V4",
	art_assemblies,
	os.path.join(out_dir, "ART_assemblies"),
	"truth_vcfs",
	"singularity",
        "artic")
#run_cte("V4",
#	badread_assemblies,
 #       os.path.join(out_dir, "Badread_assemblies"),
  #      "truth_vcfs",
   #     "singularity",
    #    "artic")
sys.stderr.write("\nDone\n")
