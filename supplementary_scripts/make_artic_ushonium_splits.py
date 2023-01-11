import glob
import os
from tqdm import tqdm

viridian_ushonium_dirs = glob.glob("viridian_ushonium_trees/*")

if not os.path.exists("artic_ushonium_trees"):
	os.mkdir("artic_ushonium_trees")
for d in tqdm(viridian_ushonium_dirs):
	if not os.path.exists(os.path.join("artic_ushonium_trees", os.path.basename(d))):
		os.mkdir(os.path.join("artic_ushonium_trees", os.path.basename(d)))
	viridianTSV = os.path.join(d, "samples.tsv")
	with open(viridianTSV, "r") as i:
		tsvContent = i.read()
	newTsvContent = tsvContent.replace("viridian_ART_assemblies", "artic_ART_assemblies").replace("viridian_Badread_assemblies", "epi2me_Badread_assemblies")
	assert not "viridian_ART_assemblies" in newTsvContent
	assert not "viridian_Badread_assemblies" in newTsvContent
	with open(os.path.join("artic_ushonium_trees", os.path.basename(d), "samples.tsv"), "w") as o:
		o.write(newTsvContent)

