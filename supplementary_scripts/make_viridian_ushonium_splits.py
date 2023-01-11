import glob
import os
import random
from tqdm import tqdm
# list and sort the assemblies
viridian_ART_assemblies = list(sorted(glob.glob("viridian_ART_assemblies/*/consensus.fa"), key = lambda a: os.path.basename(os.path.dirname(a))))
random.shuffle(viridian_ART_assemblies)
viridian_Badread_assemblies = [v.replace("viridian_ART_assemblies", "viridian_Badread_assemblies") for v in viridian_ART_assemblies]
#viridian_Badread_assemblies = list(sorted(glob.glob("viridian_Badread_assemblies/*/consensus.fa"), key = lambda a: os.path.basename(os.path.dirname(a))))
# ensure the sample order is the same
for i in tqdm(range(len(viridian_ART_assemblies))):
	sample = os.path.basename(os.path.dirname(viridian_ART_assemblies[i]))
	otherSample = os.path.basename(os.path.dirname(viridian_Badread_assemblies[i]))
	assert sample == otherSample, "The sample ordering is not consistent"
# make the output directories
directory_splits = [["0N10I", "10N0I"], ["1N9I", "9N1I"], ["2N8I", "8N2I"], ["3N7I", "7N3I"], ["4N6I", "6N4I"], ["5N5I", "5N5I"]]
for d in directory_splits:
	for s in d:
		if not os.path.exists(os.path.join("ushonium_trees", s)):
			os.mkdir(os.path.join("ushonium_trees", s))
# define how we are splitting the assemblies
proportion_by_tech = [[0, 1], [0.1, 0.9], [0.2, 0.8], [0.3, 0.7], [0.4, 0.6], [0.5, 0.5]]
# iterate through the proportions and split the assemblies
for i in tqdm(range(len(proportion_by_tech))):
	# we are going to conserve the order of the split so we can directly compare how tech changes the tree
	firstProportion = proportion_by_tech[i][0]
	secondProportion = proportion_by_tech[i][1]
	firstCount = int(8000 * firstProportion)
	secondCount = int(8000 * secondProportion)
	# split the assemblies
	firstNanoporeSplit = viridian_Badread_assemblies[:firstCount]
	firstIlluminaSplit = viridian_ART_assemblies[-secondCount:]
	secondNanoporeSplit = viridian_Badread_assemblies[-secondCount:]
	secondIlluminaSplit = viridian_ART_assemblies[:firstCount]
	assert len(firstNanoporeSplit) + len(firstIlluminaSplit) == 8000
	assert len(secondNanoporeSplit) + len(secondIlluminaSplit) == 8000
	# make the tsv file for each split
	firstTsvLines = []
	for a in firstNanoporeSplit + firstIlluminaSplit:
		thisLine = os.path.basename(os.path.dirname(a)) + "\t" + a
		firstTsvLines.append(thisLine)
	secondTsvLines = []
	for a in secondIlluminaSplit + secondNanoporeSplit:
                thisLine = os.path.basename(os.path.dirname(a)) + "\t" + a
                secondTsvLines.append(thisLine)
	# ensure the splits are the same across each tech
	for t in range(len(firstTsvLines)):
		sample = firstTsvLines[t].split("\t")[0]
		otherSample = secondTsvLines[t].split("\t")[0]
		assert sample == otherSample, "The sample ordering is not consistent"
	# write out the tsv files to the correct directory
	with open(os.path.join("ushonium_trees", directory_splits[i][0], "samples.tsv"), "w") as o:
		o.write("\n".join(firstTsvLines))
	with open(os.path.join("ushonium_trees", directory_splits[i][1], "samples.tsv"), "w") as o:
                o.write("\n".join(secondTsvLines))

