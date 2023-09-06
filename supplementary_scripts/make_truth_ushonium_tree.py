import subprocess

ushonium_command = 'singularity run /nfs/research/zi/mhunt/Containers/ushonium.20230303.9efc65b.img --ref_start_end 25,29854 --title "Truth phylogeny" --cpus 64 ushonium_trees/truth_ushonium_trees/samples.tsv ushonium_trees/truth_ushonium_trees/ushonium_output'
subprocess.run(ushonium_command, shell=True, check=True)
