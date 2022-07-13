
mafftScript=scripts/global_profile_alignment.sh
taxoniumGb=hu1.gb
refFa=reference_genome.fasta
numThreads=64

# make output dir
mkdir usher_phylogenies

echo making truth phylogeny
# Make monolithic MSA fasta -> VCF from truth assemblies:
cp /dev/null usher_phylogenies/truth_assemblies.fa
for id in masked_truth_assemblies/*; do
    echo "${id##*/}" >> usher_phylogenies/truth_assemblies.fa
    tail -n+2 $id >> usher_phylogenies/truth_assemblies.fa
    echo ""  >> usher_phylogenies/truth_assemblies.fa
done
time bash $mafftScript -i usher_phylogenies/truth_assemblies.fa -t $numThreads -o usher_phylogenies/truth_assemblies.msa.fa
scripts/faToVcf -includeNoAltN <(cat $refFa usher_phylogenies/truth_assemblies.msa.fa) \
    usher_phylogenies/truth_assemblies.vcf
# make tree from the vcf of all truth assemblies
time singularity exec singularity/usher/usher.sif usher --sort-before-placement-3 \
    --vcf usher_phylogenies/truth_assemblies.vcf \
    --tree <(echo '()') \
    -T $numThreads \
    --save-mutation-annotated-tree usher_phylogenies/truth_assemblies.pb \
    >& usher_phylogenies/truth_assemblies.usher.log
# Taxonium on unoptimised tree
usher_to_taxonium --input usher_phylogenies/truth_assemblies.pb \
    --genbank $taxoniumGb \
    --title "Daniel's faToVcf UShER (not optimized): truth assemblies" \
    --output usher_phylogenies/truth_assemblies.taxonium.jsonl.gz
# Optimise unoptimised truth tree
time singularity exec singularity/usher/usher.sif matOptimize -T $numThreads -r 8 -M 2 \
    -i usher_phylogenies/truth_assemblies.pb -o usher_phylogenies/truth_assemblies.opt.pb \
    >& usher_phylogenies/truth_assemblies.matOptimize.log
# Taxonium on optimised trees
usher_to_taxonium --input usher_phylogenies/truth_assemblies.opt.pb \
    --genbank $taxoniumGb \
    --title "Daniel's faToVcf UShER (optimized): truth assemblies" \
    --output usher_phylogenies/truth_assemblies.opt.taxonium.jsonl.gz

