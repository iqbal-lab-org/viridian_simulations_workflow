
mafftScript=scripts/global_profile_alignment.sh
taxoniumGb=scripts/hu1.gb
refFa=reference_genome.fasta
numThreads=8

# make output dir
mkdir usher_phylogenies usher_phylogenies/truth_phylogenies usher_phylogenies/viridian_phylogenies usher_phylogenies/artic_phylogenies

echo making truth phylogeny
# Make monolithic MSA fasta -> VCF from truth assemblies:
cp /dev/null usher_phylogenies/truth_phylogenies/truth_assemblies.fa
for id in masked_truth_assemblies/*; do
    echo ">${id##*/}" >> usher_phylogenies/truth_phylogenies/truth_assemblies.fa
    tail -n+2 $id >> usher_phylogenies/truth_phylogenies/truth_assemblies.fa
    echo ""  >> usher_phylogenies/truth_phylogenies/truth_assemblies.fa
done
time bash $mafftScript -i usher_phylogenies/truth_phylogenies/truth_assemblies.fa -t $numThreads -o usher_phylogenies/truth_phylogenies/truth_assemblies.msa.fa
scripts/faToVcf -includeNoAltN <(cat $refFa usher_phylogenies/truth_phylogenies/truth_assemblies.msa.fa) \
    usher_phylogenies/truth_phylogenies/truth_assemblies.vcf
# make tree from the vcf of all truth assemblies
time singularity exec singularity/usher/usher.sif usher --sort-before-placement-3 \
    --vcf usher_phylogenies/truth_phylogenies/truth_assemblies.vcf \
    --tree <(echo '()') \
    -T $numThreads \
    --save-mutation-annotated-tree usher_phylogenies/truth_phylogenies/truth_assemblies.pb \
    >& usher_phylogenies/truth_phylogenies/truth_assemblies.usher.log
# Taxonium on unoptimised tree
usher_to_taxonium --input usher_phylogenies/truth_phylogenies/truth_assemblies.pb \
    --genbank $taxoniumGb \
    --title "Daniel's faToVcf UShER (not optimized): truth assemblies" \
    --output usher_phylogenies/truth_phylogenies/truth_assemblies.taxonium.jsonl.gz
# Optimise unoptimised truth tree
time singularity exec singularity/usher/usher.sif matOptimize -T $numThreads -r 8 -M 2 \
    -i usher_phylogenies/truth_phylogenies/truth_assemblies.pb -o usher_phylogenies/truth_phylogenies/truth_assemblies.opt.pb \
    >& usher_phylogenies/truth_phylogenies/truth_assemblies.matOptimize.log
# Taxonium on optimised trees
usher_to_taxonium --input usher_phylogenies/truth_phylogenies/truth_assemblies.opt.pb \
    --genbank $taxoniumGb \
    --title "Daniel's faToVcf UShER (optimized): truth assemblies" \
    --output usher_phylogenies/truth_phylogenies/truth_assemblies.opt.taxonium.jsonl.gz

# Make monolithic MSA fasta -> VCF from viridian illumina and nanopore assemblies:
for tech in ART Badread; do
    cp /dev/null usher_phylogenies/viridian_phylogenies/viridian_${tech}.fa
    for id in viridian_${tech}_assemblies/*; do
        echo ">${id##*/}" >> usher_phylogenies/viridian_phylogenies/viridian_${tech}.fa
        tail -n+2 $id/masked.fasta >> usher_phylogenies/viridian_phylogenies/viridian_${tech}.fa
        echo ""  >> usher_phylogenies/viridian_phylogenies/viridian_${tech}.fa
    done
    time bash $mafftScript -i usher_phylogenies/viridian_phylogenies/viridian_${tech}.fa -t $numThreads -o usher_phylogenies/viridian_phylogenies/viridian_${tech}.msa.fa
    scripts/faToVcf -includeNoAltN <(cat $refFa usher_phylogenies/viridian_phylogenies/viridian_${tech}.msa.fa) \
        usher_phylogenies/viridian_phylogenies/viridian_${tech}.vcf
done

# Make monolithic MSA fasta -> VCF from artic illumina and nanopore assemblies:
for tech in ART Badread; do
    cp /dev/null usher_phylogenies/artic_phylogenies/artic_${tech}.fa
    for id in artic_${tech}_assemblies/*; do
        echo ">${id##*/}" >> usher_phylogenies/artic_phylogenies/artic_${tech}.fa
        tail -n+2 $id/consensus.fa >> usher_phylogenies/artic_phylogenies/artic_${tech}.fa
        echo ""  >> usher_phylogenies/artic_phylogenies/artic_${tech}.fa
    done
    time bash $mafftScript -i usher_phylogenies/artic_phylogenies/artic_${tech}.fa -t $numThreads -o usher_phylogenies/artic_phylogenies/artic_${tech}.msa.fa
    scripts/faToVcf -includeNoAltN <(cat $refFa usher_phylogenies/artic_phylogenies/artic_${tech}.msa.fa) \
        usher_phylogenies/artic_phylogenies/artic_${tech}.vcf
done

# Make trees the impatient way: usher on whole VCF starting from empty tree.
for assembler in artic viridian; do
    for tech in ART Badread; do
        echo ${assembler}_${tech}
        time singularity exec singularity/usher/usher.sif usher --sort-before-placement-3 \
            --vcf usher_phylogenies/${assembler}_phylogenies/${assembler}_${tech}.vcf \
            --tree <(echo '()') \
            -T $numThreads \
            --save-mutation-annotated-tree usher_phylogenies/${assembler}_phylogenies/${assembler}_${tech}.pb \
            >& usher_phylogenies/${assembler}_phylogenies/${assembler}_${tech}.usher.log
    done
done

# Taxonium on unoptimized trees
for assembler in artic viridian; do
    for tech in ART Badread; do
        usher_to_taxonium --input usher_phylogenies/${assembler}_phylogenies/${assembler}_${tech}.pb \
        --genbank $taxoniumGb \
        --title "Daniel's faToVcf UShER (not optimized): $assembler $tech" \
        --output usher_phylogenies/${assembler}_phylogenies/${assembler}_${tech}.taxonium.jsonl.gz
    done
done

# Optimize
for assembler in artic viridian; do
    for tech in ART Badread; do
        time singularity exec singularity/usher/usher.sif matOptimize -T $numThreads -r 8 -M 2 \
            -i usher_phylogenies/${assembler}_phylogenies/${assembler}_${tech}.pb -o usher_phylogenies/${assembler}_phylogenies/${assembler}_${tech}.opt.pb \
            >& usher_phylogenies/${assembler}_phylogenies/${assembler}_${tech}.matOptimize.log
    done
done

# Taxonium on optimized trees
for assembler in artic viridian; do
    for tech in ART Badread; do
        usher_to_taxonium --input usher_phylogenies/${assembler}_phylogenies/${assembler}_${tech}.opt.pb \
        --genbank $taxoniumGb \
        --title "Daniel's faToVcf UShER (optimized): $assembler $tech" \
        --output usher_phylogenies/${assembler}_phylogenies/${assembler}_${tech}.opt.taxonium.jsonl.gz
    done
done
