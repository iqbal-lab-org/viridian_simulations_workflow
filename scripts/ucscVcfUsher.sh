
mafftScript=scripts/global_profile_alignment.sh
taxoniumGb=hu1.gb
refFa=reference_genome.fasta
numThreads=64

# make output dir
mkdir usher_phylogenies

# Make monolithic MSA fasta -> VCF from viridian illumina and nanopore:
for tech in ART Badread; do
    cp /dev/null usher_phylogenies/viridian_${tech}.fa
    for id in viridian_assemblies/${tech}_assemblies/*; do
        echo ">${id##*/}" >> usher_phylogenies/viridian_${tech}.fa
        tail -n+2 $id/masked.fasta >> usher_phylogenies/viridian_${tech}.fa
        echo ""  >> usher_phylogenies/viridian_${tech}.fa
    done
    time bash $mafftScript -i usher_phylogenies/viridian_${tech}.fa -t $numThreads -o usher_phylogenies/viridian_${tech}.msa.fa
    ./faToVcf -includeNoAltN <(cat $refFa usher_phylogenies/viridian_${tech}.msa.fa) \
        usher_phylogenies/viridian_${tech}.vcf
done

# Make monolithic MSA fasta -> VCF from artic illumina and nanopore:
for tech in ART Badread; do
    cp /dev/null usher_phylogenies/artic_${tech}.fa
    for id in artic_assemblies/${tech}_assemblies/*; do
        echo ">${id##*/}" >> usher_phylogenies/artic_${tech}.fa
        tail -n+2 $id/consensus.fa >> usher_phylogenies/artic_${tech}.fa
        echo ""  >> usher_phylogenies/artic_${tech}.fa
    done
    time bash $mafftScript -i usher_phylogenies/artic_${tech}.fa -t $numThreads -o usher_phylogenies/artic_${tech}.msa.fa
    ./faToVcf -includeNoAltN <(cat $refFa usher_phylogenies/artic_${tech}.msa.fa) \
        usher_phylogenies/artic_${tech}.vcf
done

# Make trees the impatient way: usher on whole VCF starting from empty tree.
for assembler in artic viridian; do
    for tech in ART Badread; do
        echo ${assembler}_${tech}
        time singularity exec ../singularity/usher/usher.sif usher --sort-before-placement-3 \
            --vcf usher_phylogenies/${assembler}_${tech}.vcf \
            --tree <(echo '()') \
            -T $numThreads \
            --save-mutation-annotated-tree usher_phylogenies/${assembler}_${tech}.pb \
            >& usher_phylogenies/${assembler}_${tech}.usher.log
    done
done

# Taxonium on unoptimized trees
for assembler in artic viridian; do
    for tech in ART Badread; do
        singularity exec ../singularity/usher/usher.sif matUtils mask -i usher_phylogenies/${assembler}_${tech}.pb \
            -o usher_phylogenies/${assembler}_${tech}.renamed.pb >& usher_phylogenies/renaming.out
        usher_to_taxonium --input usher_phylogenies/${assembler}_${tech}.renamed.pb \
        --genbank $taxoniumGb \
        --title "Daniel's faToVcf UShER (not optimized): $assembler $tech" \
        --output usher_phylogenies/${assembler}_${tech}.taxonium.jsonl.gz
    done
done

# Optimize
for assembler in artic viridian; do
    for tech in ART Badread; do
        time singularity exec ../singularity/usher/usher.sif matOptimize -T $numThreads -r 8 -M 2 \
            -i usher_phylogenies/${assembler}_${tech}.pb -o usher_phylogenies/${assembler}_${tech}.opt.pb \
            >& usher_phylogenies/${assembler}_${tech}.matOptimize.log
    done
done

# Taxonium on optimized trees
for assembler in artic viridian; do
    for tech in ART Badread; do
        singularity exec ~../singularity/usher/usher.sif matUtils mask -i usher_phylogenies/${assembler}_${tech}.opt.pb \
            -o usher_phylogenies/${assembler}_${tech}.opt.renamed.pb >& renaming.out
        usher_to_taxonium --input usher_phylogenies/${assembler}_${tech}.opt.renamed.pb \
        --genbank $taxoniumGb \
        --title "Daniel's faToVcf UShER (optimized): $assembler $tech" \
        --output usher_phylogenies/${assembler}_${tech}.opt.taxonium.jsonl.gz
    done
done

