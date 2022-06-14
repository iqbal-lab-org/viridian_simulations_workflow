
mafftScript=scripts/global_profile_alignment.sh
taxoniumGb=hu1.gb
refFa=reference_genome.fasta
numThreads=64

# make output dir
mkdir usher_phylogenies

# Make monolithic MSA fasta -> VCF from viridian illumina:
cp /dev/null usher_phylogenies/viridian_illumina.fa
for id in $(virdidian_assemblies/ART_assemblies/*); do
    echo ">$id" >> usher_phylogenies/viridian_illumina.fa
    tail -n+2 virdidian_assemblies/ART_assemblies/$id.fasta >> usher_phylogenies/viridian_illumina.fa
    echo ""  >> usher_phylogenies/viridian_illumina.fa
done
time bash $mafftScript -i usher_phylogenies/viridian_illumina.fa -t $numThreads -o usher_phylogenies/viridian_illumina.msa.fa
./faToVcf -includeNoAltN <(cat $refFa usher_phylogenies/viridian_illumina.msa.fa) \
    usher_phylogenies/viridian_illumina.vcf

# Make monolithic MSA fasta -> VCF from viridian nanopore:
cp /dev/null usher_phylogenies/viridian_nanopore.fa
for id in $(virdidian_assemblies/Badread_assemblies/*); do
    echo ">$id" >> usher_phylogenies/viridian_nanopore.fa
    tail -n+2 virdidian_assemblies/Badread_assemblies/$id.fasta >> usher_phylogenies/viridian_nanopore.fa
    echo ""  >> usher_phylogenies/viridian_nanopore.fa
done
time bash $mafftScript -i usher_phylogenies/viridian_nanopore.fa -t $numThreads -o usher_phylogenies/viridian_nanopore.msa.fa
./faToVcf -includeNoAltN <(cat $refFa usher_phylogenies/viridian_nanopore.msa.fa) \
    viridian_nanopore.vcf

# Make monolithic MSA fasta -> VCF from artic illumina:
cp /dev/null usher_phylogenies/artic_illimina.fa
for id in $(artic_assemblies/ART_assemblies/*); do
    echo ">$id" >> usher_phylogenies/artic_illimina.fa
    tail -n+2 artic_assemblies/ART_assemblies/$id.fasta >> usher_phylogenies/artic_illimina.fa
    echo ""  >> usher_phylogenies/artic_illimina.fa
done
time bash $mafftScript -i usher_phylogenies/artic_illimina.fa -t $numThreads -o usher_phylogenies/artic_illimina.msa.fa
./faToVcf -includeNoAltN <(cat $refFa usher_phylogenies/artic_illimina.msa.fa) \
    usher_phylogenies/artic_illimina.vcf

# Make monolithic MSA fasta -> VCF from artic nanopore:
cp /dev/null usher_phylogenies/artic_nanopore.fa
for id in $(artic_assemblies/Badread_assemblies/*); do
    echo ">$id" >> usher_phylogenies/artic_nanopore.fa
    tail -n+2 artic_assemblies/Badread_assemblies/$id.fasta >> usher_phylogenies/artic_nanopore.fa
    echo ""  >> usher_phylogenies/artic_nanopore.fa
done
time bash $mafftScript -i usher_phylogenies/artic_nanopore.fa -t $numThreads -o usher_phylogenies/artic_nanopore.msa.fa
./faToVcf -includeNoAltN <(cat $refFa usher_phylogenies/artic_nanopore.msa.fa) \
    usher_phylogenies/artic_nanopore.vcf

# Make trees the impatient way: usher on whole VCF starting from empty tree.
for assembler in artic viridian; do
    for tech in ART_assemblies Badread_assemblies; do
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
    for tech in ART_assemblies Badread_assemblies; do
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
    for tech in ART_assemblies Badread_assemblies; do
        time singularity exec ../singularity/usher/usher.sif matOptimize -T $numThreads -r 8 -M 2 \
            -i usher_phylogenies/${assembler}_${tech}.pb -o usher_phylogenies/${assembler}_${tech}.opt.pb \
            >& usher_phylogenies/${assembler}_${tech}.matOptimize.log
    done
done

# Taxonium on optimized trees
for assembler in artic viridian; do
    for tech in ART_assemblies Badread_assemblies; do
        singularity exec ~../singularity/usher/usher.sif matUtils mask -i usher_phylogenies/${assembler}_${tech}.opt.pb \
            -o usher_phylogenies/${assembler}_${tech}.opt.renamed.pb >& renaming.out
        usher_to_taxonium --input usher_phylogenies/${assembler}_${tech}.opt.renamed.pb \
        --genbank $taxoniumGb \
        --title "Daniel's faToVcf UShER (optimized): $assembler $tech" \
        --output usher_phylogenies/${assembler}_${tech}.opt.taxonium.jsonl.gz
    done
done

