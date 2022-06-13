
mafftScript=scritps/global_profile_alignment.sh
taxoniumGb=hu1.gb
refFa=MN908947.3.fa
numThreads=64

# Make monolithic MSA fasta -> VCF from viridian illumina:
cp /dev/null viridian_illumina.fa
for id in $(virdidian_assemblies/ART_assemblies/*); do
    echo ">$id" >> viridian_illumina.fa
    tail -n+2 virdidian_assemblies/ART_assemblies/$id.fasta >> viridian_illumina.fa
    echo ""  >> viridian_illumina.fa
done
time bash $mafftScript -i viridian_illumina.fa -t $numThreads -o viridian_illumina.msa.fa
./faToVcf -includeNoAltN <(cat $refFa viridian_illumina.msa.fa) \
    viridian_illumina.vcf

# Make monolithic MSA fasta -> VCF from viridian nanopore:
cp /dev/null viridian_nanopore.fa
for id in $(virdidian_assemblies/Badread_assemblies/*); do
    echo ">$id" >> viridian_nanopore.fa
    tail -n+2 virdidian_assemblies/Badread_assemblies/$id.fasta >> viridian_nanopore.fa
    echo ""  >> viridian_nanopore.fa
done
time bash $mafftScript -i viridian_nanopore.fa -t $numThreads -o viridian_nanopore.msa.fa
./faToVcf -includeNoAltN <(cat $refFa viridian_nanopore.msa.fa) \
    viridian_nanopore.vcf

# Make monolithic MSA fasta -> VCF from artic illumina:
cp /dev/null artic_illimina.fa
for id in $(artic_assemblies/ART_assemblies/*); do
    echo ">$id" >> artic_illimina.fa
    tail -n+2 artic_assemblies/ART_assemblies/$id.fasta >> artic_illimina.fa
    echo ""  >> artic_illimina.fa
done
time bash $mafftScript -i artic_illimina.fa -t $numThreads -o artic_illimina.msa.fa
./faToVcf -includeNoAltN <(cat $refFa artic_illimina.msa.fa) \
    artic_illimina.vcf

# Make monolithic MSA fasta -> VCF from artic nanopore:
cp /dev/null artic_nanopore.fa
for id in $(artic_assemblies/Badread_assemblies/*); do
    echo ">$id" >> artic_nanopore.fa
    tail -n+2 artic_assemblies/Badread_assemblies/$id.fasta >> artic_nanopore.fa
    echo ""  >> artic_nanopore.fa
done
time bash $mafftScript -i artic_nanopore.fa -t $numThreads -o artic_nanopore.msa.fa
./faToVcf -includeNoAltN <(cat $refFa artic_nanopore.msa.fa) \
    artic_nanopore.vcf

# Make trees the impatient way: usher on whole VCF starting from empty tree.
for assembler in artic viridian; do
    for tech in ART_assemblies Badread_assemblies; do
        echo ${assembler}_${tech}
        time singularity exec ../singularity/usher/usher.sif usher --sort-before-placement-3 \
            --vcf ${assembler}_${tech}.vcf \
            --tree <(echo '()') \
            -T $numThreads \
            --save-mutation-annotated-tree ${assembler}_${tech}.pb \
            >& ${assembler}_${tech}.usher.log
    done
done

# Taxonium on unoptimized trees
for assembler in artic viridian; do
    for tech in ART_assemblies Badread_assemblies; do
        singularity exec ../singularity/usher/usher.sif matUtils mask -i ${assembler}_${tech}.pb -r corrected_map.tsv \
            -o ${assembler}_${tech}.renamed.pb >& renaming.out
        usher_to_taxonium --input ${assembler}_${tech}.renamed.pb \
        --genbank $taxoniumGb \
        --metadata metadata.tsv \
        --columns country,date,Nextstrain_clade,pangolin_lineage,sra_accession,pangolin_v4_lineage \
        --title "Daniel's faToVcf UShER (not optimized): $assembler $tech" \
        --output ${assembler}_${tech}.taxonium.jsonl.gz
    done
done

# Optimize
for assembler in artic viridian; do
    for tech in ART_assemblies Badread_assemblies; do
        time singularity exec ../singularity/usher/usher.sif matOptimize -T $numThreads -r 8 -M 2 \
            -i ${assembler}_${tech}.pb -o ${assembler}_${tech}.opt.pb \
            >& ${assembler}_${tech}.matOptimize.log
    done
done

# Taxonium on optimized trees
for assembler in artic viridian; do
    for tech in ART_assemblies Badread_assemblies; do
        singularity exec ~../singularity/usher/usher.sif matUtils mask -i ${assembler}_${tech}.opt.pb -r sraToGisaidFullName.tsv \
            -o ${assembler}_${tech}.opt.renamed.pb >& renaming.out
        usher_to_taxonium --input ${assembler}_${tech}.opt.renamed.pb \
        --genbank $taxoniumGb \
        --metadata metadata.tsv \
        --columns country,date,Nextstrain_clade,pangolin_lineage,sra_accession,pangolin_v4_lineage \
        --title "Daniel's faToVcf UShER (optimized): $assembler $tech" \
        --output ${assembler}_${tech}.opt.taxonium.jsonl.gz
    done
done

