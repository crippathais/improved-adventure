# SNP
/home/thais/bin/gatk/gatk SelectVariants \
    -V /home/thais/NewData/resultados/Padawan.vcf.gz \
    -select-type SNP \
    -O /home/thais/NewData/resultados/Padawan.snps.vcf.gz

/home/thais/bin/gatk/gatk VariantFiltration \
    -V /home/thais/NewData/resultados/Padawan.snps.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O /home/thais/NewData/resultados/Jedi.snps.vcf.gz

# Indel
/home/thais/bin/gatk/gatk SelectVariants \
    -V /home/thais/NewData/resultados/Padawan.vcf.gz \
    -select-type INDEL \
    -O /home/thais/NewData/resultados/Padawan.indels.vcf.gz


/home/thais/bin/gatk/gatk VariantFiltration \
-V /home/thais/NewData/resultados/Padawan.indels.vcf.gz \
-filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    -O /home/thais/NewData/resultados/Jedi.indels.vcf.gz