version 1.0

task HaplotypeCaller_GATK4_VCF {
  input {
    File input_bam
    File input_bam_index
    String vcf_basename
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Boolean is_pcr_free

  }

  String output_file_name = vcf_basename + ".g.vcf.gz"

  command <<<
    set -e
    /home/thais/bin/gatk/gatk --java-options "-Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
      HaplotypeCaller \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      -O ~{output_file_name} \
      -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation \
      -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
      ~{true="--pcr-indel-model NONE" false="" is_pcr_free} \
      -ERC GVCF
  >>>

  runtime {
    memory: "6.5 GiB"
    cpu: "2"
  }

  output {
    File output_vcf = "~{output_file_name}"
    File output_vcf_index = "~{output_file_name}.tbi"
  }
}
