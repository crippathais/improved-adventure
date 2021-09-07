version 1.0

import "GermlineVariantDiscovery.wdl" as Calling
import "BamProcessing.wdl" as BamProcessing

workflow VariantCalling {

  String pipeline_version = "1.0.1"

  input {
    String sample_name
    File input_bam
    File input_bam_index
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    Boolean is_pcr_free
  }

  call Calling.HaplotypeCaller_GATK4_VCF as HaplotypeCallerGATK4 {
    input:
      input_bam = input_bam,
      input_bam_index = input_bam_index,
      vcf_basename = sample_name,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      is_pcr_free = is_pcr_free
  }

  output {
    File output_vcf = HaplotypeCallerGATK4.output_vcf
    File output_vcf_index = HaplotypeCallerGATK4.output_vcf_index
  }

  meta {
    allowNestedInputs: true
  }
}

