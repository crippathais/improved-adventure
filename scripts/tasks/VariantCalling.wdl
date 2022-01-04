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
    File calling_interval_list
    Int scatter_count
  }

  call Calling.ScatterIntervalList as ScatterIntervalList {
    input:
      interval_list = calling_interval_list,
      scatter_count = scatter_count
  }

  scatter (scattered_interval_list in ScatterIntervalList.out) {
    call Calling.HaplotypeCaller_GATK4_VCF as HaplotypeCallerGATK4 {
      input:
        input_bam = input_bam,
        input_bam_index = input_bam_index,
        vcf_basename = sample_name,
        interval_list = scattered_interval_list,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        is_pcr_free = is_pcr_free
    }
  }

  call Calling.MergeVCFs as MergeVCFs {
    input:
      input_vcfs = HaplotypeCallerGATK4.output_vcf,
      input_vcfs_indexes = HaplotypeCallerGATK4.output_vcf_index,
      sample_name = sample_name
  }

  output {
    File output_vcf = MergeVCFs.output_vcf
    File output_vcf_index = MergeVCFs.output_vcf_index
  }

  meta {
    allowNestedInputs: true
  }
}
