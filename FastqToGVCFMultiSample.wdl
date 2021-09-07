version 1.0

import "FastqToGVCF.wdl"

workflow FastqToGVCFMultiSample {

  input {
    File sample_table

    Reference reference
    Array[File]? known_indels_sites_vcfs
    Array[File]? known_indels_sites_indices

    Reference humanReference

    # Vai aparecer apenas no BAM/CRAM final
    String library_name = "UnknownLibrary"
    String platform_name = "UnkownPlatform"

    Boolean is_pcr_free = false

    File calling_interval_list
    Int scatter_count = 12
  }

  scatter(row in read_tsv(sample_table)) {
    call FastqToGVCF.FastqToGVCF {
      input:
        sample_name = row[0],
        fastq_1 = row[1],
        fastq_2 = row[2],
        reference = reference,
        known_indels_sites_vcfs = known_indels_sites_vcfs,
        known_indels_sites_indices = known_indels_sites_indices,
        humanReference = humanReference,
        library_name = library_name,
        platform_name = platform_name,
        is_pcr_free = is_pcr_free,
        calling_interval_list = calling_interval_list,
        scatter_count = scatter_count
    }
  }

  output {
    Array[File] duplicate_metrics = FastqToGVCF.duplicate_metrics
    # Array[File] output_bqsr_reports = FastqToGVCF.output_bqsr_reports

    Array[File] output_bam = FastqToGVCF.output_bam
    Array[File] output_bam_index = FastqToGVCF.output_bam_index

    Array[File] output_vcf = FastqToGVCF.output_vcf
    Array[File] output_vcf_index = FastqToGVCF.output_vcf_index
  }
}