version 1.0

import "tasks/UnmappedBamToAlignedBam.wdl" as ToBam
import "tasks/VariantCalling.wdl" as ToGvcf
import "tasks/paired-fastq-to-unmapped-bam.wdl" as FastqToUbam
import "tasks/DNASeqStructs.wdl"
import "tasks/Alignment.wdl" as Alignment

# Esse pipeline processa uma amostra por vez. Recebe um par de arquivos FASTQ
# e reporta o BAM final (analysis-ready) e gVCF.
workflow FastqToGVCF {

  String pipeline_version = "2.3.6"

  input {
    String sample_name
    File fastq_1
    File fastq_2

    Reference reference
    Array[File]? known_indels_sites_vcfs
    Array[File]? known_indels_sites_indices

    # Vai aparecer apenas no BAM/CRAM final
    String library_name = "UnknownLibrary"
    String platform_name = "UnkownPlatform"

    Boolean is_pcr_free = false

    File calling_interval_list
    Int scatter_count = 12
  }

  # Converter o par de FASTQ para uBAM necessário para adicionar informação de
  # RG, LB, e PL
  call FastqToUbam.PairedFastQsToUnmappedBAM {
    input:
      sample_name = sample_name,
      fastq_1 = fastq_1,
      fastq_2 = fastq_2,
      # passando sample_name como RG pq é apenas um par de FASTQ por amostra
      readgroup_name = sample_name, 
      platform_name = platform_name,
      library_name = library_name
  }

  # uBAM é usado como entrada do subworkflow de alinhamento e refinamento do BAM
  call ToBam.UnmappedBamToAlignedBam {
    input:
      sample_name = sample_name,
      unmapped_bam = PairedFastQsToUnmappedBAM.output_unmapped_bam, #ExtractUnmappedReads.unmapped_bam, 
      reference = reference
  }

  # O BAM final é usado para fazer chamada de variantes
  call ToGvcf.VariantCalling as BamToGvcf {
    input:
      sample_name = sample_name,
      input_bam = UnmappedBamToAlignedBam.output_bam,
      input_bam_index = UnmappedBamToAlignedBam.output_bam_index,
      ref_fasta = reference.ref_fasta,
      ref_fasta_index = reference.ref_fasta_index,
      ref_dict = reference.ref_dict,
      is_pcr_free = is_pcr_free,
      calling_interval_list = calling_interval_list,
      scatter_count = scatter_count
  }

  output {
    File duplicate_metrics = UnmappedBamToAlignedBam.duplicate_metrics
    # File output_bqsr_reports = UnmappedBamToAlignedBam.output_bqsr_reports

    File output_bam = UnmappedBamToAlignedBam.output_bam
    File output_bam_index = UnmappedBamToAlignedBam.output_bam_index

    File output_vcf = BamToGvcf.output_vcf
    File output_vcf_index = BamToGvcf.output_vcf_index
  }
  meta {
    allowNestedInputs: true
  }
}
