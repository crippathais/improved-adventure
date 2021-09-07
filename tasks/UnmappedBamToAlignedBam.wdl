version 1.0

import "Alignment.wdl" as Alignment
import "BamProcessing.wdl" as Processing
import "DNASeqStructs.wdl"

# WORKFLOW DEFINITION
workflow UnmappedBamToAlignedBam {

  input {
    String sample_name
    File unmapped_bam
    Reference reference
    Array[File]? known_indels_sites_vcfs
    Array[File]? known_indels_sites_indices
  }

  # Map reads to reference
  call Alignment.SamToFastqAndBwaMemAndMba as SamToFastqAndBwaMemAndMba {
    input:
      unmapped_bam = unmapped_bam,
      output_bam_basename = sample_name + ".aligned.unsorted",
      reference = reference,
      hard_clip_reads = true #align the animal
  }

  # Aggregate aligned+merged flowcell BAM files and mark duplicates
  # We take advantage of the tool's ability to take multiple BAM inputs and write out a single output
  # to avoid having to spend time just merging BAM files.
  call Processing.MarkDuplicates as MarkDuplicates {
    input:
      input_bams = [SamToFastqAndBwaMemAndMba.output_bam],
      output_bam_basename = sample_name + ".aligned.unsorted.duplicates_marked",
      metrics_filename = sample_name + ".duplicate_metrics"
  }

  # Sort aggregated+deduped BAM file and fix tags
  call Processing.SortSam as SortSampleBam {
    input:
      input_bam = MarkDuplicates.output_bam,
      output_bam_basename = sample_name + ".aligned.duplicate_marked.sorted"
  }

  # call Processing.BaseRecalibrator as BaseRecalibrator {
  #   input:
  #     input_bam = SortSampleBam.output_bam,
  #     input_bam_index = SortSampleBam.output_bam_index,
  #     recalibration_report_filename = sample_name + ".recal_data.csv",
  #     ref_dict = reference.ref_dict,
  #     ref_fasta = reference.ref_fasta,
  #     ref_fasta_index = reference.ref_fasta_index,
  #     known_indels_sites_vcfs = known_indels_sites_vcfs,
  #     known_indels_sites_indices = known_indels_sites_indices
  # }

  # call Processing.ApplyBQSR as ApplyBQSR {
  #   input:
  #     input_bam = SortSampleBam.output_bam,
  #     input_bam_index = SortSampleBam.output_bam_index,
  #     output_bam_basename = sample_name,
  #     recalibration_report = BaseRecalibrator.recalibration_report,
  #     ref_dict = reference.ref_dict,
  #     ref_fasta = reference.ref_fasta,
  #     ref_fasta_index = reference.ref_fasta_index
  # }

  # Outputs that will be retained when execution is complete
  output {
    File duplicate_metrics = MarkDuplicates.duplicate_metrics
    #File output_bqsr_reports = BaseRecalibrator.recalibration_report

    File output_bam = SortSampleBam.output_bam
    File output_bam_index = SortSampleBam.output_bam_index
  }
  meta {
    allowNestedInputs: true
  }
}
