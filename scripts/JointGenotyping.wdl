version 1.0

import "tasks/JointGenotypingTasks.wdl" as Tasks


workflow JointGenotyping {

  String pipeline_version = "1.5.1"

  input {
    File unpadded_intervals_file

    String callset_name
    File sample_name_map

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    Int scatter_count = 12
  }

  call Tasks.SplitIntervalList {
    input:
      interval_list = unpadded_intervals_file,
      scatter_count = scatter_count,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict
  }

  Array[File] unpadded_intervals = SplitIntervalList.output_intervals

  scatter (idx in range(length(unpadded_intervals))) {
    # The batch_size value was carefully chosen here as it
    # is the optimal value for the amount of memory allocated
    # within the task; please do not change it without consulting
    # the Hellbender (GATK engine) team!
    call Tasks.ImportGVCFs {
      input:
        sample_name_map = sample_name_map,
        interval = unpadded_intervals[idx],
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        workspace_dir_name = "genomicsdb",
        batch_size = 50
    }

    call Tasks.GenotypeGVCFs {
      input:
        workspace_tar = ImportGVCFs.output_genomicsdb,
        interval = unpadded_intervals[idx],
        output_vcf_filename = callset_name + "." + idx + ".vcf.gz",
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict
    }
  }

  call Tasks.GatherVcfs as FinalGatherVcf {
    input:
      input_vcfs = GenotypeGVCFs.output_vcf,
      output_vcf_name = callset_name + ".vcf.gz"
  }

  output {
    File output_vcfs = FinalGatherVcf.output_vcf
    File output_vcf_indices = FinalGatherVcf.output_vcf_index
  }
  meta {
    allowNestedInputs: true
  }
}
