version 1.0

task ScatterIntervalList {
  input {
    File interval_list
    Int scatter_count
    String picard_path = "/home/thais/bin/picard.jar"
  }

  command <<<
    set -e
    mkdir out
    java -Xms1g -jar ~{picard_path} \
      IntervalListTools \
      SCATTER_COUNT=~{scatter_count} \
      SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
      UNIQUE=true \
      SORT=true \
      INPUT=~{interval_list} \
      OUTPUT=out

    python3 <<CODE
    import glob, os
    # Works around a JES limitation where multiples files with the same name overwrite each other when globbed
    intervals = sorted(glob.glob("out/*/*.interval_list"))
    for i, interval in enumerate(intervals):
      (directory, filename) = os.path.split(interval)
      newName = os.path.join(directory, str(i + 1) + filename)
      os.rename(interval, newName)
    print(len(intervals))
    CODE
  >>>
  output {
    Array[File] out = glob("out/*/*.interval_list")
  }
  runtime {
    memory: "2 GiB"
  }
}

task HaplotypeCaller_GATK4_VCF {
  input {
    File input_bam
    File input_bam_index
    String vcf_basename
    File interval_list
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
      -L ~{interval_list} \
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

# Combine multiple VCFs or GVCFs from scattered HaplotypeCaller runs
task MergeVCFs {
  input {
    Array[File] input_vcfs
    Array[File] input_vcfs_indexes
    String sample_name
    String picard_path = "/home/thais/bin/picard.jar"
  }

  # Using MergeVcfs instead of GatherVcfs so we can create indices
  # See https://github.com/broadinstitute/picard/issues/789 for relevant GatherVcfs ticket
  command {
    java -Xms2000m -jar ~{picard_path} \
      MergeVcfs \
      INPUT=~{sep=' INPUT=' input_vcfs} \
      OUTPUT=~{sample_name}.g.vcf.gz
  }
  runtime {
    memory: "3 GiB"
  }
  output {
    File output_vcf = "~{sample_name}.g.vcf.gz"
    File output_vcf_index = "~{sample_name}.g.vcf.gz.tbi"
  }
}