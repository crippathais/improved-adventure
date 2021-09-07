version 1.0

# Sort BAM file by coordinate order
task SortSam {
  input {
    File input_bam
    String output_bam_basename
    Int compression_level = 2
  }
  # SortSam spills to disk a lot more because we are only store 300000 records in RAM now because its faster for our data so it needs
  # more disk space.  Also it spills to disk in an uncompressed format so we need to account for that with a larger multiplier

  command {
    java -Dsamjdk.compression_level=~{compression_level} -Xms4000m -jar /home/thais/bin/picard.jar \
      SortSam \
      INPUT=~{input_bam} \
      OUTPUT=~{output_bam_basename}.bam \
      SORT_ORDER="coordinate" \
      CREATE_INDEX=true \
      CREATE_MD5_FILE=true \
      MAX_RECORDS_IN_RAM=300000

  }
  runtime {
    ## "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
    cpu: "1"
    memory: "5000 MiB"
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
    File output_bam_index = "~{output_bam_basename}.bai"
    File output_bam_md5 = "~{output_bam_basename}.bam.md5"
  }
}

# Mark duplicate reads to avoid counting non-independent observations
task MarkDuplicates {
  input {
    Array[File] input_bams
    String output_bam_basename
    String metrics_filename
    Int compression_level = 2

    Int memory_multiplier = 1
    Float? sorting_collection_size_ratio
  }

  Float memory_size = 7.5 * memory_multiplier
  Int java_memory_size = (ceil(memory_size) - 2)

  # Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly
  # This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
  # While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"

  command {
    java -Dsamjdk.compression_level=~{compression_level} -Xms~{java_memory_size}g -jar /home/thais/bin/picard.jar \
      MarkDuplicates \
      INPUT=~{sep=' INPUT=' input_bams} \
      OUTPUT=~{output_bam_basename}.bam \
      METRICS_FILE=~{metrics_filename} \
      VALIDATION_STRINGENCY=SILENT \
      ~{"SORTING_COLLECTION_SIZE_RATIO=" + sorting_collection_size_ratio} \
      OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
      ASSUME_SORT_ORDER="queryname" \
      CLEAR_DT="false" \
      ADD_PG_TAG_TO_READS=false
  }
  runtime {
    #docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
    memory: "~{memory_size} GiB"
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
    File duplicate_metrics = "~{metrics_filename}"
  }
}

# Generate Base Quality Score Recalibration (BQSR) model
task BaseRecalibrator {
  input {
    File input_bam
    File input_bam_index
    String recalibration_report_filename
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Array[File]? known_indels_sites_vcfs
    Array[File]? known_indels_sites_indices
  }

  command {
    gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
      -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
      -Xloggc:gc_log.log -Xms5g" \
      BaseRecalibrator \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      --use-original-qualities \
      ~{true="--known-sites" false="" defined(known_indels_sites_vcfs)} ~{sep=" --known-sites " known_indels_sites_vcfs} \
      -O ~{recalibration_report_filename}
  }
  runtime {
    memory: "6 GiB"
  }
  output {
    File recalibration_report = "~{recalibration_report_filename}"
  }
}

# Apply Base Quality Score Recalibration (BQSR) model
task ApplyBQSR {
  input {
    File input_bam
    File input_bam_index
    String output_bam_basename
    File recalibration_report
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Int compression_level = 2
    Int memory_multiplier = 1
  }

  Int memory_size = ceil(3500 * memory_multiplier)

  command {
    gatk --java-options "-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
      -XX:+PrintGCDetails -Xloggc:gc_log.log \
      -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Dsamjdk.compression_level=~{compression_level} -Xms3000m" \
      ApplyBQSR \
      --create-output-bam-md5 \
      --add-output-sam-program-record \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      --use-original-qualities \
      -O ~{output_bam_basename}.bam \
      -bqsr ~{recalibration_report} \
      --create-output-bam-index
  }
  runtime {
    memory: "~{memory_size} MiB"
  }
  output {
    File recalibrated_bam = "~{output_bam_basename}.bam"
    File recalibrated_bam_index = "~{output_bam_basename}.bai"
    File recalibrated_bam_checksum = "~{output_bam_basename}.bam.md5"
  }
}
