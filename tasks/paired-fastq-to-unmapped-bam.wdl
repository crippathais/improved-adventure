version 1.0

task PairedFastQsToUnmappedBAM {
  input {
    String sample_name
    File fastq_1
    File fastq_2
    String readgroup_name
    String library_name
    String platform_name
    Int machine_mem_gb = 7
  }
  
  command {
    /home/thais/bin/gatk/gatk --java-options "-Xmx~{machine_mem_gb - 1}g" \
    FastqToSam \
    --FASTQ ~{fastq_1} \
    --FASTQ2 ~{fastq_2} \
    --OUTPUT ~{readgroup_name}.unmapped.bam \
    --READ_GROUP_NAME ~{readgroup_name} \
    --SAMPLE_NAME ~{sample_name} \
    --LIBRARY_NAME ~{library_name} \
    --PLATFORM ~{platform_name}
  }
  runtime {
    memory: machine_mem_gb + " GB"
  }
  output {
    File output_unmapped_bam = "~{readgroup_name}.unmapped.bam"
  }
}