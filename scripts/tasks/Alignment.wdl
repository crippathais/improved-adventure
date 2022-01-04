version 1.0

import "DNASeqStructs.wdl"

# Read unmapped BAM, convert on-the-fly to FASTQ and stream to BWA MEM for alignment, then stream to MergeBamAlignment
task SamToFastqAndBwaMemAndMba {
  input {
    File unmapped_bam
    String bwa_commandline = "bwa mem -K 100000000 -M -p -v 3 -t 16 -Y $bash_ref_fasta"
    String output_bam_basename

    Reference reference

    Int compression_level = 2
    Boolean hard_clip_reads = false
  }

  command <<<
    # This is done before "set -o pipefail" because "bwa" will have a rc=1 and we don't want to allow rc=1 to succeed
    # because the sed may also fail with that error and that is something we actually want to fail on.
    BWA_VERSION=$(/home/thais/bin/bwa 2>&1 | \
    grep -e '^Version' | \
    sed 's/Version: //')

    set -o pipefail
    set -e

    if [ -z ${BWA_VERSION} ]; then
        exit 1;
    fi

    # set the bash variable needed for the command-line
    bash_ref_fasta=~{reference.ref_fasta}
    java -Xms1000m -Xmx1000m -jar /home/thais/bin/picard.jar \
      SamToFastq \
      INPUT=~{unmapped_bam} \
      FASTQ=/dev/stdout \
      INTERLEAVE=true \
      NON_PF=true | \
    /home/thais/bin/~{bwa_commandline} /dev/stdin - 2> >(tee ~{output_bam_basename}.bwa.stderr.log >&2) | \
    java -Dsamjdk.compression_level=~{compression_level} -Xms1000m -Xmx1000m -jar /home/thais/bin/picard.jar \
      MergeBamAlignment \
      VALIDATION_STRINGENCY=SILENT \
      EXPECTED_ORIENTATIONS=FR \
      ATTRIBUTES_TO_RETAIN=X0 \
      ATTRIBUTES_TO_REMOVE=NM \
      ATTRIBUTES_TO_REMOVE=MD \
      ALIGNED_BAM=/dev/stdin \
      UNMAPPED_BAM=~{unmapped_bam} \
      OUTPUT=~{output_bam_basename}.bam \
      REFERENCE_SEQUENCE=~{reference.ref_fasta} \
      SORT_ORDER="unsorted" \
      IS_BISULFITE_SEQUENCE=false \
      ALIGNED_READS_ONLY=false \
      CLIP_ADAPTERS=false \
      ~{true='HARD_CLIP_OVERLAPPING_READS=true' false="" hard_clip_reads} \
      ~{true='CLIP_OVERLAPPING_READS=true' false="" hard_clip_reads} \
      MAX_RECORDS_IN_RAM=2000000 \
      ADD_MATE_CIGAR=true \
      MAX_INSERTIONS_OR_DELETIONS=-1 \
      PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
      PROGRAM_RECORD_ID="bwamem" \
      PROGRAM_GROUP_VERSION="${BWA_VERSION}" \
      PROGRAM_GROUP_COMMAND_LINE="~{bwa_commandline}" \
      PROGRAM_GROUP_NAME="bwamem" \
      UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
      ALIGNER_PROPER_PAIR_FLAGS=true \
      ADD_PG_TAG_TO_READS=false
  >>>
  runtime {
    #docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710"
    memory: "14 GiB"
    cpu: "16"
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
    File bwa_stderr_log = "~{output_bam_basename}.bwa.stderr.log"
  }
}

task ExtractUnmappedReads {
  input {
    String sample_name
    File aligned_bam
  }

  command <<<
    /home/linuxbrew/.linuxbrew/bin/samtools view -b -f 12 ~{aligned_bam} > ~{sample_name}.unmapped.bam
  >>>

  runtime {
    #docker: "quay.io/biocontainers/samtools:1.2--0"
  }

  output {
    File unmapped_bam = "~{sample_name}.unmapped.bam"
  }
}