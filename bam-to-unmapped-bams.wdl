## Copyright Broad Institute, 2018
## 
## This WDL converts BAM  to unmapped BAMs
##
## Requirements/expectations :
## - BAM file
##
## Outputs :
## - Sorted Unmapped BAMs
##
## Cromwell version support
## - Successfully tested on v33
## - Does not work on versions < v23 due to output syntax
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation. 
## For program versions, see docker containers. 
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may 
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker 
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

# WORKFLOW DEFINITION
workflow BamToUnmappedBams {
  File input_bam

  Int? additional_disk_size
  Int additional_disk = select_first([additional_disk_size, 20])

  Float input_size = size(input_bam, "GB")
  
  String? gatk_path
  String path2gatk = select_first([gatk_path, "/gatk/gatk"])

  String? gitc_docker
  String gitc_image = select_first([gitc_docker, "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"])
  String? gatk_docker 
  String gatk_image = select_first([gatk_docker, "broadinstitute/gatk:latest"])

  String R1_adapter
  String R2_adapter
  Int trimstart
  Int minadapteroverlap

  String sample_name

  call GenerateOutputMap {
    input:
      input_bam = input_bam,
      disk_size = ceil(input_size) + additional_disk,
      docker = gitc_image
  }

  call RevertSam {
    input:
      input_bam = input_bam,
      output_map = GenerateOutputMap.output_map,
      disk_size = ceil(input_size * 3) + additional_disk,
      docker = gatk_image,
      gatk_path = path2gatk
  }

  scatter (unmapped_bam in RevertSam.unmapped_bams) {
    String output_basename = basename(unmapped_bam, ".coord.sorted.unmapped.bam")
    Float unmapped_bam_size = size(unmapped_bam, "GB")

    call SortSam {
      input:
        input_bam = unmapped_bam,
        sorted_bam_name = output_basename + ".unmapped.bam",
        disk_size = ceil(unmapped_bam_size * 10) + additional_disk,
        docker = gatk_image,
        gatk_path = path2gatk
    }

    call SamToFastq{
      input:
        input_bam = SortSam.sorted_bam
    }

    call CutAdapt{
      input:
        fastq_1 = SamToFastq.fastq_1,
        fastq_2 = SamToFastq.fastq_2,
        R1_adapter = R1_adapter,
        R2_adapter = R2_adapter,
        trimstart = trimstart,
        minadapteroverlap = minadapteroverlap,
        readgroup_name = output_basename
    }

    call PairedFastQsToUnmappedBAM{
      input:
        fastq_1 = CutAdapt.fastq1_trimmed,
        fastq_2 = CutAdapt.fastq2_trimmed,
        input_bam = SortSam.sorted_bam,
        sample_name = sample_name,
        readgroup_name = output_basename
    }
  }

  output {
    Array[File] output_bams = PairedFastQsToUnmappedBAM.output_bam
    Array[File] CutAdapt_output = CutAdapt.CutAdapt_output
  }
}

task GenerateOutputMap {
  File input_bam
  Int disk_size
  
  String docker

  command {
    set -e

    samtools view -H ${input_bam} | grep @RG | cut -f2 | sed s/ID:// > readgroups.txt

    echo -e "READ_GROUP_ID\tOUTPUT" > output_map.tsv

    for rg in `cat readgroups.txt`; do
      echo -e "$rg\t$rg.coord.sorted.unmapped.bam" >> output_map.tsv
    done
  }

  runtime {
    docker: docker
    disks: "local-disk " + disk_size + " HDD"
    preemptible: "3"
    memory: "1 GB"
  }
  output {
    File output_map = "output_map.tsv"
  }
}

task RevertSam {
  File input_bam
  File output_map
  Int disk_size

  String gatk_path

  String docker

  command {
    ${gatk_path} --java-options "-Xmx1000m" \
    RevertSam \
    --INPUT ${input_bam} \
    --OUTPUT_MAP ${output_map} \
    --OUTPUT_BY_READGROUP true \
    --VALIDATION_STRINGENCY LENIENT \
    --ATTRIBUTE_TO_CLEAR FT \
    --ATTRIBUTE_TO_CLEAR CO \
    --SORT_ORDER coordinate
  }
  runtime {
    docker: docker
    disks: "local-disk " + disk_size + " HDD"
    memory: "1200 MB"
  }
  output {
    Array[File] unmapped_bams = glob("*.bam")
  }
}

task SortSam {
  File input_bam
  String sorted_bam_name
  Int disk_size

  String gatk_path

  String docker

  command {
    ${gatk_path} --java-options "-Xmx18000m" \
    SortSam \
    --INPUT ${input_bam} \
    --OUTPUT ${sorted_bam_name} \
    --SORT_ORDER queryname \
    --MAX_RECORDS_IN_RAM 1000000
  }
  runtime {
    docker: docker
    disks: "local-disk " + disk_size + " HDD"
    memory: "16000 MB"
    preemptible: 0
  }
  output {
    File sorted_bam = "${sorted_bam_name}"
  }
}


#Convert Sam to Fastq for Cutadapt processing
task SamToFastq {
  # Command parameters
  File input_bam
  String output_basename = basename(input_bam)
  
  command {
    samtools fastq -@ 1 -n -1 ${output_basename}.R1.fastq.gz -2 ${output_basename}.R2.fastq.gz ${input_bam}
  }
  runtime {
    docker: "halllab/samtools:v1.9"
    memory: "4 GB"
    cpu: "2"
    disks: "local-disk "+ ceil(size(input_bam,"GB")*5) + " HDD"
    preemptible: 0
  }
  output {
    File fastq_1 = "${output_basename}.R1.fastq.gz"
    File fastq_2 = "${output_basename}.R2.fastq.gz"
  }
}


#Run CutAdapt to trim Illumina Adapaters
task CutAdapt {
  # Command parameters
  File fastq_1
  File fastq_2
  String R1_adapter
  String R2_adapter
  Int trimstart
  Int minadapteroverlap
  String file_output1 = basename(fastq_1,".fastq.gz") + ".trimmed.fastq.gz"
  String file_output2 = basename(fastq_2,".fastq.gz") + ".trimmed.fastq.gz"
  String readgroup_name
  String min_length_to_keep = "30"

  command {
    cutadapt -j 4 -u ${trimstart} -U ${trimstart} -O ${minadapteroverlap} \
    -m ${min_length_to_keep} -a ${R1_adapter} -A ${R2_adapter} \
    -o ${file_output1} -p ${file_output2} ${fastq_1} ${fastq_2} > ${readgroup_name}.cutadapt.out
  }
  runtime {
    docker: "kfdrc/cutadapt:latest"
    memory: "4 GB"
    cpu: "4"
    disks: "local-disk "+ ceil(size(fastq_1,"GB")*5) + " HDD"
    preemptible: 0
  }
  output {
    File fastq1_trimmed = "${file_output1}"
    File fastq2_trimmed = "${file_output2}"
    File CutAdapt_output = "${readgroup_name}.cutadapt.out"
  }
}

# Convert a pair of FASTQs to uBAM, and get readgroup, etc information from the original ubam prior to adapter trimming
task PairedFastQsToUnmappedBAM {
  # Command parameters
  File fastq_1
  File fastq_2
  File input_bam
  String readgroup_name
  String sample_name

  command {
    /gatk/gatk --java-options "-Xmx3000m" \
    FastqToSam \
    --FASTQ ${fastq_1} \
    --FASTQ2 ${fastq_2} \
    --OUTPUT ${readgroup_name}.unmapped.bam \
    --READ_GROUP_NAME `samtools view -H ${input_bam} | grep ^@RG | sed 's/.*ID://;s/\t.*//'` \
    --SAMPLE_NAME ${sample_name} \
    --LIBRARY_NAME `samtools view -H ${input_bam} | grep ^@RG | sed 's/.*LB://;s/\t.*//'` \
    --PLATFORM `samtools view -H ${input_bam} | grep ^@RG | sed 's/.*PL://;s/\t.*//'` \
    --SEQUENCING_CENTER `samtools view -H ${input_bam} | grep ^@RG | sed 's/.*CN://;s/\t.*//'`
  }
  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:latest"
    memory: "10 GB"
    cpu: "1"
    disks: "local-disk "+ ceil(size(fastq_1,"GB")*25) + " HDD"
    preemptible: 0
  }
  output {
    File output_bam = "${readgroup_name}.unmapped.bam"
  }
}