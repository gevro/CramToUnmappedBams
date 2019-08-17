## CRAM to Unmapped BAM workflow, July 2019
##
## This WDL pipeline combines the Broad GATK CRAM-to-BAM workflow with the BAM-to-Unmapped-BAM workflow.


import "https://raw.githubusercontent.com/gevro/CramToUnmappedBams/master/cram-to-bam.wdl" as CramToBam
import "https://raw.githubusercontent.com/gevro/CramToUnmappedBams/master/bam-to-unmapped-bams.wdl" as BamToUnmappedBams

# WORKFLOW DEFINITION
workflow CramtoUnmappedBam {
#CramToBam inputs
Int cram_to_bam_disk_size
Int validate_sam_file_disk_size
String cram_to_bam_mem_size
String validate_sam_file_mem_size
String? gotc_docker_override


#BamToUnmappedBams inputs   
Int? additional_disk_size
String? gatk_path
String? gitc_docker
String? gatk_docker 


call CramToBam.CramToBamFlow{
  input:
  cram_to_bam_disk_size = cram_to_bam_disk_size,
  validate_sam_file_disk_size = validate_sam_file_disk_size,
  cram_to_bam_mem_size = cram_to_bam_mem_size,
  validate_sam_file_mem_size = validate_sam_file_mem_size,
  gotc_docker_override = gotc_docker_override,
}

call BamToUnmappedBams.BamToUnmappedBams{
  input:
  input_bam = CramToBamFlow.outputBam,
  additional_disk_size = additional_disk_size,
  gatk_path = gatk_path,
  gitc_docker = gitc_docker,
  gatk_docker = gatk_docker,
}

call CreateFoFN {
  input:
    array_of_files = BamToUnmappedBams.output_bams,
    docker = "us.gcr.io/broad-gatk/gatk:latest"
}

output {
  File Bam_validation_report = CramToBamFlow.validation_report
  Array[File] output_unmappedbams = BamToUnmappedBams.output_bams
  File unmapped_bam_list = CreateFoFN.fofn_list
}
  
}

task CreateFoFN {
  # Command parameters
  String SampleName
  Array[String] array_of_files
  String fofn_name = SampleName + "unmapped.bam.list"
  
  # Runtime parameters
  String docker
  
  command {
    mv ${write_lines(array_of_files)} ${fofn_name}
  }
  output {
    File fofn_list = "${fofn_name}"
  }
  runtime {
    docker: docker
    preemptible: 3
  }
}
