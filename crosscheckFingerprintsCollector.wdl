version 1.0

import "imports/pull_bwaMem.wdl" as bwaMem
import "imports/pull_star.wdl" as star

struct InputGroup {
  File fastqR1
  File fastqR2
  String readGroup
}

workflow crosscheckFingerprintsCollector {
   input {
        File? fastqR1
        File? fastqR2
        File? bam
        File? bamIndex
        String inputType
        String aligner
        Boolean markDups
        String outputFileNamePrefix
        String refFasta
        String haplotypeMap
        Int maxReads = 0
        String sampleId
   }
   parameter_meta {
        fastqR1: "fastq file for read 1"
        fastqR2: "fastq file for read 2"
        bam: "bam file"
        bamIndex: "bam index file"
        inputType: "one of either fastq or bam"
        aligner : "aligner to use for fastq input, either bwa or star"
        markDups : "should the alignment be duplicate marked?, generally yes"
        outputFileNamePrefix: "Optional output prefix for the output"
        refFasta: "Path to the reference fasta file"
        haplotypeMap: "Path to the gzipped hotspot vcf file"
        maxReads: "The maximum number of reads to process; if set, this will sample the requested number of reads"
        sampleId : "value that will be used as the sample identifier in the vcf fingerprint"
   }

   if(inputType=="fastq" && defined(fastqR1) && defined(fastqR2)){
     if(maxReads>0){
      call downsample{
        input:
          fastqR1 = select_first([fastqR1]),
          fastqR2 = select_first([fastqR2]),
          maxReads = maxReads
      }
     }

     if(aligner=="bwa"){
       call bwaMem.bwaMem {
         input:
           fastqR1 = select_first([downsample.fastqR1mod,fastqR1]),
           fastqR2 = select_first([downsample.fastqR2mod,fastqR2]),
           outputFileNamePrefix = outputFileNamePrefix,
           readGroups = "'@RG\\tID:CROSSCHECK\\tSM:SAMPLE'",
           doTrim = false
        }
      }

      if(aligner=="star"){
       InputGroup starInput = { 
         "fastqR1": select_first([downsample.fastqR1mod,fastqR1]),
         "fastqR2": select_first([downsample.fastqR2mod,fastqR2]),
         "readGroup": "ID:CROSSCHECK SM:SAMPLE"
       }
       call star.star { 
         input:
           inputGroups = [ starInput ],
           outputFileNamePrefix = outputFileNamePrefix,
           runStar_chimOutType = "Junctions"
       }
     }
   }

  if(markDups){
    call markDuplicates {
      input :
        inputBam = select_first([bwaMem.bwaMemBam,star.starBam,bam]),
        inputBai = select_first([bwaMem.bwaMemIndex,star.starIndex,bamIndex]),
        outputFileNamePrefix = outputFileNamePrefix 
    }
  }  

   call alignmentMetrics {
     input:
        inputBam = select_first([markDuplicates.bam,bwaMem.bwaMemBam,star.starBam,bam]),
        inputBai = select_first([markDuplicates.bamIndex,bwaMem.bwaMemIndex,star.starIndex,bamIndex]),
        outputFileNamePrefix = outputFileNamePrefix,
        markDups = markDups  
   }
   
   call extractFingerprint {
     input:
        inputBam = select_first([markDuplicates.bam,bwaMem.bwaMemBam,star.starBam,bam]),
        inputBai = select_first([markDuplicates.bamIndex,bwaMem.bwaMemIndex,star.starIndex,bamIndex]),
        haplotypeMap = haplotypeMap,
        refFasta = refFasta,
        outputFileNamePrefix = outputFileNamePrefix,
        sampleId = sampleId
    }

   output {
      File outputVcf = extractFingerprint.vgz
      File outputTbi = extractFingerprint.tbi
      File json = alignmentMetrics.json
      File samstats = alignmentMetrics.samstats
     }

    meta {
     author: "Lawrence Heisler"
     email: "lawrence.heisler@oicr.on.ca"
     description: "crosscheckFingerprintsCollector, workflow that generates genotype fingerprints using gatk ExtractFingprint.  Output are vcf files that can be proccessed through gatk Crosscheck fingerprints\n##"
     dependencies: [
      {
        name: "gatk/4.2.0.0",
        url: "https://gatk.broadinstitute.org"
      },
      {
        name: "tabix/0.2.6",
        url: "http://www.htslib.org"
      }
     ]
     output_meta: {
       outputVcf: "the crosscheck fingerprint, gzipped vcf file",
       outputTbi: "index for the vcf fingerprint",
       coverage : "output from samtools coverage, with per chromosome metrics",
       json : "metrics in json format, currently only the mean coverage for the alignment",
       samstats : "output from the samstats summary"
     }
  }

}
# ==========================================
#  configure and run extractFingerprintsCollector
# ==========================================


task extractFingerprint {
input {
 File inputBam
 File inputBai
 String modules
 String refFasta
 String outputFileNamePrefix
 String haplotypeMap
 String sampleId
 Int jobMemory = 8
 Int timeout = 24
}
parameter_meta {
 inputBam: "input .bam file"
 inputBai: "index of the input .bam file"
 refFasta: "Path to reference FASTA file"
 outputFileNamePrefix: "prefix for making names for output files"
 haplotypeMap: "Hotspot SNPs are the locations of variants used for genotyping"
 sampleId : "value used as the sample identifier in the vcf fingerprint"
 jobMemory: "memory allocated for Job"
 modules: "Names and versions of modules"
 timeout: "Timeout in hours, needed to override imposed limits"
}

command <<<
  set -euo pipefail

 $GATK_ROOT/bin/gatk ExtractFingerprint \
                    -R ~{refFasta} \
                    -H ~{haplotypeMap} \
                    -I ~{inputBam} \
                    -O ~{outputFileNamePrefix}.vcf \
                    --SAMPLE_ALIAS ~{sampleId} \
                    --VALIDATION_STRINGENCY LENIENT

 $TABIX_ROOT/bin/bgzip -c ~{outputFileNamePrefix}.vcf > ~{outputFileNamePrefix}.vcf.gz
 $TABIX_ROOT/bin/tabix -p vcf ~{outputFileNamePrefix}.vcf.gz 
>>>

 runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
 }

 output {
  File vcf = "~{outputFileNamePrefix}.vcf"
  File vgz = "~{outputFileNamePrefix}.vcf.gz"
  File tbi = "~{outputFileNamePrefix}.vcf.gz.tbi"
 }
}



# ==========================================
#  downsample the fastq files to the first N reads
# ==========================================


task downsample {
 input{
  File fastqR1
  File fastqR2
  Int maxReads
  String modules
  Int jobMemory = 8
  Int timeout = 24
 }
 
 parameter_meta {
  fastqR1 : "Read1 fastq file"
  fastqR2 : "Read2 fastq file"
  maxReads : "the maximum number of reads to use"
  jobMemory: "memory allocated for Job"
  modules: "Names and versions of modules"
  timeout: "Timeout in hours, needed to override imposed limits"
 }
 
 String fastqR1m = basename(fastqR1,".fastq.gz") + ".mod.fastq"
 String fastqR2m = basename(fastqR2,".fastq.gz") + ".mod.fastq"
  
command <<<
 set -euo pipefail
 
 seqtk sample -s 100 ~{fastqR1} ~{maxReads} > ~{fastqR1m}
 gzip ~{fastqR1m}
 
 seqtk sample -s 100 ~{fastqR2} ~{maxReads} > ~{fastqR2m}
 gzip ~{fastqR2m}
>>>

 runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
 }

 output {
  File fastqR1mod = "~{fastqR1m}.gz"
  File fastqR2mod = "~{fastqR2m}.gz"
 }    
}


# ==========================================
#  Duplicate Marking
# ==========================================



task markDuplicates {
 input{
  File inputBam
  File inputBai
  String modules
  String outputFileNamePrefix
  Int jobMemory = 8
  Int timeout = 24
 }
 parameter_meta {
  inputBam: "input .bam file"
  inputBai: "index of the input .bam file"
  outputFileNamePrefix: "prefix for making names for output files"  
  jobMemory: "memory allocated for Job"
  modules: "Names and versions of modules"
  timeout: "Timeout in hours, needed to override imposed limits"
 }
 
command <<<
  set -euo pipefail
  $GATK_ROOT/bin/gatk MarkDuplicates \
                      -I ~{inputBam} \
                      --METRICS_FILE ~{outputFileNamePrefix}.dupmetrics \
                      --VALIDATION_STRINGENCY SILENT \
                      --CREATE_INDEX true \
                      -O ~{outputFileNamePrefix}.dupmarked.bam
>>>

 runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
 }

 output {
  File bam = "~{outputFileNamePrefix}.dupmarked.bam"
  File bamIndex = "~{outputFileNamePrefix}.dupmarked.bai"
 }
}

# ==========================================
#  coverage metrics from the bam file used for fingerprint analysis
# ==========================================


 task alignmentMetrics{
   input{
    File inputBam
    File inputBai
    String modules
    String outputFileNamePrefix
    Boolean markDups	
    Int jobMemory = 8
    Int timeout = 24
   }
   parameter_meta {
    inputBam: "input .bam file"
    inputBai: "index of the input .bam file"
    outputFileNamePrefix: "prefix for making names for output files"  
    jobMemory: "memory allocated for Job"
    modules: "Names and versions of modules"
    timeout: "Timeout in hours, needed to override imposed limits"
   } 

command <<<
  set -euo pipefail

  ### samtools stats
  $SAMTOOLS_ROOT/bin/samtools stats ~{inputBam} > ~{outputFileNamePrefix}.samstats.txt
  reads=`cat ~{outputFileNamePrefix}.samstats.txt | grep ^SN | grep "raw total sequences:" | cut -f3`
  mapped_reads=`cat ~{outputFileNamePrefix}.samstats.txt | grep ^SN | grep "reads mapped:" | cut -f 3`
  unmapped_reads=`cat ~{outputFileNamePrefix}.samstats.txt | grep ^SN | grep "reads unmapped:" | cut -f 3`
  mapped_bases=`cat ~{outputFileNamePrefix}.samstats.txt | grep ^SN | grep "bases mapped:" | cut -f 3`
  reads_duplicated=`cat ~{outputFileNamePrefix}.samstats.txt | grep ^SN | grep "reads duplicated:" | cut -f 3`

  ### samtools coverage, with duplicates
  $SAMTOOLS_ROOT/bin/samtools coverage --ff UNMAP,SECONDARY,QCFAIL ~{inputBam} > ~{outputFileNamePrefix}.coverage.txt
  mean_cvg=`cat ~{outputFileNamePrefix}.coverage.txt | grep -P "^chr\d+\t|^chrX\t|^chrY\t" | awk '{ space += ($3-$2)+1; bases += $7*($3-$2);} END { print bases/space }'`
  
  ### samtools coverage, deduplicated
  $SAMTOOLS_ROOT/bin/samtools coverage --ff UNMAP,SECONDARY,QCFAIL,DUP ~{inputBam} > ~{outputFileNamePrefix}.dedup.coverage.txt
  mean_dedup_cvg=`cat ~{outputFileNamePrefix}.dedup.coverage.txt | grep -P "^chr\d+\t|^chrX\t|^chrY\t" | awk '{ space += ($3-$2)+1; bases += $7*($3-$2);} END { print bases/space }'`
  
  ### json file
  echo \{\"reads\":$reads,\"mapped_reads\":$mapped_reads,\"unmapped_reads\":$unmapped_reads,\"mapped_bases\":$mapped_bases,\"reads_duplicated\":$reads_duplicated,\"mean_raw_cvg\":$mean_cvg,\"mean_dedup_cvg\":$mean_dedup_cvg\,\"markDups\":~{markDups}} > ~{outputFileNamePrefix}.json



>>>

  runtime {
   memory:  "~{jobMemory} GB"
   modules: "~{modules}"
   timeout: "~{timeout}"
  }

  output {
    File json = "~{outputFileNamePrefix}.json"
    File samstats = "~{outputFileNamePrefix}.samstats.txt"
  }
}




