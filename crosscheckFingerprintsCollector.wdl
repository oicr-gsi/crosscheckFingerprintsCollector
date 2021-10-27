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
        String outputFileNamePrefix
        String refFasta
        String haplotypeMap
   }
   parameter_meta {
        fastqR1: "fastq file for read 1"
        fastqR2: "fastq file for read 2"
        bam: "bam file"
        bamIndex: "bam index file"
	    inputType: "one of either fastq or bam"
		aligner : "aligner to use for fastq input, either bwa or star"
        outputFileNamePrefix: "Optional output prefix for the output"
        refFasta: "Path to the reference fasta file"
        haplotypeMap: "Path to the gzipped hotspot vcf file"
   }


   if(inputType=="fastq" && defined(fastqR1) && defined(fastqR2)){
     
	 if(aligner=="bwa"){
       call bwaMem.bwaMem {
         input:
           fastqR1 = select_first([fastqR1]),
           fastqR2 = select_first([fastqR2]),
           outputFileNamePrefix = outputFileNamePrefix,
           readGroups = "'@RG\\tID:ID\\tSM:SAMPLE'",
           doTrim = false
        }
	 }

	 if(aligner=="star"){
       InputGroup starInput = { 
         "fastqR1": select_first([fastqR1]),
         "fastqR2": select_first([fastqR2]),
         "readGroup": "'ID:ID\\tSM:SAMPLE'"
	   }
       call star.star { 
		  input:
            inputGroups = [ starInput ],
            outputFileNamePrefix = outputFileNamePrefix
		}


	  }
   }
  


   call extractFingerprint {
     input:
        inputBam = select_first([bwaMem.bwaMemBam,star.starBam,bam]),
        inputBai = select_first([bwaMem.bwaMemIndex,star.starIndex,bamIndex]),
        haplotypeMap = haplotypeMap,
        refFasta = refFasta,
        outputFileNamePrefix = outputFileNamePrefix
    }

    output {
      File outputVcf = extractFingerprint.vgz
      File outputTbi = extractFingerprint.tbi
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
       outputVcf: "gzipped vcf expression levels for all genes recorded in the reference",
       outputTbi: "expression levels for all isoforms recorded in the reference"
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
 Int jobMemory = 8
 Int timeout = 24
}
parameter_meta {
 inputBam: "input .bam file"
 inputBai: "index of the input .bam file"
 refFasta: "Path to reference FASTA file"
 outputFileNamePrefix: "prefix for making names for output files"
 haplotypeMap: "Hotspot SNPs are the locations of variants used for genotyping"
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
                    -O ~{outputFileNamePrefix}.vcf

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


