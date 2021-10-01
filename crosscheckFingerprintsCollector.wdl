version 1.0

import "imports/pull_bwaMem.wdl" as bwaMem

workflow crosscheckFingerprintsCollector {
   input {
        File fastqR1
        File fastqR2
        String outputFileNamePrefix = basename(fastqR1)
        String refFasta
        String haplotypeMap
   }

   call bwaMem.bwaMem {
        input:
          fastqR1 = fastqR1,
          fastqR2 = fastqR2,
          outputFileNamePrefix = outputFileNamePrefix,
          readGroups = "'@RG\\tID:ID\\tSM:SAMPLE'",
	  doTrim = false
   }

   call extractFingerprint {
        input:
 	    inputBam = bwaMem.bwaMemBam,
            inputBai = bwaMem.bwaMemIndex,
            haplotypeMap = haplotypeMap,
            refFasta = refFasta,
            outputFileNamePrefix = outputFileNamePrefix
   }

   parameter_meta {
        fastqR1: "fastq file for read 1"
        fastqR2: "fastq file for read 2"
        outputFileNamePrefix: "Optional output prefix for the output"
	refFasta: "Path to the reference fasta file"
	haplotypeMap: "Path to the gzipped hotspot vcf file"
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
      },
    ]
    output_meta: {
      outputVcf: "gzipped vcf expression levels for all genes recorded in the reference",
      outbutTbi: "expression levels for all isoforms recorded in the reference",
    }
  }

  output {
    File outputVcf = extractFingerprint.vgz
    File outputTbi = extractFingerprint.tbi
   }
}

# ==========================================
#  configure and run extractFingerprint
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


