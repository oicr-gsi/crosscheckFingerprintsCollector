version 1.0

import "imports/pull_bwaMem.wdl" as bwaMem
import "imports/pull_star.wdl" as star

struct InputGroup {
  File fastqR1
  File fastqR2
  String readGroup
}

struct GenomeResources {
    String refFasta
    String bwaRef
    String refHapMap
    String bwaMemModules
    String starRefDir
    String starModules
    String extractFingerprintModules
    String intervalsToParallelizeByString
    String alignmentMetricsModules
    String markDuplicatesModules
    String downsampleModules
    String intervalBed
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
        Boolean filterBam
        String outputFileNamePrefix
        String reference
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
        filterBam: "should use filterBamToInterval to prefiltering of the bam file to intervals? Generally true"
        intervalsToParallelizeByString : "Comma separated list of intervals to split by (e.g. chr1,chr2,chr3+chr4)."
        outputFileNamePrefix: "Optional output prefix for the output"
        reference : "the reference genome for input sample"
        maxReads: "The maximum number of reads to process; if set, this will sample the requested number of reads"
        sampleId : "value that will be used as the sample identifier in the vcf fingerprint"
   }

Map[String,GenomeResources] resources = {
  "hg38": {
    "refFasta" : "$HG38_ROOT/hg38_random.fa",
    "bwaRef" : "$HG38_BWA_INDEX_ROOT/hg38_random.fa",
    "refHapMap" : "$CROSSCHECKFINGERPRINTS_HAPLOTYPE_MAP_ROOT/oicr_hg38_chr.map",
    "intervalBed": "$CROSSCHECKFINGERPRINTS_HAPLOTYPE_MAP_ROOT/oicr_hg38_intervals.bed",
    "bwaMemModules" : "samtools/1.9 bwa/0.7.17 hg38-bwa-index-with-alt/0.7.17",
    "starRefDir" : "$HG38_STAR_INDEX100_ROOT",
    "starModules" :"star/2.7.3a hg38-star-index100/2.7.3a",
    "extractFingerprintModules" : "gatk/4.2.0.0 tabix/0.2.6 hg38/p12 crosscheckfingerprints-haplotype-map/20230324",
    "alignmentMetricsModules" : "samtools/1.15",
    "markDuplicatesModules" : "gatk/4.2.0.0 samtools/1.15",
    "downsampleModules" :  "seqtk/1.3",
    "intervalsToParallelizeByString" : "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM"
  },
  "hg19": {
    "refFasta" : "$HG19_ROOT/hg19_random.fa",
    "bwaRef" : "$HG19_BWA_INDEX_ROOT/hg19_random.fa",
    "refHapMap" : "$CROSSCHECKFINGERPRINTS_HAPLOTYPE_MAP_ROOT/oicr_hg19_chr.map",
    "intervalBed": "$CROSSCHECKFINGERPRINTS_HAPLOTYPE_MAP_ROOT/oicr_hg19_intervals.bed",
    "bwaMemModules" : "samtools/1.9 bwa/0.7.17 hg19-bwa-index/0.7.17",
    "starRefDir" : "$HG19_STAR_INDEX100_ROOT",
    "starModules" : "star/2.7.3a  hg19-star-index100/2.7.3a",
    "extractFingerprintModules" : "gatk/4.2.0.0 tabix/0.2.6 hg19/p13 crosscheckfingerprints-haplotype-map/20230324",
    "alignmentMetricsModules" : "samtools/1.15",
    "markDuplicatesModules" : "gatk/4.2.0.0 samtools/1.15",
    "downsampleModules" :  "seqtk/1.3",
    "intervalsToParallelizeByString" : "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM"
  }}

   if(inputType=="fastq" && defined(fastqR1) && defined(fastqR2)){
     if(maxReads>0){
      call downsample{
        input:
          fastqR1 = select_first([fastqR1]),
          fastqR2 = select_first([fastqR2]),
          maxReads = maxReads,
          modules = resources [ reference ].downsampleModules
      }
     }

     if(aligner=="bwa"){
       call bwaMem.bwaMem {
         input:
           fastqR1 = select_first([downsample.fastqR1mod,fastqR1]),
           fastqR2 = select_first([downsample.fastqR2mod,fastqR2]),
           outputFileNamePrefix = outputFileNamePrefix,
           readGroups = "'@RG\\tID:CROSSCHECK\\tSM:SAMPLE'",
           doTrim = false,
           runBwaMem_bwaRef = resources [ reference ].bwaRef,
           runBwaMem_modules = resources [ reference ].bwaMemModules
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
           runStar_chimOutType = "Junctions",
           runStar_genomeIndexDir = resources [ reference].starRefDir,
           runStar_modules = resources [ reference].starModules
       }
     }
   }

  if(filterBam){
    call filterBamToIntervals {
      input:
         inputBam = select_first([bwaMem.bwaMemBam,star.starBam,bam]),
         inputBai = select_first([bwaMem.bwaMemIndex,star.starIndex,bamIndex]),
         intervalBed = resources [ reference].intervalBed,
         outputFileNamePrefix = outputFileNamePrefix,
         modules = "crosscheckfingerprints-haplotype-map/20230324 samtools/1.14"     
    }
  }


  if(markDups){

    call splitStringToArray {
      input:
        str = resources [ reference].intervalsToParallelizeByString
    }

    Array[Array[String]] intervalsToParallelizeBy = splitStringToArray.out

    scatter (intervals in intervalsToParallelizeBy){
      call markDuplicates {
        input :
          inputBam = select_first([filterBamToIntervals.bam,bwaMem.bwaMemBam,star.starBam,bam]),
          inputBai = select_first([filterBamToIntervals.bamIndex,bwaMem.bwaMemIndex,star.starIndex,bamIndex]),
          outputFileNamePrefix = outputFileNamePrefix,
          intervals = intervals,
          modules = resources [ reference ].markDuplicatesModules 
      }
    }
    Array[File] markDuplicatedBams = markDuplicates.bam
    
	
	call mergeBams {
      input:
      bams = markDuplicatedBams,
      outputFileName = outputFileNamePrefix,
      suffix = ""
    }
  }

  
   call alignmentMetrics {
     input:
        inputBam = select_first([mergeBams.mergedBam,filterBamToIntervals.bam,bwaMem.bwaMemBam,star.starBam,bam]),
        inputBai = select_first([mergeBams.mergedBamIndex,filterBamToIntervals.bamIndex,bwaMem.bwaMemIndex,star.starIndex,bamIndex]),
        outputFileNamePrefix = outputFileNamePrefix,
        markDups = markDups,
        maxReads = maxReads,
        modules = resources [ reference ].alignmentMetricsModules  
   }
   
   call extractFingerprint {
     input:
        inputBam = select_first([mergeBams.mergedBam,filterBamToIntervals.bam,bwaMem.bwaMemBam,star.starBam,bam]),
        inputBai = select_first([mergeBams.mergedBamIndex,filterBamToIntervals.bamIndex,star.starIndex,bamIndex]),
        haplotypeMap = resources [ reference ].refHapMap,
        refFasta =  resources [ reference ].refFasta,
        outputFileNamePrefix = outputFileNamePrefix,
        sampleId = sampleId,
        modules = resources [ reference ].extractFingerprintModules
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
      },
      {
        name: "samtools/1.1",
        url: "http://www.htslib.org/"
      },
      {
        name: "seqtk/1.3",
        url: "https://github.com/lh3/seqtk"
      },
      { name: "gsi crosscheckfingerprints-haplotype-map module",
        url: "https://gitlab.oicr.on.ca/ResearchIT/modulator"
      },
      { 
        name: "gsi hg38 modules : hg38/p12 hg38-bwa-index-with-alt/0.7.17 hg38-star-index100/2.7.3a",
        url: "https://gitlab.oicr.on.ca/ResearchIT/modulator"
      },
      {
        name: "gsi hg19 modules : hg19/p13 hg19-bwa-index/0.7.17 hg19-star-index100/2.7.3a",
        url: "https://gitlab.oicr.on.ca/ResearchIT/modulator"
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
#  Filter Bam to Intervals
# ==========================================



task filterBamToIntervals {
 input{
  File inputBam
  File inputBai
  String intervalBed
  String modules
  String outputFileNamePrefix
  Int jobMemory = 16
  Int overhead = 6
  Int timeout = 24
 }
 parameter_meta {
  inputBam: "input .bam file"
  inputBai: "index of the input .bam file"
  outputFileNamePrefix: "prefix for making names for output files"  
  jobMemory: "memory allocated for Job"
  overhead: "memory allocated to overhead of the job other than used in markDuplicates command"
  modules: "Names and versions of modules"
  timeout: "Timeout in hours, needed to override imposed limits"
 }
 
command <<<
  set -euo pipefail
  samtools view -b ~{inputBam} -L ~{intervalBed} > ~{outputFileNamePrefix}.filtered.bam
  samtools index ~{outputFileNamePrefix}.filtered.bam
>>>

 runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
 }

 output {
  File bam = "~{outputFileNamePrefix}.filtered.bam"
  File bamIndex = "~{outputFileNamePrefix}.filtered.bam.bai"
 }
}





# ==========================================
#  Split a string to array
# ==========================================
task splitStringToArray {
  input {
    String str
    String lineSeparator = ","
    String recordSeparator = "+"

    Int jobMemory = 1
    Int cores = 1
    Int timeout = 1
    String modules = ""
  }

  command <<<
    set -euo pipefail

    echo "~{str}" | tr '~{lineSeparator}' '\n' | tr '~{recordSeparator}' '\t'
  >>>

  output {
    Array[Array[String]] out = read_tsv(stdout())
  }

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules}"
  }

  parameter_meta {
    str: "Interval string to split (e.g. chr1,chr2,chr3+chr4)."
    lineSeparator: "Interval group separator - these are the intervals to split by."
    recordSeparator: "Interval interval group separator - this can be used to combine multiple intervals into one group."
    jobMemory: "Memory allocated to job (in GB)."
    cores: "The number of cores to allocate to the job."
    timeout: "Maximum amount of time (in hours) the task can run for."
    modules: "Environment module name and version to load (space separated) before command execution."
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
                    --SAMPLE_ALIAS ~{sampleId}

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
  Array[String] intervals
  Int jobMemory = 16
  Int overhead = 6
  Int timeout = 24
 }
 parameter_meta {
  inputBam: "input .bam file"
  inputBai: "index of the input .bam file"
  outputFileNamePrefix: "prefix for making names for output files"  
  jobMemory: "memory allocated for Job"
  overhead: "memory allocated to overhead of the job other than used in markDuplicates command"
  modules: "Names and versions of modules"
  timeout: "Timeout in hours, needed to override imposed limits"
 }
 
command <<<
  set -euo pipefail
  ln -s ~{inputBam} inputBam.bam
  ln -s ~{inputBai} inputBam.bam.bai
  samtools view -b \
        inputBam.bam \
        ~{sep=" " intervals} > intervalBam.bam
  samtools index intervalBam.bam intervalBam.bam.bai

  $GATK_ROOT/bin/gatk --java-options "-Xmx~{jobMemory - overhead}G" MarkDuplicates \
                      -I intervalBam.bam \
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
 }
}

# ==========================================
#  Merge Bam files
# ==========================================

task mergeBams {
  input {
    Array[File] bams
    String outputFileName
    String suffix = ".merge"
    String? additionalParams

    Int jobMemory = 24
    Int overhead = 6
    Int cores = 1
    Int timeout = 6
    String modules = "gatk/4.1.6.0"
  }

  command <<<
    set -euo pipefail

    gatk --java-options "-Xmx~{jobMemory - overhead}G" MergeSamFiles \
    ~{sep=" " prefix("--INPUT=", bams)} \
    --OUTPUT="~{outputFileName}~{suffix}.bam" \
    --CREATE_INDEX=true \
    --SORT_ORDER=coordinate \
    --ASSUME_SORTED=false \
    --USE_THREADING=true \
    --VALIDATION_STRINGENCY=SILENT \
    ~{additionalParams}
  >>>

  output {
    File mergedBam = "~{outputFileName}~{suffix}.bam"
    File mergedBamIndex = "~{outputFileName}~{suffix}.bai"
  }

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules}"
  }

  parameter_meta {
    bams: "Array of bam files to merge together."
    outputFileName: "Output files will be prefixed with this."
    additionalParams: "Additional parameters to pass to GATK MergeSamFiles."
    jobMemory: "Memory allocated to job (in GB)."
    overhead: "Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory."
    cores: "The number of cores to allocate to the job."
    timeout: "Maximum amount of time (in hours) the task can run for."
    modules: "Environment module name and version to load (space separated) before command execution."
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
    String maxReads
   }
   parameter_meta {
    inputBam: "input .bam file"
    inputBai: "index of the input .bam file"
    outputFileNamePrefix: "prefix for making names for output files"  
    jobMemory: "memory allocated for Job"
    modules: "Names and versions of modules"
    maxReads : "the maximum number of reads to use"
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
  echo \{\"reads\":$reads,\"mapped_reads\":$mapped_reads,\"unmapped_reads\":$unmapped_reads,\"mapped_bases\":$mapped_bases,\"reads_duplicated\":$reads_duplicated,\"mean_raw_cvg\":$mean_cvg,\"mean_dedup_cvg\":$mean_dedup_cvg\,\"markDups\":~{markDups},\"maxReads\":~{maxReads}} > ~{outputFileNamePrefix}.json



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