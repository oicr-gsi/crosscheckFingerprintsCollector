# crosscheckFingerprintsCollector

crosscheckFingerprintsCollector, workflow that generates genotype fingerprints using gatk ExtractFingprint.  Output are vcf files that can be proccessed through gatk Crosscheck fingerprints
##

## Overview

## Dependencies

* [gatk 4.2.0.0](https://gatk.broadinstitute.org)
* [tabix 0.2.6](http://www.htslib.org)
* [samtools 1.1](http://www.htslib.org/)
* [seqtk 1.3](https://github.com/lh3/seqtk)
* [gsi crosscheckfingerprints-haplotype-map module](https://gitlab.oicr.on.ca/ResearchIT/modulator)
* [gsi hg38 modules : hg38 p12 hg38-bwa-index-with-alt 0.7.17 hg38-star-index100 2.7.3a](https://gitlab.oicr.on.ca/ResearchIT/modulator)
* [gsi hg19 modules : hg19 p13 hg19-bwa-index 0.7.17 hg19-star-index100 2.7.3a](https://gitlab.oicr.on.ca/ResearchIT/modulator)


## Usage

### Cromwell
```
java -jar cromwell.jar run crosscheckFingerprintsCollector.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`inputType`|String|one of fastq, bam, or cram
`aligner`|String|aligner to use for fastq input, either bwa or star
`markDups`|Boolean|should the alignment be duplicate marked?, generally yes
`filterBam`|Boolean|should use filterBamToInterval to prefiltering of the bam/cram file to intervals? Generally true
`outputFileNamePrefix`|String|Optional output prefix for the output
`reference`|String|the reference genome for input sample
`sampleId`|String|value that will be used as the sample identifier in the vcf fingerprint


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`fastqR1`|File?|None|fastq file for read 1
`fastqR2`|File?|None|fastq file for read 2
`bam`|File?|None|bam file
`bamIndex`|File?|None|bam index file
`cram`|File?|None|cram file (may be a merged-lanes cram)
`cramIndex`|File?|None|index for the cram file
`is_lane_level`|Boolean|true|true if the input bam/cram is already at lane level; false if it is a merged-lanes file that needs to be split before processing
`maxReads`|Int|0|The maximum number of reads to process; if set, this will sample the requested number of reads


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`downsample.jobMemory`|Int|8|memory allocated for Job
`downsample.timeout`|Int|24|Timeout in hours, needed to override imposed limits
`bwaMem.adapterTrimmingLog_timeout`|Int|48|Hours before task timeout
`bwaMem.adapterTrimmingLog_jobMemory`|Int|12|Memory allocated indexing job
`bwaMem.indexBam_timeout`|Int|48|Hours before task timeout
`bwaMem.indexBam_modules`|String|"samtools/1.9"|Modules for running indexing job
`bwaMem.indexBam_jobMemory`|Int|12|Memory allocated indexing job
`bwaMem.bamMerge_timeout`|Int|72|Hours before task timeout
`bwaMem.bamMerge_modules`|String|"samtools/1.9"|Required environment modules
`bwaMem.bamMerge_jobMemory`|Int|32|Memory allocated indexing job
`bwaMem.runBwaMem_timeout`|Int|96|Hours before task timeout
`bwaMem.runBwaMem_jobMemory`|Int|32|Memory allocated for this job
`bwaMem.runBwaMem_threads`|Int|8|Requested CPU threads
`bwaMem.runBwaMem_addParam`|String?|None|Additional BWA parameters
`bwaMem.adapterTrimming_timeout`|Int|48|Hours before task timeout
`bwaMem.adapterTrimming_jobMemory`|Int|16|Memory allocated for this job
`bwaMem.adapterTrimming_addParam`|String?|None|Additional cutadapt parameters
`bwaMem.adapterTrimming_modules`|String|"cutadapt/1.8.3"|Required environment modules
`bwaMem.slicerR2_timeout`|Int|48|Hours before task timeout
`bwaMem.slicerR2_jobMemory`|Int|16|Memory allocated for this job
`bwaMem.slicerR2_modules`|String|"slicer/0.3.0"|Required environment modules
`bwaMem.slicerR1_timeout`|Int|48|Hours before task timeout
`bwaMem.slicerR1_jobMemory`|Int|16|Memory allocated for this job
`bwaMem.slicerR1_modules`|String|"slicer/0.3.0"|Required environment modules
`bwaMem.countChunkSize_timeout`|Int|48|Hours before task timeout
`bwaMem.countChunkSize_jobMemory`|Int|16|Memory allocated for this job
`bwaMem.numChunk`|Int|1|number of chunks to split fastq file [1, no splitting]
`bwaMem.trimMinLength`|Int|1|minimum length of reads to keep [1]
`bwaMem.trimMinQuality`|Int|0|minimum quality of read ends to keep [0]
`bwaMem.adapter1`|String|"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"|adapter sequence to trim from read 1 [AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC]
`bwaMem.adapter2`|String|"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"|adapter sequence to trim from read 2 [AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT]
`star.indexBam_timeout`|Int|48|hours before task timeout
`star.indexBam_modules`|String|"picard/2.19.2"|modules for running indexing job
`star.indexBam_jobMemory`|Int|12|Memory allocated indexing job
`star.runStar_timeout`|Int|72|hours before task timeout
`star.runStar_jobMemory`|Int|64|Memory allocated for this job
`star.runStar_threads`|Int|6|Requested CPU threads
`star.runStar_peOvMMp`|Float|0.1|maximum proportion of mismatched bases in the overlap area
`star.runStar_chimSegmentReadGapMax`|Int|3|maximum gap in the read sequence between chimeric segments
`star.runStar_peOvNbasesMin`|Int|10|minimum number of overlap bases to trigger mates merging and realignment
`star.runStar_chimOutJunForm`|Int?|None|flag to add metadata to chimeric junction output for functionality with starFusion - 1 for metadata, 0 for no metadata
`star.runStar_chimNonchimScoDMin`|Int|10|to trigger chimeric detection, the drop in the best non-chimeric alignment score with respect to the read length has to be greater than this value
`star.runStar_chimMulmapNmax`|Int|50|maximum number of chimeric multi-alignments
`star.runStar_chimScoreSeparation`|Int|1|minimum difference (separation) between the best chimeric score and the next one
`star.runStar_chimScoJunNonGTAG`|Int|0|penalty for a non-GTAG chimeric junction
`star.runStar_chimMulmapScoRan`|Int|3|the score range for multi-mapping chimeras below the best chimeric score
`star.runStar_alignIntMax`|Int|100000|maximum intron size
`star.runStar_alignMatGapMax`|Int|100000|maximum gap between two mates
`star.runStar_alignSJDBOvMin`|Int|10|minimum overhang for annotated spliced alignments
`star.runStar_chimJunOvMin`|Int|10|minimum overhang for a chimeric junction
`star.runStar_chimSegmin`|Int|10|minimum length of chimeric segment length
`star.runStar_multiMax`|Int|-1|multiMax parameter for STAR
`star.runStar_saSparsed`|Int|2|saSparsed parameter for STAR
`star.runStar_uniqMAPQ`|Int|255|Score for unique mappers
`star.runStar_chimScoreDropMax`|Int|30|max drop (difference) of chimeric score (the sum of scores of allchimeric segments) from the read length
`star.runStar_outFilterMultimapNmax`|Int|50|max number of multiple alignments allowed for a read: if exceeded, the read is considered unmapped
`star.runStar_addParam`|String?|None|Additional STAR parameters
`star.runStar_genereadSuffix`|String|"ReadsPerGene.out"|ReadsPerGene file suffix
`star.runStar_chimericjunctionSuffix`|String|"Chimeric.out"|Suffix for chimeric junction file
`star.runStar_transcriptomeSuffix`|String|"Aligned.toTranscriptome.out"|Suffix for transcriptome-aligned file
`star.runStar_starSuffix`|String|"Aligned.sortedByCoord.out"|Suffix for sorted file
`filterBamPreSplit.jobMemory`|Int|16|memory allocated for Job
`filterBamPreSplit.overhead`|Int|6|memory allocated to overhead of the job other than used in markDuplicates command
`filterBamPreSplit.timeout`|Int|24|Timeout in hours, needed to override imposed limits
`splitLanes.jobMemory`|Int|16|memory allocated for Job
`splitLanes.timeout`|Int|24|Timeout in hours, needed to override imposed limits
`splitStringToArray.lineSeparator`|String|","|Interval group separator - these are the intervals to split by.
`splitStringToArray.recordSeparator`|String|"+"|Interval interval group separator - this can be used to combine multiple intervals into one group.
`splitStringToArray.jobMemory`|Int|1|Memory allocated to job (in GB).
`splitStringToArray.cores`|Int|1|The number of cores to allocate to the job.
`splitStringToArray.timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`splitStringToArray.modules`|String|""|Environment module name and version to load (space separated) before command execution.
`filterBamLane.jobMemory`|Int|16|memory allocated for Job
`filterBamLane.overhead`|Int|6|memory allocated to overhead of the job other than used in markDuplicates command
`filterBamLane.timeout`|Int|24|Timeout in hours, needed to override imposed limits
`markDuplicates.jobMemory`|Int|16|memory allocated for Job
`markDuplicates.overhead`|Int|6|memory allocated to overhead of the job other than used in markDuplicates command
`markDuplicates.timeout`|Int|24|Timeout in hours, needed to override imposed limits
`mergeIntervalBams.additionalParams`|String?|None|Additional parameters to pass to GATK MergeSamFiles.
`mergeIntervalBams.jobMemory`|Int|24|Memory allocated to job (in GB).
`mergeIntervalBams.overhead`|Int|6|Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory.
`mergeIntervalBams.cores`|Int|1|The number of cores to allocate to the job.
`mergeIntervalBams.timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`mergeIntervalBams.modules`|String|"gatk/4.1.6.0"|Environment module name and version to load (space separated) before command execution.
`alignmentMetrics.jobMemory`|Int|8|memory allocated for Job
`alignmentMetrics.timeout`|Int|24|Timeout in hours, needed to override imposed limits
`extractFingerprint.jobMemory`|Int|8|memory allocated for Job
`extractFingerprint.timeout`|Int|24|Timeout in hours, needed to override imposed limits


### Outputs

Output | Type | Description | Labels
---|---|---|---
`outputVcf`|Pair[Array[File]+,Map[String,String]]|per-lane crosscheck fingerprint vcf.gz files, file names carry read group|
`outputTbi`|Pair[Array[File]+,Map[String,String]]|per-lane vcf.gz.tbi index files, file names carry read group|
`json`|Pair[Array[File]+,Map[String,String]]|per-lane alignment metrics json files, file names carry read group|
`samstats`|Pair[Array[File]+,Map[String,String]]|per-lane samstats summary files, file names carry read group|


## Commands
This section lists command(s) run by crosscheckFingerprintsCollector workflow

* Running crosscheckFingerprintsCollector

```
    set -euo pipefail
    EXT=$(basename ~{inputBam} | rev | cut -d. -f1 | rev)
    if [ "$EXT" = "cram" ]; then
      # samtools split does not support -T; convert CRAM to BAM first
      ln -s ~{inputBam} input.cram
      ln -s ~{inputBai} input.cram.crai
      samtools view -b -T ~{refFasta} -o input_converted.bam input.cram
      samtools index input_converted.bam
      samtools split -f "~{outputFileNamePrefix}_%!.bam" input_converted.bam
    else
      ln -s ~{inputBam} input.bam
      ln -s ~{inputBai} input.bam.bai
      samtools split -f "~{outputFileNamePrefix}_%!.bam" input.bam
    fi
    for f in ~{outputFileNamePrefix}_*.bam; do samtools index "$f"; done
```
```
  set -euo pipefail
  EXT=$(basename ~{inputBam} | rev | cut -d. -f1 | rev)
  ln -s ~{inputBam} input.$EXT
  if [ "$EXT" = "cram" ]; then ln -s ~{inputBai} input.cram.crai
  else                          ln -s ~{inputBai} input.bam.bai
  fi
  samtools view -b -T ~{refFasta} -L ~{intervalBed} input.$EXT > ~{outputFileNamePrefix}.filtered.bam
  samtools index ~{outputFileNamePrefix}.filtered.bam
```
```
    set -euo pipefail

    echo "~{str}" | tr '~{lineSeparator}' '\n' | tr '~{recordSeparator}' '\t'
```
```
  set -euo pipefail

 $GATK_ROOT/bin/gatk ExtractFingerprint \
                    -R ~{refFasta} \
                    -H ~{haplotypeMap} \
                    -I ~{inputBam} \
                    -O ~{outputFileNamePrefix}.vcf \
                    --SAMPLE_ALIAS ~{sampleId}

 $TABIX_ROOT/bin/bgzip -c ~{outputFileNamePrefix}.vcf > ~{outputFileNamePrefix}.vcf.gz
 $TABIX_ROOT/bin/tabix -p vcf ~{outputFileNamePrefix}.vcf.gz
```
```
 set -euo pipefail

 seqtk sample -s 100 ~{fastqR1} ~{maxReads} > ~{fastqR1m}
 gzip ~{fastqR1m}

 seqtk sample -s 100 ~{fastqR2} ~{maxReads} > ~{fastqR2m}
 gzip ~{fastqR2m}
```
```
  set -euo pipefail
  EXT=$(basename ~{inputBam} | rev | cut -d. -f1 | rev)
  ln -s ~{inputBam} input.$EXT
  if [ "$EXT" = "cram" ]; then ln -s ~{inputBai} input.cram.crai
  else                          ln -s ~{inputBai} input.bam.bai
  fi
  samtools view -b -T ~{refFasta} input.$EXT \
        ~{sep=" " intervals} > intervalBam.bam
  samtools index intervalBam.bam intervalBam.bam.bai

  $GATK_ROOT/bin/gatk --java-options "-Xmx~{jobMemory - overhead}G" MarkDuplicates \
                      -I intervalBam.bam \
                      --METRICS_FILE ~{outputFileNamePrefix}.dupmetrics \
                      --VALIDATION_STRINGENCY SILENT \
                      --CREATE_INDEX true \
                      -O ~{outputFileNamePrefix}.dupmarked.bam
```
```
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
```
```
  set -euo pipefail

  ### samtools stats
  $SAMTOOLS_ROOT/bin/samtools stats --reference ~{refFasta} ~{inputBam} > ~{outputFileNamePrefix}.samstats.txt
  reads=`cat ~{outputFileNamePrefix}.samstats.txt | grep ^SN | grep "raw total sequences:" | cut -f3`
  mapped_reads=`cat ~{outputFileNamePrefix}.samstats.txt | grep ^SN | grep "reads mapped:" | cut -f 3`
  unmapped_reads=`cat ~{outputFileNamePrefix}.samstats.txt | grep ^SN | grep "reads unmapped:" | cut -f 3`
  mapped_bases=`cat ~{outputFileNamePrefix}.samstats.txt | grep ^SN | grep "bases mapped:" | cut -f 3`
  reads_duplicated=`cat ~{outputFileNamePrefix}.samstats.txt | grep ^SN | grep "reads duplicated:" | cut -f 3`

  ### samtools coverage, with duplicates
  $SAMTOOLS_ROOT/bin/samtools coverage --ff UNMAP,SECONDARY,QCFAIL --reference ~{refFasta} ~{inputBam} > ~{outputFileNamePrefix}.coverage.txt
  mean_cvg=`cat ~{outputFileNamePrefix}.coverage.txt | grep -P "^chr\d+\t|^chrX\t|^chrY\t" | awk '{ space += ($3-$2)+1; bases += $7*($3-$2);} END { print bases/space }'`

  ### samtools coverage, deduplicated
  $SAMTOOLS_ROOT/bin/samtools coverage --ff UNMAP,SECONDARY,QCFAIL,DUP --reference ~{refFasta} ~{inputBam} > ~{outputFileNamePrefix}.dedup.coverage.txt
  mean_dedup_cvg=`cat ~{outputFileNamePrefix}.dedup.coverage.txt | grep -P "^chr\d+\t|^chrX\t|^chrY\t" | awk '{ space += ($3-$2)+1; bases += $7*($3-$2);} END { print bases/space }'`

  ### json file
  echo \{\"reads\":$reads,\"mapped_reads\":$mapped_reads,\"unmapped_reads\":$unmapped_reads,\"mapped_bases\":$mapped_bases,\"reads_duplicated\":$reads_duplicated,\"mean_raw_cvg\":$mean_cvg,\"mean_dedup_cvg\":$mean_dedup_cvg\,\"markDups\":~{markDups},\"maxReads\":~{maxReads}} > ~{outputFileNamePrefix}.json


```
## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
