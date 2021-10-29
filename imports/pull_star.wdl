version 1.0

struct InputGroup {
  File fastqR1
  File fastqR2
  String readGroup
}

workflow star {
  input {
    Int indexBam_timeout = 48
    String indexBam_modules = "picard/2.19.2"
    Int indexBam_jobMemory = 12
    Int runStar_timeout = 72
    Int runStar_jobMemory = 64
    Int runStar_threads = 6
    Float runStar_peOvMMp = 0.1
    Int runStar_chimSegmentReadGapMax = 3
    Int runStar_peOvNbasesMin = 10
    Int? runStar_chimOutJunForm
    Int runStar_chimNonchimScoDMin = 10
    Int runStar_chimMulmapNmax = 50
    Int runStar_chimScoreSeparation = 1
    Int runStar_chimScoJunNonGTAG = 0
    Int runStar_chimMulmapScoRan = 3
    Int runStar_alignIntMax = 100000
    Int runStar_alignMatGapMax = 100000
    Int runStar_alignSJDBOvMin = 10
    Int runStar_chimJunOvMin = 10
    Int runStar_chimSegmin = 10
    Int runStar_multiMax = -1
    Int runStar_saSparsed = 2
    Int runStar_uniqMAPQ = 255
    Int runStar_chimScoreDropMax = 30
    Int runStar_outFilterMultimapNmax = 50
    String runStar_chimOutType = "WithinBAM HardClip Junctions"
    String runStar_modules = "star/2.7.6a hg38-star-index100/2.7.6a"
    String? runStar_addParam
    String runStar_genereadSuffix = "ReadsPerGene.out"
    String runStar_chimericjunctionSuffix = "Chimeric.out"
    String runStar_transcriptomeSuffix = "Aligned.toTranscriptome.out"
    String runStar_starSuffix = "Aligned.sortedByCoord.out"
    String runStar_genomeIndexDir = "$HG38_STAR_INDEX100_ROOT/"
    Array[InputGroup] inputGroups
    String outputFileNamePrefix
  }

  scatter (ig in inputGroups) {
    File read1s       = ig.fastqR1
    File read2s       = ig.fastqR2
    String readGroups = ig.readGroup
  }

  parameter_meta {
      indexBam_timeout: "hours before task timeout"
      indexBam_modules: "modules for running indexing job"
      indexBam_jobMemory: "Memory allocated indexing job"
      runStar_timeout: "hours before task timeout"
      runStar_jobMemory: "Memory allocated for this job"
      runStar_threads: "Requested CPU threads"
      runStar_peOvMMp: "maximum proportion of mismatched bases in the overlap area"
      runStar_chimSegmentReadGapMax: "maximum gap in the read sequence between chimeric segments"
      runStar_peOvNbasesMin: "minimum number of overlap bases to trigger mates merging and realignment"
      runStar_chimOutJunForm: "flag to add metadata to chimeric junction output for functionality with starFusion - 1 for metadata, 0 for no metadata"
      runStar_chimNonchimScoDMin: "to trigger chimeric detection, the drop in the best non-chimeric alignment score with respect to the read length has to be greater than this value"
      runStar_chimMulmapNmax: "maximum number of chimeric multi-alignments"
      runStar_chimScoreSeparation: "minimum difference (separation) between the best chimeric score and the next one"
      runStar_chimScoJunNonGTAG: "penalty for a non-GTAG chimeric junction"
      runStar_chimMulmapScoRan: "the score range for multi-mapping chimeras below the best chimeric score"
      runStar_alignIntMax: "maximum intron size"
      runStar_alignMatGapMax: "maximum gap between two mates"
      runStar_alignSJDBOvMin: "minimum overhang for annotated spliced alignments"
      runStar_chimJunOvMin: "minimum overhang for a chimeric junction"
      runStar_chimSegmin: "minimum length of chimeric segment length"
      runStar_multiMax: "multiMax parameter for STAR"
      runStar_saSparsed: "saSparsed parameter for STAR"
      runStar_uniqMAPQ: "Score for unique mappers"
      runStar_chimScoreDropMax: "max drop (difference) of chimeric score (the sum of scores of allchimeric segments) from the read length"
      runStar_outFilterMultimapNmax: "max number of multiple alignments allowed for a read: if exceeded, the read is considered unmapped"
      runStar_chimOutType: "Indicate where chimeric reads are to be written"
      runStar_modules: "modules for running STAR"
      runStar_addParam: "Additional STAR parameters"
      runStar_genereadSuffix: "ReadsPerGene file suffix"
      runStar_chimericjunctionSuffix: "Suffix for chimeric junction file"
      runStar_transcriptomeSuffix: "Suffix for transcriptome-aligned file"
      runStar_starSuffix: "Suffix for sorted file"
      runStar_genomeIndexDir: "Path to STAR index"
    inputGroups: "Array of fastq files to align with STAR and the merged filename"
    outputFileNamePrefix: "Prefix for filename"
  }

  call runStar {
    input:
    timeout = runStar_timeout,
    jobMemory = runStar_jobMemory,
    threads = runStar_threads,
    peOvMMp = runStar_peOvMMp,
    chimSegmentReadGapMax = runStar_chimSegmentReadGapMax,
    peOvNbasesMin = runStar_peOvNbasesMin,
    chimOutJunForm = runStar_chimOutJunForm,
    chimNonchimScoDMin = runStar_chimNonchimScoDMin,
    chimMulmapNmax = runStar_chimMulmapNmax,
    chimScoreSeparation = runStar_chimScoreSeparation,
    chimScoJunNonGTAG = runStar_chimScoJunNonGTAG,
    chimMulmapScoRan = runStar_chimMulmapScoRan,
    alignIntMax = runStar_alignIntMax,
    alignMatGapMax = runStar_alignMatGapMax,
    alignSJDBOvMin = runStar_alignSJDBOvMin,
    chimJunOvMin = runStar_chimJunOvMin,
    chimSegmin = runStar_chimSegmin,
    multiMax = runStar_multiMax,
    saSparsed = runStar_saSparsed,
    uniqMAPQ = runStar_uniqMAPQ,
    chimScoreDropMax = runStar_chimScoreDropMax,
    outFilterMultimapNmax = runStar_outFilterMultimapNmax,
    chimOutType = runStar_chimOutType,
    modules = runStar_modules,
    addParam = runStar_addParam,
    genereadSuffix = runStar_genereadSuffix,
    chimericjunctionSuffix = runStar_chimericjunctionSuffix,
    transcriptomeSuffix = runStar_transcriptomeSuffix,
    starSuffix = runStar_starSuffix,
    genomeIndexDir = runStar_genomeIndexDir,
    read1s = read1s,
    read2s = read2s,
    readGroups = readGroups,
    outputFileNamePrefix = outputFileNamePrefix
  }

  call indexBam {
   input:
   timeout = indexBam_timeout,
   modules = indexBam_modules,
   jobMemory = indexBam_jobMemory,
   inputBam = runStar.outputBam }

  meta {
   author: "Peter Ruzanov, Alexander Fortuna"
   email: "peter.ruzanov@oicr.on.ca, alexander.fortuna@oicr.on.ca"
   description: "STAR 2.1"
   dependencies: [
      {
        name: "star/2.7.6a",
        url: "https://github.com/alexdobin/STAR"
      },
      {
        name: "picard/2.19.2",
        url: "https://broadinstitute.github.io/picard/"
      }
    ]
}

output {
  File starBam          = runStar.outputBam
  File starIndex        = indexBam.outputBai
  File transcriptomeBam = runStar.transcriptomeBam
  File starChimeric     = runStar.outputChimeric
  File geneReadFile     = runStar.geneReads
 }
}

# ==========================================
#  TASK 1 of 2: run STAR aligner
# ==========================================
task runStar {
input {
  Array[File]+ read1s
  Array[File]+ read2s
  Array[String]+ readGroups
  String genomeIndexDir = "$HG38_STAR_INDEX100_ROOT/"
  String outputFileNamePrefix
  String starSuffix = "Aligned.sortedByCoord.out"
  String transcriptomeSuffix = "Aligned.toTranscriptome.out"
  String chimericjunctionSuffix = "Chimeric.out"
  String genereadSuffix = "ReadsPerGene.out"
  String? addParam
  String modules = "star/2.7.6a hg38-star-index100/2.7.6a"
  String chimOutType = "WithinBAM HardClip Junctions"
  Int outFilterMultimapNmax = 50
  Int chimScoreDropMax = 30
  Int uniqMAPQ = 255
  Int saSparsed = 2
  Int multiMax = -1
  Int chimSegmin = 10
  Int chimJunOvMin = 10
  Int alignSJDBOvMin = 10
  Int alignMatGapMax = 100000
  Int alignIntMax = 100000
  Int chimMulmapScoRan = 3
  Int chimScoJunNonGTAG = 0
  Int chimScoreSeparation = 1
  Int chimMulmapNmax = 50
  Int chimNonchimScoDMin = 10
  Int? chimOutJunForm
  Int peOvNbasesMin = 10
  Int chimSegmentReadGapMax = 3
  Float peOvMMp = 0.1
  Int threads = 6
  Int jobMemory = 64
  Int timeout = 72
}

parameter_meta {
 read1s: "array of read1s"
 read2s: "array of read2s"
 readGroups: "array of readgroup lines"
 starSuffix: "Suffix for sorted file"
 genomeIndexDir: "Path to STAR index"
 transcriptomeSuffix: "Suffix for transcriptome-aligned file"
 chimericjunctionSuffix: "Suffix for chimeric junction file"
 genereadSuffix: "ReadsPerGene file suffix"
 addParam: "Additional STAR parameters"
 modules: "modules for running STAR"
 uniqMAPQ: "Score for unique mappers"
 saSparsed: "saSparsed parameter for STAR"
 multiMax: "multiMax parameter for STAR"
 chimSegmin: "minimum length of chimeric segment length"
 chimJunOvMin: "minimum overhang for a chimeric junction"
 chimOutType: "Indicate where chimeric reads are to be written"
 alignSJDBOvMin: "minimum overhang for annotated spliced alignments"
 alignMatGapMax: "maximum gap between two mates"
 alignIntMax: "maximum intron size"
 chimMulmapScoRan: "the score range for multi-mapping chimeras below the best chimeric score"
 chimScoJunNonGTAG: "penalty for a non-GTAG chimeric junction"
 outFilterMultimapNmax: "max number of multiple alignments allowed for a read: if exceeded, the read is considered unmapped"
 chimMulmapNmax: "maximum number of chimeric multi-alignments"
 chimNonchimScoDMin: "to trigger chimeric detection, the drop in the best non-chimeric alignment score with respect to the read length has to be greater than this value"
 chimOutJunForm: "flag to add metadata to chimeric junction output for functionality with starFusion - 1 for metadata, 0 for no metadata"
 peOvNbasesMin: "minimum number of overlap bases to trigger mates merging and realignment"
 peOvMMp: "maximum proportion of mismatched bases in the overlap area"
 chimScoreDropMax: "max drop (difference) of chimeric score (the sum of scores of allchimeric segments) from the read length"
 chimScoreSeparation: "minimum difference (separation) between the best chimeric score and the next one"
 chimSegmentReadGapMax: "maximum gap in the read sequence between chimeric segments"
 threads: "Requested CPU threads"
 jobMemory: "Memory allocated for this job"
 timeout: "hours before task timeout"
}

# missing --clip3pAdapterSeq $adaptors
command <<<
 set -euo pipefail

 STAR --twopassMode Basic \
      --genomeDir ~{genomeIndexDir} \
      --readFilesIn ~{sep="," read1s} ~{sep="," read2s} \
      --readFilesCommand zcat \
      --outFilterIntronMotifs RemoveNoncanonical \
      --outFileNamePrefix ~{outputFileNamePrefix}. \
      --outSAMmultNmax ~{multiMax} \
      --outSAMattrRGline ~{sep=" , " readGroups} \
      --outSAMstrandField intronMotif \
      --outSAMmapqUnique  ~{uniqMAPQ} \
      --outSAMunmapped Within KeepPairs \
      --genomeSAsparseD ~{saSparsed} \
      --outSAMtype BAM SortedByCoordinate \
      --quantMode TranscriptomeSAM GeneCounts \
      --chimSegmentMin ~{chimSegmin} \
      --chimJunctionOverhangMin ~{chimJunOvMin} \
      --alignSJDBoverhangMin ~{alignSJDBOvMin} \
      --alignMatesGapMax ~{alignMatGapMax} \
      --alignIntronMax ~{alignIntMax} \
      --alignSJstitchMismatchNmax 5 -1 5 5 \
      --chimMultimapScoreRange ~{chimMulmapScoRan} \
      --chimScoreJunctionNonGTAG ~{chimScoJunNonGTAG} \
      --chimMultimapNmax ~{chimMulmapNmax} \
      --chimNonchimScoreDropMin ~{chimNonchimScoDMin} \
      ~{"--chimOutJunctionFormat " + chimOutJunForm} \
      --peOverlapNbasesMin ~{peOvNbasesMin} \
      --peOverlapMMp ~{peOvMMp} \
      --outFilterMultimapNmax ~{outFilterMultimapNmax} \
      --runThreadN ~{threads} --chimOutType ~{chimOutType} \
      --chimScoreDropMax ~{chimScoreDropMax} \
      --chimScoreSeparation ~{chimScoreSeparation} \
      --chimSegmentReadGapMax ~{chimSegmentReadGapMax} ~{addParam}

 awk 'NR<2{print $0;next}{print $0| "sort -V"}' ~{outputFileNamePrefix}.~{chimericjunctionSuffix}.junction \
 > tmp && mv tmp ~{outputFileNamePrefix}.~{chimericjunctionSuffix}.junction

>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  cpu:     "~{threads}"
  timeout: "~{timeout}"
}

output {
 File outputBam        = "~{outputFileNamePrefix}.~{starSuffix}.bam"
 File outputChimeric   = "~{outputFileNamePrefix}.~{chimericjunctionSuffix}.junction"
 File transcriptomeBam = "~{outputFileNamePrefix}.~{transcriptomeSuffix}.bam"
 File geneReads        = "~{outputFileNamePrefix}.~{genereadSuffix}.tab"
}

meta {
  output_meta: {
    outputBam:        "Output bam aligned to genome",
      outputChimeric:   "Output chimeric junctions file",
    transcriptomeBam: "Output bam aligned to transcriptome",
    geneReads:        "Output raw read counts per transcript"
  }
}

}

# ==========================================
#  TASK 2 of 2: index bam file with picard
# ==========================================
task indexBam {
input {
    File  inputBam
  Int   jobMemory = 12
  String modules  = "picard/2.19.2"
  Int timeout     = 48
}

parameter_meta {
 inputBam:  "Input bam file"
 jobMemory: "Memory allocated indexing job"
 modules:   "modules for running indexing job"
 timeout:   "hours before task timeout"
}

command <<<
 java -Xmx~{jobMemory-6}G -jar $PICARD_ROOT/picard.jar BuildBamIndex \
                              VALIDATION_STRINGENCY=LENIENT \
                              OUTPUT="~{basename(inputBam, '.bam')}.bai" \
                              INPUT=~{inputBam}
>>>

runtime {
   memory: "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
}

output {
  File outputBai = "~{basename(inputBam, '.bam')}.bai"
}

meta {
  output_meta: {
    outputBai: "Output index file for bam aligned to genome"
  }
}

}
