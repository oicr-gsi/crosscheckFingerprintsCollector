# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0] - 2024-06-25
### Added
- [GRD-797](https://jira.oicr.on.ca/browse/GRD-797) - add vidarr labels to outputs (changes to medata only)

## [1.0.0] - 2021-10-01
### Added
- Initiated repo
- Added WDL, subworkflows and imports
- Added Vidarr files

## [1.0.1]
## Added
- no additions, retagged when working out some issues

## [1.0.2] - 2021-10-29
### Added
- ability to align with STAR for whole transcriptome, pulled in STAR workflow
- option to select aligner, bwa or star
- modifications to regression tests to put WT libraries through STAR
- modification to inputs to use star aligned bams for WT libraries, in bam mode

## [1.0.3] - 2021-11-19
### Added
- added input parameter maxReads, with default 0
- if maxReads > 0, then a new task is run, downsample, using seqtk to downsample the fastq files prior to alignment


## [1.0.4] - 2021-05-12
### Added
- added input markDups, if set to true will run the task markDuplicates on the bam file prior to fingerprinting.  ExtractFingerprint accounts for marked duplicates
- added assessCoverage task to generate a samtools coverage report on the bam file (with -ff default filters), then generated a json with the mean depth
- added option sampleId : the provided sampleId is given to the ExtractFingeprints --SAMPLE_ALIAS argument, and will appear in the vcf fingeprint as the sample identifier.  This allows crosscheck to be run as by=SAMPLE, and will help to provide uniformity to how the fingerprint comparisons are run
- added duplicate marking and sampleID to each of the regression tests as these are both expected to be used regularl


## [1.0.5] - 2021-05-31
### Added
- added samtools stats to generate basic stats on the final bam file
### Changed
- renamed assessCoverage task to alignmentMetrics
- running samtools coverage with and withouth -ff DUP, to assess coverage before and after DUP removal
- capturing multiple stats in the json structure

## [1.0.8] - 2022-07-19
### Added
- Adding java option Xmx to make sure markDupilcates caommand running has enough memory
