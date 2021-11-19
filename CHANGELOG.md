# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2021-10-01
### Added
- Initiated repo
- Added WDL, subworkflows and imports
- Added Vidarr files

## [1.0.1]
## Added
- no additions, retagged when working out some issues

## [1.0.2] - 2021-10-29
## Added
- ability to align with STAR for whole transcriptome, pulled in STAR workflow
- option to select aligner, bwa or star
- modifications to regression tests to put WT libraries through STAR
- modification to inputs to use star aligned bams for WT libraries, in bam mode

## [1.0.3] - 2021-11-19
## Added
- added input parameter maxReads, with default 0
- if maxReads > 0, then a new task is run, downsample, using seqtk to downsample the fastq files prior to alignment

