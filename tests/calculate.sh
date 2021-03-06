#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

#enter the workflow's final output directory ($1)
cd $1

#find all files, return their md5sums to std out
find . -name *.vcf.gz -xtype f -exec sh -c "zcat {} | grep -v ^# | md5sum" \;
ls | sed 's/.*\.//' | sort | uniq -c
