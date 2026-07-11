#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

#enter the workflow's final output directory ($1)
cd $1


ls | sed 's/.*\.//' | sort | uniq -c
