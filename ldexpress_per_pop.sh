#!/bin/bash

while read population; do
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Processing: $population"
    Rscript ldexpress_script_appendage.R "$population"
done < populations.txt