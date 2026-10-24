#!/bin/bash
cd "$(dirname "$0")"
while [ $(pgrep -f "run_lava.R info_corrected" | wc -l) -ge 3 ]; do sleep 20; done
Rscript run_lava.R info_corrected.tsv ov_0.1724.tsv p1_c3 1873 2495 > log_p1_c3.txt 2>&1
