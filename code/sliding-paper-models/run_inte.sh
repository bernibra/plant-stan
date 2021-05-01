#!/usr/bin/env bash

--args \$LSB_JOBINDEX $1 <
bsub -R "rusage[mem=2500]" -W 24:00 -n 35 "R --vanilla --slave --args $1 < 1d-baseline.R"
