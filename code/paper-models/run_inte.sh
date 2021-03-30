#!/usr/bin/env bash

bsub -R "rusage[mem=2840]" -W 2:00 -n 45 "R --vanilla --slave < 1d-baseline.R"
