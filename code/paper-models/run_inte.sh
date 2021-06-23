#!/usr/bin/env bash

bsub -R "rusage[mem=3500]" -W 10:00 -n 35 "R --vanilla --slave < 1d-baseline.R"
