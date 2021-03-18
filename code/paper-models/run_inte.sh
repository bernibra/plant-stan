#!/usr/bin/env bash

bsub -R "rusage[mem=4300]" -W 200:00 -n 30 "R --vanilla --slave < 1d-baseline.R"
