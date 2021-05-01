#!/usr/bin/env bash

bsub -R "rusage[mem=2840]" -W 200:00 -n 45 "R --vanilla --slave < 1d-baseline.R"
