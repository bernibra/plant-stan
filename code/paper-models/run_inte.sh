#!/usr/bin/env bash

bsub -R "rusage[mem=4548]" -W 120:00 -n 15 "R --vanilla --slave < 1d-baseline.R"
