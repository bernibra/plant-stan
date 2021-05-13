#!/usr/bin/env bash

bsub -R "rusage[mem=3500]" -W 200:00 -n 35 "R --vanilla --slave < 1d-categorical.R"
