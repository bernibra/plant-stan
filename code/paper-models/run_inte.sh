#!/usr/bin/env bash

bsub -R "rusage[mem=4548]" -W 200:00 -n 15 "R --vanilla --slave < 1d-categorical.R"
