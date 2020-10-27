#!/usr/bin/env bash

bsub -R "rusage[mem=20048]" -W 24:00 -n 3 "R --vanilla --slave < main.R"

