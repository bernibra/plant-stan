#!/usr/bin/env bash

bsub -R "rusage[mem=20048]" -W 4:00 -n 15 "R --vanilla --slave < main.R"

