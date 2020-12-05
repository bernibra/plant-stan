#!/usr/bin/env bash

bsub -R "rusage[mem=7048]" -W 100:00 -n 15 "R --vanilla --slave < main.R"
