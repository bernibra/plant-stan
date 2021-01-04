#!/usr/bin/env bash

bsub -R "rusage[mem=4548]" -W 100:00 -n 30 "R --vanilla --slave < main.R"
