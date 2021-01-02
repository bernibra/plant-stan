#!/usr/bin/env bash

bsub -R "rusage[mem=4548]" -W 100:00 -n 15 "R --vanilla --slave < main.R"
