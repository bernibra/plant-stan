#!/usr/bin/env bash

bsub -R "rusage[mem=20048]" -W 8:00 -n 3 "R --vanilla --slave < model-simulated-data.R"

