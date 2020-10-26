#!/usr/bin/env bash

bsub -R "rusage[mem=10048]" -W 8:00 "R --vanilla --slave < model-simulated-data.R"

