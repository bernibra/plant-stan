#!/usr/bin/env bash

bsub -R "rusage[mem=68220]" "R --vanilla --slave < compare-models.R"
