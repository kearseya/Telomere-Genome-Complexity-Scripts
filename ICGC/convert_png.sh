#!/bin/bash

for i in {1..284}; do gs -o high_names/${i}.png -sDEVICE=pngalpha -r720 -g6000x20000 -dFirstPage=${i} -dLastPage=${i} -f high/names/high_br_names.pdf; done
# -c "<</Install {-118 -530 translate}>> setpagedevice" 
