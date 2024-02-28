#!/bin/bash
for i in *; do tesseract ${i} $(basename ${i} .png) --psm 4; done
