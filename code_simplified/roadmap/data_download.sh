#!/bin/bash

DL_PATH=/fast/AG_Ohler/hyliou/mtproject/raw/roadmap/

# -i:read url from a file
# -P the prefix where all files are saved
# -nc: no clobber, overwrite the existed files with the same file names
wget -i ${DL_PATH}url-gappedPeak.txt -P ${DL_PATH} -nc