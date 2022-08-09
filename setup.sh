#!/bin/sh
python3 setup_compression.py build_ext --inplace
python3 setup_decompression.py build_ext --inplace
