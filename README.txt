README

NAME
PIC Compression and Decompression Algorithm

DESCRIPTION
Compression and decompression algorithm for 3D Cartesian coordinates contained in PDB and mmCIF files. See the article "PIC: an image centric compression algorithm for structural protein data" for further details. 

VISUALS
See the aforementioned article for sample outputs from the compression algorithm. 

INSTALATION
Place the PIC.py, PICcompression.py, PICdecompression.py, and PICstatistics.py files in your working directory.

REQUIREMENTS
The aforementioned files that make up the PIC algorithm are written in Python V3. This implementation makes use of the Python packages Bio.PDB, math, PIL.{Image, ImageOps, ImageEnhance}, os, time, sys, gzip, shutil, matplotlib.pyplot, and numpy. Python V3 and the aforementioned libraries must be installed on your system before using the PIC algorithm.

USAGE
Run "python3 PIC.py" in your command terminal and follow the prompts to compress, decompress, or gather statistics on the compression and decompression of a PDB or mmCIF file, as well as to reproduce the results found in the aforementioned paper. The PDB or mmCIF file being compressed should be located at your working directory, or in a directory downstream. All outputs will be saved at the location of the PDB or mmCIF file, with the exception of the figures produced when recreating the results found in the aforementioned article, which will be saved at the working directory.

Note that it is possible to change the epsilon value of the algorithm, as described in the aforementioned article, as well as the printout frequency by changing the epsilon and notificationFrequency variables in the final function definition in each of the files PICcompression.py, PICdecompression.py, and PICstatistics.py. Note that use of values of epsilon other than those discussed in the article may cause errors in the execution of the algorithm.

SUPPORT
Please email luke.staniscia@mail.utoronto.ca if you have any questions or comments.

AUTHORS AND ACKNOWLEDGMENT
Code written by Luke Staniscia. Luke would like to thank Prof. Yun William Yu for his guidance throughout the course of this project, as well as his research group for their thought provoking questions. 

LICENSE
Freely available to non-commercial users.

PROJECT STATUS
Complete.