README

NAME
PIC Compression and Decompression Algorithm

DESCRIPTION
Compression and decompression algorithm for 3D Cartesian coordinates contained in PDB and mmCIF files. See the article "PIC: an image-centric compression algorithm for structural protein data" for further details. 

VISUALS
See the aforementioned article for sample outputs from the compression algorithm. 

INSTALATION
Place the PIC.py, PICcompression.py, PICdecompression.py, and PICstatistics.py files in your working directory.

REQUIREMENTS
The aforementioned files that make up the PIC algorithm are written in Python V3. This implementation makes use of the Python packages math, PIL.{Image, ImageOps, ImageEnhance}, os, time, sys, gzip, shutil, urllib.request, ssl, matplotlib.pyplot, and numpy. Python V3 and the aforementioned libraries must be installed on your system before using the PIC algorithm.

USAGE
Similar to the gZip command line tool but accommodates additional options. Run "python3 PIC.py" in your command terminal followed by one of the following options and one or more path(s) or director(y/ies).

			Compress .pdb/.cif files at the following path(s) and delete original files.
	-k		Compress .pdb/.cif files at the following path(s) and keep original files.
	-v 		Compress .pdb/.cif files at the following path(s) and get compression statistics (compression savings, compression time, image space used).
	-kv		Compress .pdb/.cif files at the following path(s), keep original files, and get compression statistics (compression savings, compression time, image space used).
	-d		Decompress .pdb/.cif files at the following director(y/ies) and delete compressed files.
	-dk		Decompress .pdb/.cif files at the following director(y/ies) and keep compressed files.
	-dv 	Decompress .pdb/.cif files at the following director(y/ies) and get decompression statistics (decompression time).
	-dkv	Decompress .pdb/.cif files at the following director(y/ies), keep compressed files, and get decompression statistics (decompression time).
	-s		Get advanced compression and decompression statistics for .pdb/.cif files at the following path(s). Keeps all original and compressed files.
	-r 		Compresses all .pdb/.cif files in the following directory. 
	-rk		Compresses all .pdb/.cif files in the following directory and keep original files.
	-rv 	Compresses all .pdb/.cif files in the following directory and get compression statistics (compression savings, compression time, image space used).
	-rkv 	Compresses all .pdb/.cif files in the following directory, keep original files, and get compression statistics (compression savings, compression time, image space used).
	-rd 	Decompress all .pdb/.cif files in the following directory. 
	-rdk	Decompress all .pdb/.cif files in the following directory and keep compressed files.
	-rdv 	Decompress all .pdb/.cif files in the following directory and get decompression statistics (decompression time).
	-rdkv 	Decompress all .pdb/.cif files in the following directory, keep compressed files, and get decompression statistics (decompression time).
	-e		Reproduce experiments in the aforementioned article. No further arguments should be given, they will be ignored. Keeps all original and compressed files. 

All paths and directories should be given with respect to the working directory. All outputs will be saved at the location of the original or compressed file(s), with the exception of the figures produced when using the option "-e", which will be saved at the working directory. For options beginning with "-d", only director(y/ies), not file(s), should be specified. The director(y/ies) should be up to the filename of the compressed .pdb/.cif file(s) but not including the filename extension(s). Advanced compression and decompression statistics produced when running using the option "-s" are for the 3D Cartesian coordinates only, not the whole file. Whole file compression and decompression statistics are outputted when running using the options "-v" or "-kv" and "-dv" or "-dkv" respectively. 

Note that it is possible to change the epsilon value of the algorithm, as described in the aforementioned article, as well as the printout frequency by changing the epsilon and notificationFrequency variables in the final function definition in each of the files PICcompression.py, PICdecompression.py, and PICstatistics.py. Note that use of values of epsilon other than those discussed in the article may cause errors in the execution of the algorithm.

SUPPORT
Please email luke.staniscia@mail.utoronto.ca if you have any questions or comments.

AUTHORS AND ACKNOWLEDGMENT
Code written by Luke Staniscia. Luke would like to thank Prof. Yun William Yu for his guidance throughout the course of this project, as well as Jim Shaw, Ziye Tao, and Alex Leighton for their thought-provoking questions that inspired parts of the aforementioned article. 

LICENSE
Freely available to non-commercial users.

PROJECT STATUS
Complete.