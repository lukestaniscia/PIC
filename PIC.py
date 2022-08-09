# PIC Algorithm Command Line Executor/Wrapper
# By: Luke Staniscia


import sys
import os
from PICcompression import *
from PICdecompression import *
from PICstatistics import *

def compress(options, paths):
	getCompressionStatistics = False
	if options in ["-v", "-kv"]:
		getCompressionStatistics = True
	try:
		for path in paths:
			compressionStatistics = PICcompress(path, returnStatistics = getCompressionStatistics)
			print("##########$$$$$$$$$$##########$$$$$$$$$$##########")
			if getCompressionStatistics == True:
				compressionTimeMin = math.floor(compressionStatistics[1]/60)
				compressionTimeSec = round(compressionStatistics[1] % 60,1)
				print("Compression Savings: " + str(compressionStatistics[0]) + "%")
				print("Compression Time: " + str(compressionTimeMin) + ":" + str(compressionTimeSec) + " min:sec")
				for i in range(len(compressionStatistics[2])):
					print("Image Space used in Image #" + str(i + 1) + ": " + str(compressionStatistics[2][i]) + "%")
			if options not in ["-k", "-kv"]:
				os.remove(path)
	except: 
		print("Error. " + path + " does not exist. Try again.")
		sys.exit()

def decompress(options, paths):
	getDecompressionStatistics = False
	if options in ["-dv", "-dkv"]:
		getDecompressionStatistics = True
	try:
		for path in paths:
			decompressionStatistics = PICdecompress(path, returnStatistics = True)
			print("##########$$$$$$$$$$##########$$$$$$$$$$##########")
			if getDecompressionStatistics == True:
				decompressionTimeMin = math.floor(decompressionStatistics[1]/60)
				decompressionTimeSec = round(decompressionStatistics[1] % 60,1)
				print("Decompression Time: " + str(decompressionTimeMin) + ":" + str(decompressionTimeSec) + " min:sec")
			if options not in ["-dk", "-dkv"]:
				os.remove(path + "_meta.txt")
				os.remove(path + "_parameters.bin")
				for i in range(decompressionStatistics[0]):
					os.remove(path + "_img_" + str(i + 1) + ".png")
	except:
		print("Error. " + path + " does not exist. Try again.")
		sys.exit()

def advancedStatistics(paths, getStatistics):
	try: 
		if getStatistics == False:
			for path in paths:
				compressionStatistics = PICcompress(path, returnStatistics = True)
				compressionStatistics = compressionStatistics[1:] #remove compression savings for the statistics function
				print("##########$$$$$$$$$$##########$$$$$$$$$$##########")
				decompressionStatistics = PICdecompress(path[:len(path)-4], returnStatistics = True)
				decompressionStatistics = decompressionStatistics[1:] #return number of images for the statistics function
				print("##########$$$$$$$$$$##########$$$$$$$$$$##########")
				PICstatistics(path, compressionStatistics = compressionStatistics, decompressionStatistics = decompressionStatistics)
				print("##########$$$$$$$$$$##########$$$$$$$$$$##########")
		else:
			compressionStatistics = PICcompress(paths, returnStatistics = True, constructViewableImage = True)
			compressionStatistics = compressionStatistics[1:]
			print("##########$$$$$$$$$$##########$$$$$$$$$$##########")
			decompressionStatistics = PICdecompress(paths[:len(paths)-4], returnStatistics = True)
			decompressionStatistics = decompressionStatistics[1:]
			print("##########$$$$$$$$$$##########$$$$$$$$$$##########")
			return PICstatistics(paths, compressionStatistics = compressionStatistics, decompressionStatistics = decompressionStatistics)
	except:
	 	print("Error. One or more paths are invalid. Try again.")
	 	sys.exit()

args = sys.argv[1:]
if args[0][0] == "-": #options specificed by user
	options = args[0]
	paths = args[1:]
else: #no options specified by the user - assume compression
	options = ""
	paths = args
for path in paths:
	if options in ["", "-k", "-v", "-kv", "-s"]:
		if path[len(path)-4:] != ".pdb" and path[len(path)-4:] != ".cif":
			print("One or more paths are invalid. Try again.")
			sys.exit()
	if options[:2] == "-r":
		if path[len(path)-1] != "/":
			print("Invalid path. Try again.")
			sys.exit()
if options in ["", "-k", "-v", "-kv"]:
	compress(options, paths)
elif options in ["-d", "-dk", "-dv", "-dkv"]:
	decompress(options, paths)
elif options == "-s":
	advancedStatistics(paths, False)
elif options in ["-r", "-rk", "-rv", "-rkv", "-rd", "-rdk", "-rdv", "-rdkv"]: 
	foundPaths = []
	if len(options) == 2 or options[:3] != "-rd": #compress .pdb and/or .cif file(s) in the specified directory
		options = options[2:]
		if len(options) > 0: #additional options specified besides simple compress
			options = "-" + options
		for item in os.listdir(paths[0]):
			if item[len(item)-4:] == ".pdb" or item[len(item)-4:] == ".cif":
				foundPaths = foundPaths + [paths[0] + item]
		if foundPaths == []:
			print("No .pdb or .cif files to be compressed.")
			sys.exit()
		else:
			compress(options, foundPaths)
	else: #decompress .pdb and/or .cif file(s) in the specified directory
		options = "-" + options[2:]
		for item in os.listdir(paths[0]):
			if len(item) >= 15 and item[len(item)-15:] == "_parameters.bin":
				foundPaths = foundPaths + [paths[0] + item[:len(item)-15]]
		if foundPaths == []:
			print("No files to be decompressed.")
			sys.exit()
		else: 
			try:
				decompress(options, foundPaths)
			except:
				print("Some compression files are missing. Re-compress the .pdb and/or .cif file(s) and try again.")
				sys.exit()
elif options == "-e":
	import urllib.request
	import shutil
	import ssl
	import matplotlib.pyplot as plt
	import numpy as np
	import math

	print("CREATING TEST PROTEIN DIRECTORY AND DOWNLOADING NECESSARY FILES FROM THE PROTEIN DATA BANK")
	ssl._create_default_https_context = ssl._create_unverified_context #allow the download of files from websites with no certificates
	newPath = 'Test Proteins/' 
	if not os.path.exists(newPath):
		print("Creating " + "'" + newPath + "'" + " path")
		os.makedirs(newPath)
	proteinIDs = ["2ja9", "2jan", "2jbp", "2ja8", "2ign", "2jd8", "2ja7", "2fug", "2b9v", "2j28", "6hif", "3j7q", "3j9m", "6gaw", "5t2a", "4ug0", "4v60", "4wro", "6fxc", "4wq1"]
	for i in range(10):
		newPath = 'Test Proteins/' + proteinIDs[i] + '/'
		if not os.path.exists(newPath):
			print("Creating " + "'" + newPath + "'" + " path")
			os.makedirs(newPath)
		if not os.path.exists(newPath + proteinIDs[i] + ".pdb"):
			url = "https://files.rcsb.org/download/" + proteinIDs[i].upper() + ".pdb"
			print("Downloading " + url)
			with urllib.request.urlopen(url) as response, open(newPath + proteinIDs[i] + ".pdb", 'wb') as out_file:
				shutil.copyfileobj(response, out_file)
	for i in range(10):
		newPath = 'Test Proteins/' + proteinIDs[i+10] + '/'
		if not os.path.exists(newPath):
			print("Creating " + "'" + newPath + "'" + " path")
			os.makedirs(newPath)
		if not os.path.exists(newPath + proteinIDs[i+10] + ".cif"):
			url = "https://files.rcsb.org/download/" + proteinIDs[i+10].upper() + ".cif"
			print("Downloading " + url)
			with urllib.request.urlopen(url) as response, open(newPath + proteinIDs[i+10] + ".cif", 'wb') as out_file:
				shutil.copyfileobj(response, out_file)

	print("COMPRESSING TEST PROTEIN FILES")
	filenameExtensions = [".pdb" for i in range(10)] + [".cif" for i in range(10)]
	LATEXOutputs = []
	numAtoms = []
	txtSizes = []
	gZipSizes = []
	pngSizes = []
	for i in range(len(proteinIDs)):
		path = "Test Proteins/" + proteinIDs[i] + "/" + proteinIDs[i] + filenameExtensions[i]
		results = advancedStatistics(path, True)
		print("##########$$$$$$$$$$##########$$$$$$$$$$##########")
		LATEXOutputs = LATEXOutputs + [results[0]]
		numAtoms = numAtoms + [results[1][0]]
		txtSizes = txtSizes + [results[1][1]]
		gZipSizes = gZipSizes + [results[1][2]]
		pngSizes = pngSizes + [results[1][3]]

	print("TABLE ONE DATA")
	for latex in LATEXOutputs:
		print(latex)

	print("COMPUTING FIGURES 3 AND 4")
	gZipCRs = []
	pngCRs = []
	maxpngCR = 0
	maxSavings = 0
	for i in range(len(proteinIDs)):
		gZipCR = round(txtSizes[i]/gZipSizes[i],3) #computing compression ratios
		pngCR = round(txtSizes[i]/pngSizes[i],3)
		if pngCR > maxpngCR:
			maxpngCR = pngCR
		gZipCRs = gZipCRs + [gZipCR]
		pngCRs = pngCRs + [pngCR]
		savings = round((1 - pngSizes[i]/gZipSizes[i])*100,1)
		if savings > maxSavings:
			maxSavings = savings

	print("PLOTING FIGURE 3")
	plt.plot(numAtoms, gZipCRs, "r", label = "gZip")
	plt.plot(numAtoms, pngCRs, "b", label = "PIC")
	plt.xlabel("Number of Atoms")
	plt.ylabel("Compression Ratio")
	plt.title(" Compression Ratios vs. Protein Size")
	plt.legend()
	plt.savefig("CRvsnumAtoms")
	plt.close()

	print("PLOTTING FIGURE 4")
	plt.plot(gZipCRs, pngCRs, "o")
	plt.plot([0., maxpngCR], [0., maxpngCR], "black") #plot diagonal line 
	plt.xlabel("gZip Compression Ratio")
	plt.ylabel("PIC Compression Ratio")
	plt.title("PIC vs. gZip Compression Ratios")
	plt.savefig("pngCRvsgZipCR")
	plt.close()

	print("Max Savings: " + str(maxSavings) + "%")
else:
	print("Invalid option. Try again.")
