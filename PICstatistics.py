# PIC Algorithm Compression and Decompression Statistics Gatherer
# By: Luke Staniscia

#import used libraries/packages
from Bio.PDB import *
import math
import sys
from PIL import Image, ImageOps, ImageEnhance
import gzip
import shutil
import os
import time

def cart2sph(cord, precision = 1):
	x = cord[0]
	y = cord[1]
	z = cord[2]
	XsqPlusYsq = x**2 + y**2
	r = round(math.sqrt(XsqPlusYsq + z**2) * (10**precision)) #integer encoding spherical coordinates
	az = math.atan2(y,x) * (180/math.pi)   
	if az < 0:
		az = az + 360
	az = round(az * (10**precision)) % (360 * (10**precision))
	elev = round(abs(math.atan2(math.sqrt(XsqPlusYsq),z) * (180/math.pi)) * (10**precision))
	return [az, elev, r]

def translate(cords, precision = 3):
	newCord = [0,0,0]
	for j in range(3):
		newCord[j] = round(cords[j] - globalCentroid[j],precision)
	return newCord

def cart2plan(cartCords):
	sphCords = cart2sph(translate(cartCords))
	return [sphCords[2], round(bytesPerTenthDegree * sphCords[0] * 8), sphCords[1]]

def toBits(x):
	binX = "10" #encoding data
	newBin = bin(abs(x))[2:] #eliminate "0b" prefix
	while len(newBin) < math.ceil(math.log(1801,2)):
		newBin = "0" + newBin
	binX = binX + newBin
	return binX

def keyy(x):
	maxBits = round(3600*bytesPerTenthDegree*8)
	return (x[0]*maxBits + x[1])*1800 + int(x[2][2:],2)

def reverseTextLine(s):
	i = 0
	floats = []
	for j in range(3):
		currentFloat = ""
		while s[i] != "\t" and s[i] != "\n":
			currentFloat = currentFloat + s[i]
			i = i + 1
		floats = floats + [float(currentFloat)]
		i = i + 1
	return floats

def PICstatistics(proteinID, directory, isPDBFile, epsilon = 0.25, compressionStatistics = [], decompressionStatistics = [], notificationFrequency = 500):

	print("##########$$$$$$$$$$########## COMPUTING STATISTICS ON THE COMPRESSION/DECOMPRESSION OF " + proteinID + " ##########$$$$$$$$$$##########")

	global bytesPerTenthDegree
	bytesPerTenthDegree = epsilon

	print("READING PARAMETER FILE")
	parameterFile = open(directory + "_parameters.bin", "rb")
	parameterBitsRead = ""
	byte = parameterFile.read(1)
	while byte:
		newBits = bin(byte[0])[2:]
		while len(newBits) < 8:
			newBits = "0" + newBits
		parameterBitsRead = parameterBitsRead + newBits
		byte = parameterFile.read(1)
	parameterFile.close()

	firstByte = True
	tracker = 1
	parameters = []
	currentParameterBits = ""
	while len(parameterBitsRead) > 0:
		currentBits = parameterBitsRead[:8]
		parameterBitsRead = parameterBitsRead[8:]
		if currentBits[0] == "1" and firstByte == False: #new value started
			print("Reading Parameter " + str(tracker + 1))
			parameter = int(currentParameterBits[1:],2)
			if currentParameterBits[0] == "1":
				parameter = (-1)*parameter
			parameters = parameters + [parameter]
			currentParameterBits = currentBits[1:]
			tracker = tracker + 1
		else:
			currentParameterBits = currentParameterBits + currentBits[1:]
			firstByte = False
	print("Reading Parameter " + str(tracker + 1))
	parameter = int(currentParameterBits[1:],2)
	if currentParameterBits[0] == "1":
		parameter = (-1)*parameter
	parameters = parameters + [parameter]

	print("Adjusting and Assiging Parameters")
	for i in range(3):
		parameters[len(parameters)-1-i] = parameters[len(parameters)-1-i]/1000 #converting integer encoded values back to floats
	numImages = parameters[0]
	croppingParameters = parameters[1:len(parameters)-3]
	global globalCentroid
	globalCentroid = parameters[len(parameters)-3:len(parameters)]

	print("READING PROTEIN FILES")
	tracker = 0
	maxCord = 0
	numAtoms = 0
	if isPDBFile == True:
		parser = PDBParser()
		structure = parser.get_structure(proteinID, directory + ".pdb")
		orgSize = os.path.getsize(directory + ".pdb")
	else:
		parser = MMCIFParser()
		structure = parser.get_structure(proteinID, directory + ".cif")
		orgSize = os.path.getsize(directory + ".cif")
	for atom in structure.get_atoms():
		if (tracker + 1) % notificationFrequency == 0:
			print("Reading Atom " + str(tracker + 1))
		cartCords = atom.get_coord()
		for j in range(3):
			if abs(cartCords[j]) > maxCord:
				maxCord = abs(cartCords[j])
		numAtoms = numAtoms + 1
		tracker = tracker + 1 
	maxCord = math.ceil(maxCord*10)

	print("CONSTRUCTING ORIGINAL COORDINATES AND WRITING BINARY FILE")

	print("Initalizing Binary File")
	binFile = open(directory + '_bin.bin','wb')
	binFileBits = ""
	newBits = bin(maxCord)[2:] #variably packing maxCord so the binary file can be decoded
	while len(newBits) % 7 != 0:
		newBits = "0" + newBits
	while len(newBits) > 7:
		binFileBits = binFileBits + "0" + newBits[:7]
		newBits = newBits[7:]
	binFileBits = binFileBits + "1" + newBits

	tracker = 0
	toTextQueue = []
	for atom in structure.get_atoms():
		if (tracker + 1) % notificationFrequency == 0:
			print("Queueing and Writing Binary Coordinates for Atom " + str(tracker + 1) + " of " + str(numAtoms))
		cartCords = atom.get_coord()
		for i in range(3):
			cartCords[i] = round(cartCords[i],3)
		for j in range(3):
			newBits = bin(math.trunc(round(cartCords[j]*10)))
			if newBits[0] == "-":
				newBits = newBits[3:]
				while len(newBits) < math.ceil(math.log(maxCord + 1,2)):
					newBits = "0" + newBits
				newBits = "1" + newBits
			else:
				newBits = newBits[2:]
				while len(newBits) < math.ceil(math.log(maxCord + 1,2)):
					newBits = "0" + newBits
				newBits = "0" + newBits
			binFileBits = binFileBits + newBits
		if len(binFileBits) % 8 == 0:
			Bytes = []
			while len(binFileBits) > 0:
				Bytes = Bytes + [int(binFileBits[:8],2)]
				binFileBits = binFileBits[8:]
			binFileByteArray = bytearray(Bytes)
			binFile.write(binFileByteArray)
		cords = cart2plan(cartCords)
		toTextQueue = toTextQueue + [[cords[0],cords[1],toBits(cords[2]),cartCords]]
		tracker = tracker + 1

	print("Finishing Wiriting Binary File")
	Bytes  = []
	while len(binFileBits) > 8:
		Bytes = Bytes + [int(binFileBits[:8],2)]
		binFileBits = binFileBits[8:]
	while len(binFileBits) % 8 != 0:
		binFileBits = binFileBits + "0" #adding zero bits to the end as entries are of binary(maxCord + 1) length 
	if len(binFileBits) > 0:
		Bytes = Bytes + [int(binFileBits,2)]
	binFileByteArray = bytearray(Bytes)
	binFile.write(binFileByteArray)
	binFile.close()

	print("CONSTRUCING TEXT FILE")
	toTextQueue.sort(key = keyy) #sort queue so data points are in the same order as how written in the dcompressed file
	
	tracker = 0
	textFile = open(directory + '_txt.txt','w+')
	for data in toTextQueue:
		if (tracker + 1) % notificationFrequency == 0:
			print("Writing Data Point " + str(tracker + 1) + " of " + str(len(toTextQueue)) + " to Text File")
		textFile.write(str(round(data[3][0],1)) + '\t' + str(round(data[3][1],1)) + '\t' + str(round(data[3][2],1)) + '\n')
		tracker = tracker + 1
	textFile.close()

	print("COMPARING DECOMPRESSED FILE")
	tracker = 0
	decompressedFile = open(directory + "_decompressed_txt.txt", "r")
	numberMismatches = 0
	for data in toTextQueue:
		if (tracker + 1) % notificationFrequency == 0:
			print("Decompressing Data Point " + str(tracker + 1) + " of " + str(len(toTextQueue)))
		decompressedLine = decompressedFile.readline()
		originalCords = data[3]
		decompressedCords = reverseTextLine(decompressedLine)
		euclideanDistance = 0
		for i in range(3):
			euclideanDistance = euclideanDistance + (originalCords[i] - decompressedCords[i])**2
		euclideanDistance = math.sqrt(euclideanDistance)
		if euclideanDistance > 0.2*math.sqrt(3):
			print("MISMATCH!")
			print(originalCords)
			print(decompressedCords)
			numberMismatches = numberMismatches + 1
		tracker = tracker + 1
	decompressedFile.close()
	if numberMismatches == 0:
		print("DECOMPRESSION SUCESSFUL")
	else:
		print("DECOMPRESSION NOT SUCESSFUL; Number of Mismatches: " + str(numberMismatches))

	print("COMPRESSING BINARY FILE AND RETREVING FILE SIZES")
	print("Retreving Text, Binary, and PNG File Sizes")
	textSize = os.path.getsize(directory + '_txt.txt')
	binSize = os.path.getsize(directory + '_bin.bin')
	pngSize = 0
	for i in range(numImages):
		pngSize = pngSize + os.path.getsize(directory + "_img_" + str(i + 1) + ".png")
	pngSize = pngSize + os.path.getsize(directory + "_parameters.bin")
	print("Compressing Binary File with gZip")
	with open(directory + '_bin.bin', 'rb') as f_in, gzip.open(directory + '_bin.bin.gz', 'wb') as f_out:
	    shutil.copyfileobj(f_in, f_out)
	gZipSize = os.path.getsize(directory + '_bin.bin.gz')

	print("COMPUTING COMPRESSION RATIOS")
	gZipCR = round(textSize/gZipSize,3)
	pngCR = round(textSize/pngSize,3)

	if len(compressionStatistics) > 0:
		print("COMPUTING COMPRESSION TIME STATISTICS")
		compressionTimeMin = math.floor(compressionStatistics[0]/60)
		compressionTimeSec = round(compressionStatistics[0] % 60,1)

	if len(decompressionStatistics) > 0:
		print("COMPUTING DECOMPRESSION TIME STATISTICS")
		decompressionTimeMin = math.floor(decompressionStatistics[0]/60)
		decompressionTimeSec = round(decompressionStatistics[0] % 60,1)

	print("###### COMPRESSION STATISTICS START #####")
	print("Protein ID: " + proteinID)
	print("Atom Count: " + str(numAtoms))
	print("###")
	print("Original File Size: " + str(round(orgSize/1000,1)) + " KB")
	print("Rounded Coordinate Text File Size: " + str(round(textSize/1000,1)) + " KB")
	print("Rounded Coordinate Binary File Size: " + str(round(binSize/1000,1)) + " KB")
	print("###")
	print("gZip Compression")
	print("Size: " + str(round(gZipSize/1000,1)) + " KB")
	print("Ratio: " + str(gZipCR))
	print("###")
	print("PIC Compression")
	print("Size: " + str(round(pngSize/1000,1)) + " KB")
	print("Ratio: " + str(pngCR))
	print("###")
	if len(compressionStatistics) > 0:
		print("Compression Time: " + str(compressionTimeMin) + ":" + str(compressionTimeSec) + " min:sec")
	if len(decompressionStatistics) > 0:
		print("Deompression Time: " + str(decompressionTimeMin) + ":" + str(decompressionTimeSec) + " min:sec")
	print("Number of Images used in Compression: " + str(numImages))
	for i in range(numImages):
		print("###")
		print("For Image " + str(i + 1))
		if len(compressionStatistics) > 0:
			print("Proportion of Available Image Space used: " + str(compressionStatistics[1][i]) + "%")
		print("Horizontal Cropping Parameter: " + str(croppingParameters[2*i]))
		print("Vertical Cropping Parameter: " + str(croppingParameters[2*i+1]))
	print("##### COMPRESSION STATISTICS END #####")

	print("DATA LATEX PRINTOUT")

	print("Formatting LATEX Output")
	if pngCR > gZipCR:
		gZipSizeString =  str(round(gZipSize/1000,1))
		gZipCRString =  str(gZipCR)
		pngSizeString = "\\textbf{" + str(round(pngSize/1000,1)) + "}"
		pngCRString = "\\textbf{" + str(pngCR) + "}"
	else: 
		gZipSizeString = "\\textbf{" + str(round(gZipSize/1000,1)) + "}"
		gZipCRString = "\\textbf{" + str(gZipCR) + "}"
		pngSizeString = str(round(pngSize/1000,1))
		pngCRString = str(pngCR)
	latex = proteinID + " & " + str(numAtoms) + " & " + str(round(orgSize/1000,1)) + " & " + str(round(textSize/1000,1)) + " & " + str(round(binSize/1000,1)) +  " & " + gZipSizeString + " & " + gZipCRString + " & " + pngSizeString + " & " + pngCRString
	if len(compressionStatistics) > 0:
		latex = latex +  " & " + str(compressionTimeMin) + ":" + str(compressionTimeSec) + " & " + str(numImages) + " & " + str(compressionStatistics[1]) 
	else: 
		latex = latex + " & " + str(numImages)
	if len(decompressionStatistics) > 0:
		latex = latex + " & " + str(decompressionTimeMin) + ":" + str(decompressionTimeSec) + " \\\\"
	else: 
		latex = latex + " \\\\"
	print(latex)
	return [latex, [numAtoms, textSize, gZipSize, pngSize]]