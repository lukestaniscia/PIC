# PIC Compressor
# By: Luke Staniscia

#import used libraries/packages
from Bio.PDB import *
import math
from PIL import Image, ImageOps, ImageEnhance
import os
import time

def cart2sph(cord, precision = 1):
	x = cord[0]
	y = cord[1]
	z = cord[2]
	XsqPlusYsq = x**2 + y**2
	r = round(math.sqrt(XsqPlusYsq + z**2) * (10**precision))  #integer encoding spherical coordinates
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

def numPixels2Cords(r, numBits):
	y = math.floor(numBits/8) % round(bytesPerTenthDegree * 3600) 
	return[r, y]

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

def wrapList(aList, start, stop):
	if stop < start:
		return [1 for i in range(len(aList))] #will force the calling loop to continue
	else:
		return aList[start:stop]

def updateImageBoundries(cords):
	global firstDataPoint
	global minX
	global maxX
	global minY
	global maxY
	if firstDataPoint == True:
		minX = cords[0]
		maxX = cords[0]
		minY = cords[1]
		maxY = cords[1]
		firstDataPoint = False
	else:
		if cords[0] < minX:
			minX = cords[0]
		if cords[0] > maxX:
			maxX = cords[0]
		if cords[1] < minY:
			minY = cords[1]
		if cords[1] > maxY:
			maxY = cords[1]

def keyy(x):
	maxBits = round(3600*bytesPerTenthDegree*8)
	return (x[0]*maxBits + x[1])*1800 + int(x[2][2:],2)

def PICcompress(proteinID, directory, isPDBFile, epsilon = 0.25, returnStatistics = False, notificationFrequency = 500):

	print("##########$$$$$$$$$$########## COMPRESSING " + proteinID + " ##########$$$$$$$$$$##########")

	startTime = time.time() #start tracking compression time
	global bytesPerTenthDegree
	bytesPerTenthDegree = epsilon

	print("DATA PRE-PROCESSING")

	print("Computing Global Centroid") 
	global globalCentroid
	globalCentroid = [0,0,0]
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
		cartCords = atom.get_coord()
		for j in range(3):
			globalCentroid[j] = globalCentroid[j] + cartCords[j]
			if abs(cartCords[j]) > maxCord:
				maxCord = abs(cartCords[j])
		numAtoms = numAtoms + 1
	for i in range(3):
		globalCentroid[i] = round(globalCentroid[i]/numAtoms,3)
	maxCord = math.ceil(maxCord*10)

	print("Computing Maximum Global Radius")
	global maxGlobalR
	maxGlobalR = 0
	for atom in structure.get_atoms():
		cords = cart2sph(translate(atom.get_coord()))
		if cords[2] > maxGlobalR:
			maxGlobalR = cords[2]
	global imageHeight
	imageHeight = round(3600*bytesPerTenthDegree)
	global imageWidth
	imageWidth = maxGlobalR + 1

	print("QUEUEING DATA")
	tracker = 0
	queue = [[radius] for radius in range(maxGlobalR+1)]
	for atom in structure.get_atoms():
		if (tracker + 1) % notificationFrequency == 0:
			print("Queing Data for Atom " + str(tracker + 1) + " of " + str(numAtoms))
		cartCords = atom.get_coord()
		for i in range(3):
			cartCords[i] = round(cartCords[i],3)
		cords = cart2plan(cartCords)
		queue[cords[0]] = queue[cords[0]] + [[cords[1], toBits(cords[2])]]
		tracker = tracker + 1

	print("Filtering Queue")
	filteredQueue = []
	for radius in queue:
		if len(radius) > 1:
			filteredQueue = filteredQueue + [radius]
	queue = filteredQueue

	print("MAPPING AND ENCODING DATA")
	Images = []
	occupiedImageSpace = []
	global maxBits
	maxBits = round(3600*bytesPerTenthDegree*8)
	croppingParameters = []
	while len(queue) > 0:
		Images = Images + [[[0 for row in range(imageHeight)] for col in range(imageWidth)]]
		currentImage = len(Images) - 1
		print("Mapping Image " + str(currentImage + 1))
		occupiedImageSpace = occupiedImageSpace + [0]
		global firstDataPoint
		firstDataPoint = True

		tracker = 0
		nextQueue = []
		for radius in queue:
			if (tracker + 1) % notificationFrequency == 0:
				print("Image #" + str(currentImage + 1) + " Mapping: Radius Set " + str(tracker + 1) + " of " + str(len(queue)) + " processed; Image Space used: " + str(round(occupiedImageSpace[currentImage]/(imageHeight*imageWidth*8)*100,1)) + "%; Length of Next Queue: " + str(len(nextQueue)))
			r = radius[0]
			radius = radius[1:]
			nextQueueRadius = [r]
			usedBits = [0 for i in range(maxBits)]
			nextDataPosition = radius[0][0]
			nextData = radius[0][1]
			while len(radius) > 0:
				if wrapList(usedBits, nextDataPosition, (nextDataPosition + len(nextData)) % maxBits) == [0 for i in range(len(nextData))]: #see related article for further details on this nested if-else statment
					newBits = nextData
					radius = radius[1:]
				else:
					adjustedNextDataPosition = nextDataPosition
					traversedBits = 0
					while adjustedNextDataPosition < (nextDataPosition + 8*bytesPerTenthDegree) and wrapList(usedBits, adjustedNextDataPosition, (adjustedNextDataPosition + len(nextData)) % maxBits) != [0 for i in range(len(nextData))]:
						adjustedNextDataPosition = adjustedNextDataPosition + 1
						traversedBits = traversedBits + 1
					if adjustedNextDataPosition < (nextDataPosition + 8*bytesPerTenthDegree):
						nextDataPosition = adjustedNextDataPosition
						newBits = nextData
						radius = radius[1:]
					else:
						delta = 0 #initalize pointer
						while wrapList(usedBits,(adjustedNextDataPosition + delta) % maxBits,(adjustedNextDataPosition + delta + len(nextData) + math.ceil(math.log(maxBits,2))) % maxBits) != [0 for i in range(len(nextData) + math.ceil(math.log(maxBits,2)))] and traversedBits < maxBits:
							delta = delta + 1
							traversedBits = traversedBits + 1
						if  wrapList(usedBits,(adjustedNextDataPosition + delta) % maxBits,(adjustedNextDataPosition + delta + len(nextData) + math.ceil(math.log(maxBits,2))) % maxBits) == [0 for i in range(len(nextData) + math.ceil(math.log(maxBits,2)))] and traversedBits < maxBits:
							nextDataPosition = (adjustedNextDataPosition + delta) % maxBits
							newBin = bin(delta)[2:]
							while len(newBin) < math.ceil(math.log(maxBits,2)):
								newBin = "0" + newBin
							newBits = nextData[0] + "1" + newBin + nextData[2:]
							radius = radius[1:]
						else: 
							nextQueueRadius = nextQueueRadius + radius 
							radius = []
				if len(newBits) > 0:
					for i in range(len(newBits)):
						usedBits[nextDataPosition + i] = 1
					cords = numPixels2Cords(r, nextDataPosition)
					existingStartingBits = bin(Images[currentImage][cords[0]][cords[1]])[2:]
					while len(existingStartingBits) < 8:
						existingStartingBits = "0" + existingStartingBits
					offset = nextDataPosition % 8
					nextDataPositionEnd = nextDataPosition + len(newBits)
					endCords = numPixels2Cords(r, nextDataPositionEnd)
					existingEndingBits = bin(Images[currentImage][endCords[0]][endCords[1]])[2:]
					while len(existingEndingBits) < 8:
						existingEndingBits = "0" + existingEndingBits
					endSet = nextDataPositionEnd % 8
					newBits = existingStartingBits[:offset] + newBits + existingEndingBits[endSet:]
					while len(newBits) > 0:
						Images[currentImage][cords[0]][cords[1]] = int(newBits[:8],2)
						updateImageBoundries(cords)
						newBits = newBits[8:]
						cords[1] = cords[1] + 1
				if len(radius) > 0:
					nextDataPosition = radius[0][0]
					nextData = radius[0][1]
			for bit in usedBits:
				if bit == 1:
					occupiedImageSpace[currentImage] = occupiedImageSpace[currentImage] + 1
			if len(nextQueueRadius) > 1:
				nextQueue = nextQueue + [nextQueueRadius]
			tracker = tracker + 1
		queue = nextQueue

		print("Cropping Coordinates and Saving Cropping Parameters")
		#ensuring cropped image is large enough to compress
		global maxX
		global maxY
		if maxX - minX + 1 < 32:
			maxX = minX + 31
		if maxY - minY + 1 < 32:
			maxY = minY + 31

		croppingParameters = croppingParameters + [[minX, minY]]
		Images[currentImage] = Images[currentImage][minX:maxX+1] #cropping out all black columns
		for i in range(len(Images[currentImage])):
			Images[currentImage][i] = Images[currentImage][i][minY:maxY+1] #cropping out all black rows

	print("CONSTRUCTING IMAGES")
	k = 0 #image tacker
	for image in Images:
		print("Constructing and Saving Image " + str(k + 1) + " of " + str(len(Images)))
		newImageObject = Image.new("L",[len(image),len(image[0])])
		pixels = newImageObject.load()
		for i in range(len(image)):
			for j in range(len(image[0])):
				pixels[i,j] = image[i][j]
		newImageObject.save(directory + "_img_" + str(k + 1) + ".png", optimize = True)
		ImageInvert = ImageOps.invert(newImageObject) #invert and increase contrast on image for easier viewing
		ImageInvertContrasted = ImageEnhance.Contrast(ImageInvert).enhance(5)
		ImageInvertContrasted.save(directory + "_img_" + str(k + 1) + "_forViewing.png", optimize = True)
		k = k + 1

	print("CONSTRUCTING PARAMETER FILE")
	parameters = [len(Images)]
	for croppingParameter in croppingParameters:
		parameters = parameters + [croppingParameter[0], croppingParameter[1]]
	for component in globalCentroid:
		parameters = parameters + [math.trunc(round(component*1000))]
	tracker = 0
	parameterBits = ""
	parameterFile = open(directory + '_parameters.bin','wb')
	for parameter in parameters: #variably pack parameters
		print("Writing Parameter " + str(tracker + 1) + " of " + str(len(parameters)) + " to Binary File")
		if parameter < 0:
			newBits = "1"
		else:
			newBits = "0"
		newBits = newBits + bin(abs(parameter))[2:]
		while len(newBits) % 7 != 0:
			newBits = newBits[0] + "0" + newBits[1:]
		parameterBits = parameterBits + "1" + newBits[:7]
		newBits = newBits[7:]
		while len(newBits) > 0:
			parameterBits = parameterBits + "0" + newBits[:7]
			newBits = newBits[7:]
		Bytes = []
		while len(parameterBits) > 0:
			Bytes = Bytes + [int(parameterBits[:8],2)]
			parameterBits = parameterBits[8:]
		binFileByteArray = bytearray(Bytes)
		parameterFile.write(binFileByteArray)
		tracker = tracker + 1
	parameterFile.close()

	endTime = time.time() #stop tracking compression time
	if returnStatistics == True:
		print("COMPUTING COMPRESSION TIME")
		compressionTime = endTime - startTime

		print("COMPUTING COMPRESSION STATISTICS")
		i = 0
		for image in Images:
			imageWidth = len(image)
			imageHeight = len(image[0])
			occupiedImageSpace[i] = round(occupiedImageSpace[i]/(imageWidth*imageHeight*8)*100,1)
			i = i + 1

		return [compressionTime, occupiedImageSpace]