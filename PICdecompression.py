# PIC Decompressor
# By: Luke Staniscia

#import used libraries/packages
import math
from PIL import Image, ImageOps, ImageEnhance
import time

def keyy(x):
	return (x[0]*maxBits + x[1])*1800 + int(x[2][2:],2)

def sph2Cart(cord, precision = 3):
	az = cord[0]
	elev = cord[1]
	r = cord[2]
	az = az*(math.pi/180)
	elev = elev*(math.pi/180)
	x = round(r*math.sin(elev)*math.cos(az),precision)
	y = round(r*math.sin(elev)*math.sin(az),precision)
	z = round(r*math.cos(elev),precision)
	return [x, y, z]

def reverseTranslate(cord, precision = 3):
	newCord = [0,0,0]
	for i in range(3):
		newCord[i] = round(cord[i] + globalCentroid[i],precision)
	return newCord

def PICdecompress(proteinID, directory, epsilon = 0.25, returnStatistics = False, notificationFrequency = 500):

	print("##########$$$$$$$$$$########## DECOMPRESSING " + proteinID + " ##########$$$$$$$$$$##########")

	startTime = time.time() #start tracking decompression time
	global bytesPerTenthDegree
	bytesPerTenthDegree = epsilon

	print("READING PARAMETER FILE")
	parametersFile = open(directory + "_parameters.bin", "rb")
	parameterBitsRead = ""
	byte = parametersFile.read(1)
	while byte:
		newBits = bin(byte[0])[2:]
		while len(newBits) < 8:
			newBits = "0" + newBits
		parameterBitsRead = parameterBitsRead + newBits
		byte = parametersFile.read(1)
	parametersFile.close()

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

	print("READING IMAGE FILES")
	Images = []
	for k in range(numImages):
		print("Reading Image " + str(k + 1) + " of " + str(numImages))
		imageObject = Image.open(directory + "_img_" + str(k+1) + ".png")
		imageWidth = imageObject.width
		imageHeight = imageObject.height
		Images = Images + [[[0 for row in range(imageHeight)] for col in range(imageWidth)]]
		pixels = imageObject.load()
		for i in range(imageWidth):
			for j in range(imageHeight):
				Images[k][i][j] = pixels[i,j]

	print("REVERSING CROPPING")
	for i in range(numImages):
		print("Reversing Horizontal Cropping on Image " + str(i + 1) + " of " + str(numImages))
		imageHeight = len(Images[i][0])
		Images[i] = [[0 for row in range(imageHeight)] for j in range(croppingParameters[2*i])] + Images[i] #adding back all black columns on the left and right sides of the image(s)
		if i == 0:
			imageWidth = len(Images[0])
		else:
			Images[i] = Images[i] + [[0 for row in range(imageHeight)] for j in range(imageWidth - len(Images[i]))]
		print("Reversing Vertical Cropping on Image " + str(i + 1) + " of " + str(numImages))
		for j in range(imageWidth):
			Images[i][j] = [0 for row in range(croppingParameters[2*i+1])] + Images[i][j] + [0 for row in range(math.floor(3600*bytesPerTenthDegree - imageHeight - croppingParameters[2*i+1]))] #adding back all black rows on the top and bottom sides of the image(s)

	print("EXTRACTING DATA FROM IMAGE")
	tracker = 0
	queue = []
	global maxBits
	maxBits = len(Images[0][0])*8
	for image in Images:
		for i in range(len(image)):
			if (i + 1) % notificationFrequency == 0:
				print("Extracting Data from Image " + str(tracker + 1) + " of " + str(numImages) + "; Column " + str(i + 1) + " of " + str(len(image)))
			data = ""
			for j in range(len(image[0])):
				newBits = bin(image[i][j])[2:]
				while len(newBits) < 8:
					newBits = "0" + newBits
				data = data + newBits
			j = 0
			while len(data) > 0:
				if data[0] == "0":
					data = data[1:]
					j = j + 1
				else:
					if data[1] == "0": #short data point
						queue = queue + [[i, j, data[0:2 + math.ceil(math.log(1801,2))]]]
						data = data[2 + math.ceil(math.log(1801,2)):]
						j = j + 2 + math.ceil(math.log(1801,2))
					else: #long data point
						queue = queue + [[i, j, data[0:2 + math.ceil(math.log(1801,2)) + math.ceil(math.log(maxBits,2))]]]
						data = data[2 + math.ceil(math.log(1801,2)) + math.ceil(math.log(maxBits,2)):]
						j = j + 2 + math.ceil(math.log(1801,2)) + math.ceil(math.log(maxBits,2))
		tracker = tracker + 1

	print("ADJUSTING EXTRACTED DATA")
	tracker = 0
	for data in queue:
		if (tracker + 1) % notificationFrequency == 0:
			print("Adjusting Data Point " + str(tracker + 1) + " of " + str(len(queue)))
		if data[2][1] == "0": #short data point
			data[1] = math.trunc(data[1] - data[1] % (bytesPerTenthDegree*8)) #small adjustment to primary intended position
		else: #long data point
			data[1] = math.trunc((data[1] - int(data[2][2:2 + math.ceil(math.log(maxBits,2))],2) - bytesPerTenthDegree*8) % maxBits) #adjust to primary intented position
			data[2] = data[2][0] + "0" + data[2][2 + math.ceil(math.log(maxBits,2)):] #remove pointer
		tracker = tracker + 1
	queue.sort(key = keyy) #sort queue so data points are in a standardized order

	print("WRITING DECOMPRESSED FILE")
	tracker = 0
	decompressedFile = open(directory + '_decompressed_txt.txt','w+')
	for data in queue:
		if (tracker + 1) % notificationFrequency == 0:
			print("Writing Data Point " + str(tracker + 1) + " of " + str(len(queue)))
		r = round(data[0]/10,1)
		elev = round(int(data[2][2:],2)/10,1) #first two bits are encoding data
		az = round(data[1]/(8*10*bytesPerTenthDegree),1)
		cartCords = reverseTranslate(sph2Cart([az,elev,r]))
		decompressedFile.write(str(round(cartCords[0],3)) + '\t' + str(round(cartCords[1],3)) + '\t' + str(round(cartCords[2],3)) + '\n')
		tracker = tracker + 1
	decompressedFile.close()

	endTime = time.time() #stop tracking decompression time
	if returnStatistics == True:
		print("COMPUTING DECOMPRESSION TIME")
		decompressionTime = endTime - startTime

		return [decompressionTime]