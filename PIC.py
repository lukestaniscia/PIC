# PIC Algorithm Command Line Executor/Wrapper
# By: Luke Staniscia

print("Welcome to the PIC Compressor!")
print("Would you like to compress (C), decompress (D), get advanced statistics of a protein's compression and decompression (S), or reconstruct results in the article (R)? Enter one of the four previous operation modes.")
mode = ""
while mode == "":
	mode = str(input("Mode: "))
	if mode == "C":
		print("What protein would you like to compress? Enter a single four character protein ID.")
		proteinID = ""
		while proteinID == "":
			proteinID = str(input("Protein ID: "))
			if len(proteinID) != 4:
				print("Invalid protein ID. Please enter a single four character protein ID.")
				proteinID = ""
		print("Is the structural protein data formmated as a PDB (P) or mmCIF file (C)?")
		fileType = ""
		while fileType == "":
			fileType = str(input("File Type: "))
			if fileType == "P":
				isPDBFile = True
				fileType = "pdb"
				print("Enter the directory of the " + proteinID + ".pdb file wrt the current directory. The directory should end in /.")
				directory = str(input("Directory: "))
			elif fileType == "C":
				isPDBFile = False
				fileType = ".cif"
				print("Enter the directory of the " + proteinID + ".cif file wrt the current directory. The directory should end in /.")
				directory = str(input("Directory: "))
			else: 
				print("Invalid file type. Enter P if the structural protein data is formattted as a PDB file or C if the data is formmated as a mmCIF file.")
				fileType = ""
		directory = directory + proteinID
		print("Would you like to know compression statistics? Enter Y if yes or N if no.")
		getCompressionStatistics = ""
		while getCompressionStatistics == "":
			getCompressionStatistics = str(input("Get Compression Statistics?: "))
			if getCompressionStatistics == "Y":
				getCompressionStatistics = True
			elif getCompressionStatistics == "N":
				getCompressionStatistics = False
			else:
				print("Invalid selection. Enter Y if you would like compression statistics or N if not.")
				getCompressionStatistics = ""
		from PICcompression import *
		try: 
			compressionStatistics = PICcompress(proteinID, directory, isPDBFile, returnStatistics = getCompressionStatistics)
			print("##########$$$$$$$$$$##########$$$$$$$$$$##########")
			if getCompressionStatistics == True:
				compressionTimeMin = math.floor(compressionStatistics[0]/60)
				compressionTimeSec = round(compressionStatistics[0] % 60,1)
				print("Compression Time: " + str(compressionTimeMin) + ":" + str(compressionTimeSec) + " min:sec")
				for i in range(len(compressionStatistics[1])):
					print("Image Space used in Image #" + str(i + 1) + ": " + str(compressionStatistics[1][i]) + "%")
		except: 
			print("Error. Make sure you entered all information correctly and the " + proteinID + fileType + " file is in the approriate place then try again.")
	elif mode == "D":
		print("What protein would you like to decompress? Enter a single four character protein ID.")
		proteinID = ""
		while proteinID == "":
			proteinID = str(input("Protein ID: "))
			if len(proteinID) != 4:
				print("Invalid protein ID. Please enter a single four character protein ID.")
				proteinID = ""
		print("Enter the directory of the compressed " + proteinID + " files wrt the current directory. The directory should end in /.")
		directory = str(input("Directory: "))
		directory = directory + proteinID
		print("Would you like to know decompression statistics? Enter Y if yes or N if no.")
		getDecompressionStatistics = ""
		while getDecompressionStatistics == "":
			getDecompressionStatistics = str(input("Get Decompression Statistics?: "))
			if getDecompressionStatistics == "Y":
				getDecompressionStatistics = True
			elif getDecompressionStatistics == "N":
				getDecompressionStatistics = False
			else:
				print("Invalid selection. Enter Y if you would like decompression statistics or N if not.")
				getDecompressionStatistics = ""
		from PICdecompression import *
		try: 
			decompressionStatistics = PICdecompress(proteinID, directory, returnStatistics = getDecompressionStatistics)
			print("##########$$$$$$$$$$##########$$$$$$$$$$##########")
			if getDecompressionStatistics == True:
				decompressionTimeMin = math.floor(decompressionStatistics[0]/60)
				decompressionTimeSec = round(decompressionStatistics[0] % 60,1)
				print("Decompression Time: " + str(decompressionTimeMin) + ":" + str(decompressionTimeSec) + " min:sec")
		except:
			print("Error. Make sure you entered all information correctly and the compressed " + proteinID + " files are in the approriate place then try again.")
	elif mode == "S":
		print("What protein would you like to get advanced compression and decompression statistics on? Enter a single four character protein ID.")
		proteinID = ""
		while proteinID == "":
			proteinID = str(input("Protein ID: "))
			if len(proteinID) != 4:
				print("Invalid protein ID. Please enter a single four character protein ID.")
				proteinID = ""
		print("Is the structural protein data formmated as a PDB (P) or mmCIF file (C)?")
		fileType = ""
		while fileType == "":
			fileType = str(input("File Type: "))
			if fileType == "P":
				isPDBFile = True
				print("Enter the directory of the " + proteinID + ".pdb file and compressed " + proteinID + " files wrt the current directory. The directory should end in /.")
				directory = str(input("Directory: "))
			elif fileType == "C":
				isPDBFile = False
				print("Enter the directory of the " + proteinID + ".cif file amd compressed " + proteinID + " files wrt the current directory. The directory should end in /.")
				directory = str(input("Directory: "))
			else: 
				print("Invalid file type. Enter P if the structural protein data is formattted as a PDB file or C if the data is formmated as a mmCIF file.")
				fileType = ""
		directory = directory + proteinID
		from PICstatistics import *
		try: 
			PICstatistics(proteinID, directory, isPDBFile)
			print("##########$$$$$$$$$$##########$$$$$$$$$$##########")
		except: 
			print("Error. Make sure you entered all information correctly, all " + proteinID + " files are in the approriate place, and the protein has been compressed and decompressed then try again.")
	elif mode == "R":
		print("Enter the directory of the folder that contains the twenty protein files, each in their own folder. The directory should end in /.")
		directory = str(input("Directory: "))
		try:
			from PICcompression import *
			from PICdecompression import *
			from PICstatistics import *
			import matplotlib.pyplot as plt
			import numpy as np
			import math

			proteinIDs = ["2ja9", "2jan", "2jbp", "2ja8", "2ign", "2jd8", "2ja7", "2fug", "2b9v", "2j28", "6hif", "3j7q", "3j9m", "6gaw", "5t2a", "4ug0", "4v60", "4wro", "6fxc", "4wq1"]
			isPDBFiles = [True for i in range(10)] + [False for i in range(10)]
			LATEXOutputs = []
			numAtoms = []
			txtSizes = []
			gZipSizes = []
			pngSizes = []
			for i in range(len(proteinIDs)):
				proteinDirectory = directory + proteinIDs[i] + "/" + proteinIDs[i]
				compressionStatistics = PICcompress(proteinIDs[i], proteinDirectory, isPDBFiles[i], returnStatistics = True)
				print("##########$$$$$$$$$$##########$$$$$$$$$$##########")
				decompressionStatistics = PICdecompress(proteinIDs[i], proteinDirectory, returnStatistics = True)
				print("##########$$$$$$$$$$##########$$$$$$$$$$##########")
				results = PICstatistics(proteinIDs[i], proteinDirectory, isPDBFiles[i], compressionStatistics = compressionStatistics, decompressionStatistics = decompressionStatistics)
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
		except:
			print("Error. Make sure the twenty protein files are in the approriate places then try again.")
	else:
		print("Invalid option. Select C to compress, D to decompress, S to get advanced statistics of a protein's compression and decompression, or R to reconstruct the results in the article.")
		mode = ""