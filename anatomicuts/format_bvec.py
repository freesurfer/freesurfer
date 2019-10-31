import sys

if( len(sys.argv)<2):

	print(sys.argv[0] + " input_bvec output_bvec")
else:

	inFile = open(sys.argv[1], 'r')
	outFile = open(sys.argv[2], 'w+')	
	format3N = sys.argv[3] #3N or  N3
	bvecs=[]
	for line in inFile:
		splitted =  line.replace("  "," ").replace("\n","").split(' ')
		bvecs.append(splitted)
		
	if format3N is "3N":
		for i in range(len(bvecs[0])):
			if len(bvecs) >2:
				outFile.write(bvecs[0][i]+" "+bvecs[1][i]+" "+bvecs[2][i]+ "\n")
			else:
				outFile.write(bvecs[0][i]+"\n")
	else:
		for j in range(len(bvecs[0])):
			for i in range(len(bvecs)):
				outFile.write(bvecs[i][j]+" ")
			outFile.write("\n")
	


	outFile.close()
	inFile.close()

