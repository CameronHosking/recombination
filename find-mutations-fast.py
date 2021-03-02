import sys

def main(argv):
    nextLineIsReference = False
    referenceStrain = ""
    for line in open(argv[0]):
        if nextLineIsReference:
            referenceStrain = line
            break
        if line=='>Reference_Sequence\n':
            nextLineIsReference = True

    strain = ""
    for line in open(argv[0]):
        if line[0]=='>':
            strain = line[1:-1]
        else:
            lastMatch = -1
            for i,c in enumerate(line):
                if c == referenceStrain[i]:
                    if(lastMatch!=i-1):
                        start = lastMatch+1
                        end = i
                        print(strain,referenceStrain[start:end],line[start:end],str(start + 1) + ":" + str(end),sep="\t")
                    lastMatch = i
if __name__ == '__main__':
	main(sys.argv[1:])
