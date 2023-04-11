# Using readlines()
filein = open('aligned.fasta', 'r')
Lines = filein.readlines()
 
# writing to file
fileout = open('aligned.onlySpike.fasta', 'w')

count = 0
# Strips the newline character


for line in Lines:

    if count%2==0:#das ist ein Header
        fileout.write(line) 
    if count%2==1:#das ist ein Read
        fileout.write(line[21563:25384]) 
        fileout.write("\n")

    count = count+1