# Using readlines()

filename = "aligned.onlySpike.max20pIndel.noN.medium1000.fasta"

filein = open(filename, 'r')
lines = filein.readlines()
 
# writing to file
fileout = open(filename+".selDates.fasta", 'w')

i = 0
# Strips the newline character


while i<len(lines):

    if lines[i].count('/')==2:
        fileout.write(lines[i]) 
        fileout.write(lines[i+1]) 
    i = i+2
