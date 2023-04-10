# Using readlines()

filename = "aligned.onlySpike.max20pIndel.noN.supershort100.fasta"

filein = open(filename, 'r')
lines = filein.readlines()
 
# writing to file
fileout = open(filename+"dates.txt", 'w')

count = 0
# Strips the newline character

while count<len(lines):

    line = lines[count][-5:]

    out = lines[count][1:-1]+"\t"+line
    print(out)
    fileout.write(out) 

    count = count+2
