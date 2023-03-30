# Using readlines()
filein = open('aligned.onlySpike.fasta', 'r')
lines = filein.readlines()
 
# writing to file
fileout = open('aligned.onlySpike.max20pIndel.noN.fasta', 'w')


# Strips the newline character
i=0
summe = 0
while i<len(lines):
    percN = lines[i+1].count('N')/len(lines[i+1])
    if(percN ==0): #Nehme nur solche Reads, die keine Ns drin haben
        print(percN)
        summe = summe+1

        fileout.write(lines[i])
        fileout.write(lines[i+1])


    i=i+2

print(summe, " von ",len(lines))