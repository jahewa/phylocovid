# Using readlines()
filename = "aligned.fasta"

filein = open(filename, 'r')
lines = filein.readlines()
 
# writing to file
fileout = open(filename + "max20pIndel_noN_withDate.fasta", 'w')

# Strips the newline character
i=0
summe = 0
while i<len(lines):
    percInDel = lines[i+1].count('-')/len(lines[i+1])
    percN = lines[i+1].count('N')/len(lines[i+1])

    date = lines[i][-5:-1]
    slashCount = lines[i].count('/')

    if(percInDel<=0.2 and percN ==0): #Filter for percInDel <=0.2 and no N; and format fits
        if(slashCount==2):
        #if(date == "2019" or date == "2020" or date == "2021" or date == "2022" or date=="2023"):
            print(i, percInDel, percN)
            summe = summe+1
            fileout.write(lines[i])
            fileout.write(lines[i+1])
    i=i+2

print(summe, " von ",len(lines))