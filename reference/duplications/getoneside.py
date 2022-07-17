# This script extracts duplications observed only on 5' or 3' flanks of pairs of insertions, for manual inspection.
dups=[]
refins=[]
ref113=''
# Load reference herv-k113 sequences
with open("LTR113.fasta","r") as a:
    for line in a:
        ref113+=line
a.close()
# Load source file with sequences of extracted viral loci
with open("allhervk.fasta","r") as refi:
    for line in refi:
        if line.startswith(">"):
            insert=[line]
        else:
            insert.append(line)
            refins.append(tuple(insert))
# Load duplication list
with open("duplications.all.txt","r") as f:
    for line in f:
        for i in tuple(line.rstrip("\n").split(",")):
            dups.append(i)
f.close()
dups=set(dups)
pairs=[]
# Load list of duplications on 5' flank
with open("allreference.5.duplications.txt","r") as f:
    start=1
    for line in f:
        if tuple(line.split(","))[2].split("_")[0] in dups and tuple(line.split(","))[3].split("_")[0] in dups:
            pass
        else:
            if start==1:
                start=0
                #print("Single 5' ends:")
            #print(line.rstrip("\n"))
            pairs.append(tuple([line.split(",")[2],line.split(",")[3]]))
f.close()
counter=0
# Save 5' flank only duplications
for i in tuple(pairs):
    found=0
    five=[]
    for j in refins:
        if i[0]==j[0].lstrip(">").rstrip("\n") or i[1]==j[0].lstrip(">").rstrip("\n"):
            five.append(j[0]+j[1])
            found+=1
        if found==2:
            break
    if len(five)>1:
        counter+=1
        with open("5end.duplications"+str(counter)+".fasta",'a') as fiveend:
            fiveend.write(ref113)
            fiveend.write(five[0])
            fiveend.write(five[1])
fiveend.close()
pairs=[]
# Load list of duplications on 3' flank
with open("allreference.3.duplications.txt","r") as f:
    start=1
    for line in f:
#        print(line)
        #print(tuple(line.split(",")))

        if tuple(line.split(","))[2].split("_")[0] in dups and tuple(line.split(","))[3].split("_")[0] in dups:
            pass
        else:
            if start==1:
                start=0
                #print("Single 3' ends:")
            #print(line.rstrip("\n"))
            pairs.append(tuple([line.split(",")[2],line.split(",")[3]]))
counter=0
# Save 3' flank only duplications
for i in tuple(pairs):
    found=0
    three=[]
    for j in refins:
        if i[0]==j[0].lstrip(">").rstrip("\n") or i[1]==j[0].lstrip(">").rstrip("\n"):
            three.append(j[0]+j[1])
            found+=1
        if found==2:
            break
    if len(three)>1:
        counter+=1
        with open("3end.duplications"+str(counter)+".fasta",'a') as threeend:
            threeend.write(ref113)
            threeend.write(three[0])
            try:
                threeend.write(three[1])
            except IndexError:
                print(three)
                break
        threeend.close()
f.close()
