# This script selects groups from BLAST matches, that match within the 10bp immidiate flank, to find possible segmental duplications with matching TSDs.

# 5' end selection

"""
with open("allreference.5.grouped.txt","r") as f:
    for line in f:
        good="n"
        line=line.rstrip("\n").split(",")
        if line[-1]=="":
            line=line[:-1]
        counter=0

        for j in line[2:]:
            if counter==0 or counter==2:
                counter+=1
            elif counter==1:
                firstnum=int(j)
                counter+=1
            else:
                secondnum=int(j)
                if firstnum>=490 and secondnum>=490:
                    print(",".join(line))
                counter=0
f.close()
"""
# 3' end selection

with open("allreference.3.grouped.txt","r") as f:
    for line in f:
        good="n"
        line=line.rstrip("\n").split(",")
        if line[-1]=="":
            line=line[:-1]
        counter=0
# extracts the matches and looks for the range
        for j in line[2:]:
            if counter==0:
                firstnum=int(j)
                counter+=1
            elif counter==2:
                secondnum=int(j)
# if the range encompasses first 10bp, select it 
                if firstnum<=10 and secondnum<=10:
                    print(",".join(line))
                counter+=1
            elif counter==1:
                counter+=1
            else:
                counter=0


f.close()
