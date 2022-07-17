# This script loads lists of detected duplication on 5' and 3' flanks and creates a list of proper duplications, that occur on both ends of a pair of insertions.
five=[]
three=[]
# Load 5' duplication list
with open("allreference.5.duplications.txt","r") as f:
    for line in f:
        five.append(line.rstrip("\n").split(","))
f.close()
# Load 3' duplication list
with open("allreference.3.duplications.txt","r") as f:
    for line in f:
        three.append(line.rstrip("\n").split(","))
f.close()
# Compile pairs of insertions that form duplications on 5' and 3' flanks
for i in tuple(five):
    for j in tuple(three):
        if (i[2].split("_")[0]==j[2].split("_")[0] and i[3].split("_")[0]==j[3].split("_")[0]) or (i[2].split("_")[0]==j[3].split("_")[0] and i[3].split("_")[0]==j[2].split("_")[0]):
            a=",".join(i).rstrip(",")
            b=",".join(j).rstrip(",")
            print(a)
            print(b)
            print("")
