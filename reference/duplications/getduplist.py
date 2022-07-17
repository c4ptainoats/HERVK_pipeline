pairs=[]
pairnames=[]
flag=0

# This script loads pairs of duplicated sequences and creates a list of insertions that are involved in segmental duplication, grouping them into appropriate groups.

# Open file with listed duplications.
with open('duplications.txt','r') as f:
    for line in f:
        if line=="\n":
            pairs.append(entry)
        else:
            line=",".join(line.split(",")[2:])
            if "\t" in line:
                line=line.replace("\t",",")
            if flag==0:
                entry=[line.rstrip("\n").lstrip("\ufeff")]
                flag=1
            else:
                flag=0
                entry.append(line.rstrip("\n").lstrip("\ufeff"))
f.close()
# Load duplications into memory.

for i in tuple(pairs):
    pairnames.append(tuple([i[0].split(",")[0].split("_")[0],i[0].split(",")[1].split("_")[0]]))

listofdup=[list(pairnames[0])]

# Group insertions and print them.

for i in tuple(pairnames[1:]):
    found=0
    counter=0
    for j in listofdup:
        if i[0] in j and i[1] in j:
            found=1
            break
        if i[0] in j and i[1] not in j:
            listofdup[counter].append(i[1])
            found=1
            break
        if i[1] in j and i[0] not in j:
            listofdup[counter].append(i[0])
            found=1
            break
        counter+=1
    if found==0:
        listofdup.append(list(i))
counter=0
for i in listofdup:
    for j in i[0:-1]:
        print(j,end=',')
        counter+=1
    print(i[-1])
    counter+=1
print(counter)
