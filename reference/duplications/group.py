# This script groups all of the fragmented results between subject-query pairs and prints the grouped results.

allleft=[]

with open("allreference.5.clean.txt","r") as f:
    for line in f:
        line=line.rstrip("\n").split("\t")
        if len(allleft)==0:
            allleft.append(line[0:2]+line[6:10])
        else:
            found=0
            counter=0
            for i in allleft:
                if line[0]==i[0] and line[1]==i[1]:
                    allleft[counter].append(line[6])
                    allleft[counter].append(line[7])
                    allleft[counter].append(line[8])
                    allleft[counter].append(line[9])
                    found=1
                counter+=1
            if found==0:
                allleft.append(line[0:2]+line[6:10])

f.close()
for i in allleft:
    print(i[0]+","+i[1],end=',')
    numbers=i[2:]
    a=0
    for j in numbers:
        print(j,end=",")

    print("")
