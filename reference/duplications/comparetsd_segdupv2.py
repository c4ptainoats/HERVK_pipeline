
lines=[]
tsds=[]
# This script loads TSD infromation from the table and grouped results to produce possible pairs of duplications that share the same TSD and produce positive BLAST results.

# Load TSD status list
with open("tsd_wooble_final.csv",'r') as f:
    for line in f:
        line=tuple(line.split(","))
        if "TRUNCATED" in line[18].rstrip("\n") and "3" in line[18].rstrip("\n"):
            tsds.append(tuple([line[0],line[16].rstrip("\n")]))
        elif "TRUNCATED" in line[18].rstrip("\n") and "5" in line[18].rstrip("\n"):
            tsds.append(tuple([line[0],line[17].rstrip("\n")]))
        elif "-" in line[18].rstrip("\n"):
            tsds.append(tuple([line[0],(line[16].rstrip("\n")+"/"+line[17].rstrip("\n"))]))
        else:
            tsds.append(tuple([line[0],line[18].rstrip("\n")]))
f.close()
# Open 3' ends selected for duplications
with open("allreference.3.selected.txt",'r') as f:
    for line in f:
        line=line.lstrip("\ufeff")
        tsd1=""
        tsd2=""
        liner=line.rstrip("\n").split(",")
        liner[0]=liner[0].split("_")[0]
        liner[1]=liner[1].split("_")[0]
        liner=tuple(liner)
# Search for TSDs in the status list
        for i in tsds:
            if i[0]==liner[0]:
                tsd1=i[1].rstrip("\n")
            if i[0]==liner[1]:
                tsd2=i[1].rstrip("\n")


        if ("/" in tsd1 and "/" not in tsd2 and (tsd2==tsd1.split("/")[0] or tsd2==tsd1.split("/")[1])) or ("/" in tsd2 and "/" not in tsd1 and (tsd1==tsd2.split("/")[0] or tsd1==tsd2.split("/")[1])):
            print(tsd1+","+tsd2+","+line,end='')
        else:
            if tsd1==tsd2:
                print(tsd1+","+tsd2+","+line,end='')
            else:
#### FOR 3' ENDS
## compare BLAST match, select if a pair on 3' end matches at worst from 2nd basepair of TSD/flank
                lines=line.rstrip("\n").split(",")
                if lines[-1]=="":
                    lines=lines[:-1]
                counter=0

                for j in lines[2:]:
                    if counter==0:
                        firstnum=int(j)
                        counter+=1
                    elif counter==2:
                        secondnum=int(j)
                        if firstnum<=2 and secondnum<=2:
                            print(tsd1+","+tsd2+","+line,end='')
                        counter+=1
                    elif counter==1:
                        counter+=1
                    else:
                        counter=0


### FOR 5' ENDS
## compare BLAST match, select if a pair on 3' end matches at worst from 2nd basepair of TSD/flank
"""
                lines=line.rstrip("\n").split(",")
                if lines[-1]=="":
                    line=line[:-1]
                counter=0

                for j in lines[2:]:
                    if counter==0 or counter==2:
                        counter+=1
                    elif counter==1:
                        firstnum=int(j)
                        counter+=1
                    else:
                        secondnum=int(j)
                        if firstnum>=499 and secondnum>=499:
                            print(tsd1+","+tsd2+","+line,end='')
                        counter=0
"""
f.close()
