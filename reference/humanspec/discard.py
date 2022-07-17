# Mark segmental duplications in the Human vs pan blast results.
def discarddup(filename,duplist):
    dupl=[]
    removed=[]
    with open(duplist,"r") as f:
        for line in f:
            if line=="\n":
                continue
            line=line.rstrip("\n").split("_")
            if len(dupl)==0:
                dupl.append([line[1],line[0]])
            else:
                found=0
                counter=0
                for i in dupl:
                    if i[0]==line[1]:
                        found=1
                        dupl[counter]=dupl[counter]+[line[0]]
                        break
                    counter+=1
                if found==0:
                    dupl.append([line[1],line[0]])
    f.close()
    counter=0
    for i in dupl:
        dupl[counter]=tuple(i)
        counter+=1
    dupl=tuple(dupl)
    found=0
    marker=0
    with open(filename,"r") as f:
        for line in f:
            for d in dupl:
                if line.split("\t")[1].split("_")[0] in set(d[2:]):
                    removed.append(line.split("\t")[1].split("_")[0]+"_"+d[0])
                    found=1
                    break
                if line.split("\t")[1].split("_")[0]==d[1]:
                    print(line.split("\t")[0]+"\t"+"SD"+d[0]+"_"+line.split("\t")[1]+"\t"+"\t".join(line.split("\t")[2:]))
                    marker=1
                    break
            if found==1 or marker==1:
                found=0
                marker=0
                continue
            else:
                print(line,end='')
    f.close()
    return tuple(set(removed))
def main():
    with open("removedp3.txt","w") as reso:
        str1=""
        for i in discarddup("blast.panp3.out",'dups.csv'):
            str1+=i+"\n"
        reso.write(str1)
    reso.close()
if __name__ == "__main__":
    main()
