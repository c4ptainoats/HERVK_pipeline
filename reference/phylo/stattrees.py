from os import listdir
import re
# Produces statistic on posterior probabilities in % values for every range, ranges according to resolution set below
resolution=0.05
results=[]
for i in listdir():
    if i.endswith("fixed.tree"):
        results.append([i])
        with open(i,"r") as f:
            digits=[]
            for line in f:
                if line.startswith("tree"):
                    while True:
                        p=re.search("posterior", line)
                        if p==None:
                            break
                        else:
                            e=1
                            dig=''
                            while True:
                                if line[p.span()[1]+e].isdigit()==True or line[p.span()[1]+e]=="." or line[p.span()[1]+e]=="E" or line[p.span()[1]+e]=="-":
                                    dig+=line[p.span()[1]+e]
                                else:
                                    line=line[p.span()[1]+e:]
                                    break
                                e+=1
                            digits.append(float(dig))
            targetval=1
            while targetval-resolution>=-0.1:
                counter=0

                for d in digits:
                    if targetval>=d>targetval-resolution:
                        counter+=1
                results[-1].append(str(round(targetval,2))+'-'+str(round(targetval-resolution,2))+"\t"+str(round((float(counter)/len(digits))*100,2)))
                targetval-=resolution

        f.close()

results=sorted(results,key=lambda a:int(a[0].split(".")[0].split("frag")[1]))
for j in range(len(results[0])-1):
    if j==0:
        for i in results:
            print(i[j],end='\t\t')

    else:
        for i in results:
            print(i[j],end='\t')
    print("")
