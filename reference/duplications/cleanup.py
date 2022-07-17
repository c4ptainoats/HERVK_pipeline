lines=[]
# This script removes duplicated results from BLAST, including same subject - query pairs (A vs A) and inverted results (A-B vs B-A)

# Open blast results
with open("allreference.5.txt",'r') as f:
    for line in f:
        line=tuple(line.rstrip("\n").split("\t"))
        # Disregard A vs A results.
        if line[0]==line[1]:
            pass
        else:
            if len(lines)==0:
                lines.append(line)
            else:
                found=0
                # Disregard B vs A if A vs B is already there.
                for i in lines:
                    if (i[0]==line[1] and i[1]==line[0]):
                        found=1
                        break

                if found==0:
                    lines.append(line)
f.close()
# Print filtered results.
for i in lines:
    print("\t".join(i))
