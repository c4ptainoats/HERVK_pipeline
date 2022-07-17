# This script scans all the insertions and looks for single basepair differences between flanking TSDs. If it finds one, it introduces it into the table for future use. The single basepair difference between flanking TSDs is allowed due to random mutation rate and the fact that the polymerase that introduces TSDs has no repair mechanism which might make the mutations more prone to occur in TSD zone.

# Open file with TSD result.
with open('solotsdv31.csv','r') as f:
    for line in f:
        line1=line.rstrip("\n")
        line=line.rstrip("\n").rstrip(",").split(",")
# If there is no matching TSD or insertion is truncated print the line.
        if line[-1]=="-" or "truncated" in line[-1]:
            print(line1)
        else:
# If there is just a matching TSD print the line.
            if (line[-1] in line[-3][-len(line[-1]):]) and (line[-1] in line[-2][:len(line[-1])]):
                print(line1)
            else:
# If the TSD matches with a mutation, print both variants.
                if line[-1] in line[-3][-len(line[-1]):]:
                    print(",".join(line)+"/"+line[-2][:len(line[-1])])
                else:
                    print(",".join(line)+"/"+line[-3][-len(line[-1]):])
                    
f.close()