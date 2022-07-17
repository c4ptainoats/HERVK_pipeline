# This script opens extracted sequences and creates separate files for 500bp of 5' flanking sequence and 3' flanking sequences.

with open('allhervk.fasta','r') as f:
    for line in f:
        if line.startswith(">"):
            name=line
        else:
            line=line.rstrip("\n")
            with open('allreference.5.fasta','a') as res1:
                res1.write(name)
                res1.write(line[:500]+"\n")
            res1.close()
            with open('allreference.3.fasta','a') as res2:
                res2.write(name)
                res2.write(line[-500:]+"\n")
            res2.close()
f.close()
