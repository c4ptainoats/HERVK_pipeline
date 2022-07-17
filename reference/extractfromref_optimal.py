# This script extracts ranges from reference genome files and prints the sequences on screen. Instructions for usage below.

# Dictionary for sequences in reverse order.
reverseatcg={"A":"T",
"T":"A",
"C":"G",
"G":"C",
"N":"N",
'g':'c',
'c':'g',
'a':'t',
't':'a',
'n':'n'
}
# Dictionary for chromosome names in the official files.
# There are differences in the names of the chromosomes between build 37 and 38 of the human genome. Tbe following lines of code standardise the chromosome name to just the numerical identifier of the chromosome (i.e. gi|224589800|ref|NC_000001.10| (build 37) and chr1 (build 38) are now both recognised by "1")
chrname={
"gi|224589800|ref|NC_000001.10|":"1",
"gi|224589811|ref|NC_000002.11|":"2",
"gi|224589815|ref|NC_000003.11|":"3",
"gi|224589816|ref|NC_000004.11|":"4",
"gi|224589817|ref|NC_000005.9|":"5",
"gi|224589818|ref|NC_000006.11|":"6",
"gi|224589819|ref|NC_000007.13|":"7",
"gi|224589820|ref|NC_000008.10|":"8",
"gi|224589821|ref|NC_000009.11|":"9",
"gi|224589801|ref|NC_000010.10|":"10",
"gi|224589802|ref|NC_000011.9|":"11",
"gi|224589803|ref|NC_000012.11|":"12",
"gi|224589804|ref|NC_000013.10|":"13",
"gi|224589805|ref|NC_000014.8|":"14",
"gi|224589806|ref|NC_000015.9|":"15",
"gi|224589807|ref|NC_000016.9|":"16",
"gi|224589808|ref|NC_000017.10|":"17",
"gi|224589809|ref|NC_000018.9|":"18",
"gi|224589810|ref|NC_000019.9|":"19",
"gi|224589812|ref|NC_000020.10|":"20",
"gi|224589813|ref|NC_000021.8|":"21",
"gi|224589814|ref|NC_000022.10|":"22",
"gi|224589822|ref|NC_000023.10|":"X",
"gi|224589823|ref|NC_000024.9|":"Y",
"chr1":"1",
"chr2":"2",
"chr3":"3",
"chr4":"4",
"chr5":"5",
"chr6":"6",
"chr7":"7",
"chr8":"8",
"chr9":"9",
"chr10":"10",
"chr11":"11",
"chr12":"12",
"chr13":"13",
"chr14":"14",
"chr15":"15",
"chr16":"16",
"chr17":"17",
"chr18":"18",
"chr19":"19",
"chr20":"20",
"chr21":"21",
"chr22":"22",
"chrX":"X",
"chrY":"Y",
"1":"1",
"2":"2",
"3":"3",
"4":"4",
"5":"5",
"6":"6",
"7":"7",
"8":"8",
"9":"9",
"10":"10",
"11":"11",
"12":"12",
"13":"13",
"14":"14",
"15":"15",
"16":"16",
"17":"17",
"18":"18",
"19":"19",
"20":"20",
"21":"21",
"22":"22",
"X":"X",
"Y":"Y",
"chr1_KI270765v1_alt":"chr1_KI270765v1_alt",
"chr2_KI270768v1_alt":"chr2_KI270768v1_alt",
"chr3_KI270779v1_alt":"chr3_KI270779v1_alt",
"chr3_KI270780v1_alt":"chr3_KI270780v1_alt",
"chr3_KI270895v1_alt":"chr3_KI270895v1_alt",
"chr3_KI270924v1_alt":"chr3_KI270924v1_alt",
"chr3_KI270934v1_alt":"chr3_KI270934v1_alt",
"chr3_KI270935v1_alt":"chr3_KI270935v1_alt",
"chr3_KI270936v1_alt":"chr3_KI270936v1_alt",
"chr3_KI270937v1_alt":"chr3_KI270937v1_alt",
"chrUn_KI270749v1":"chrUn_KI270749v1",
"chrUn_gl000219":"chrUn_gl000219",
"chrUn_gl000222":"chrUn_gl000222",
"chrUn_gl000223":"chrUn_gl000223",
"chrUn_gl000231":"chrUn_gl000231",
"chrUn_gl000212":"chrUn_gl000212",
"chrUn_gl000219":"chrUn_gl000219",
"chrUn_gl000232":"chrUn_gl000232",
"chr17_ctg5_hap1":"chr17_ctg5_hap1",
"chrUn_GL000219v1":"chrUn_GL000219v1"
}



# List for extracted locations.
# The script requires a list of the locations of the insertions, as a text file in the format:
# [cytogenetic identifier]  [chromosome(numeric)] [orientation] [5' position in the genome] [3' position in the genome] <- these fields are tab seperated
# The output format combines the above fields to create a new name for each extracted sequence, followed ny the DNA sequence extracted ***
locs=[]

# The following section opens TSV file with viral location to extract and loads it into memory.
with open('toextrac.left.txt','r') as f:
    for line in f:
        locs.append(tuple(line.rstrip("\n").split("\t")))
f.close()
locs=tuple(locs)
# Create variables for sequences and headers.
genome=[]
header=''
sequence=''


# This section opens fasta reference genome sequence and loads individual chromosomes into memory, also scanning the list of known HERV-K insertions and extracting sequences of insertions from each analysed chromosome.
with open('hg38.fa','r') as g:
    for line in g:
# If line starts with header sign, get the chromosome name.
        if line.startswith(">"):
            if header=='':
                analysed=[]
                header=chrname[line.rstrip("\n").lstrip(">")]
# Pick insertions that belong to the analysed chromosome
                for i in locs:

                    if header in set(i) or ("chr"+str(header)) in set(i):
                        analysed.append(i)

            else:
                analysed=tuple(analysed)
# Extract insertions from chromosome loaded to the memory before.
                for i in analysed:
                    start=int(i[3])
                    end=int(i[4])
                    print(">"+i[0]+"_chr"+chrname[i[1]]+":"+str(i[3])+"-"+str(i[4]))
# Print sequences on screen, change orientation if necessary.
                    if i[2]=="+":
                        print(sequence[start-1:end])
                    else:
                        seq=sequence[start-1:end][::-1]
                        seq1=''
                        for char in seq:
                            seq1+=reverseatcg[char]
                        print(seq1)
                analysed=[]
                sequence=''
# If line starts with header sign, get the chromosome name. Pass on chromosomes that you don't know according to dictionary at the beginning.
                try:
                    header=chrname[line.rstrip("\n").lstrip(">")]
                except KeyError:
                    continue
# Pick insertions that belong to the analysed chromosome.
                for i in locs:
                    if header in set(i) or ("chr"+str(header)) in set(i):
                        analysed.append(i)
# If line doesn't have a header sign, assemble it into sequence.
        else:
            sequence+=line.rstrip("\n")


g.close()
analysed=tuple(analysed)
# Extract insertions from chromosome loaded to the memory before.
for i in analysed:
    print(i)
    start=int(i[3])
    end=int(i[4])
    print(">"+i[0]+"_chr"+chrname[i[1]]+":"+str(i[3])+"-"+str(i[4]))
# Print sequences on screen, change orientation if necessary
    if i[2]=="+":
        print(sequence[start-1:end])
    else:
        seq=sequence[start-1:end][::-1]
        seq1=''
        for char in seq:
            seq1+=reverseatcg[char]
        print(seq1)
