# This script analyses TSD data of all the insertions and outputs the numbers of insertions involved in duplication and recombination events. It compares TSDs and flanking sequences of insertions to create lists of inserts that have the same TSD and tries to categorize if the TSDs come from recombinations or segmental duplications. It requires output of the 'aligner.py' script to get the duplication data.


# Lists for insertions that are not truncated but have different flanks, 5' truncated ones, 3' truncated ones and normal inserts with proper TSDs.
recombinators=[]
five=[]
three=[]
normal=[]
fiventhree=[]
fiven=0
threen=0
fntn=0
notsdn=0
tsdn=0

# Added after the report - in results table there were inserts with non-matching TSDs (8) and 3' (8) or 5' (6) truncated that i reported as no activity. How can this be? The function below takes the 8bp flank in each of these cases, and systematically introduces 1 mutation at each position and compares this with the established list of TSDs (provided in table MikeDBv43.grch.csv, big table from report - (back in the report had mixed up locations - MikeDBv41 - in google drive)).
def singlemut(seq):
    if type(seq)==str:
        nucl=tuple(["A","T","C","G"])
        seq=seq.upper()
        result=[]
        count=0
        for i in seq:
            for j in nucl:
                seq2=list(seq)
                seq2[count]=j
                result.append("".join(seq2))
            count+=1
        result=tuple(set(result))
    elif type(seq)==list or type(seq)==tuple:
            nucl=tuple(["A","T","C","G"])
            result=[]
            for z in seq:
                count=0
                z=z.upper()
                for i in z:
                    for j in nucl:
                        seq2=list(z)
                        seq2[count]=j
                        result.append("".join(seq2))
                    count+=1
            result=tuple(set(result))
    else:
        result="Unknown input type, please provide a flank as a string or a list of potential tsds"
    return result


# Open data source with insertion information, along with TSDs, in format : name, presence in HG19 genome, presence in HG38 geneome, chromosome, orientation, hg19 location, length, famili, age, age ,range, comment, hg38 location, 5' end, 3' end, 5' flank, 3' flank, TSD.
with open('tsd_wooble_final.csv','r') as f:
    for line in f:
        lines=line.rstrip("\n").lstrip("\ufeff").split(",")

# Get insertion name
        ident=tuple([lines[0],line.rstrip('\n')])
# Get left and right flanks
        ltsd0=tuple([lines[16][-4:],lines[16][-5:],lines[16][-6:],lines[16][-7:],lines[16][-8:]])
        rtsd0=tuple([lines[17][:4],lines[17][:5],lines[17][:6],lines[17][:7],lines[17][:8]])
        ltsd=singlemut(tuple([lines[16][-4:],lines[16][-5:],lines[16][-6:],lines[16][-7:],lines[16][-8:]]))
        rtsd=singlemut(tuple([lines[17][:4],lines[17][:5],lines[17][:6],lines[17][:7],lines[17][:8]]))
# If no proper TSD, add to recombinators.
        if lines[18]=="-":
            recombinators.append(tuple([ident,ltsd0,rtsd0,ltsd,rtsd]))
            notsdn+=1
# If 5' in TSD section, add to 5' truncated
        elif "5'" in lines[18] and "3'" not in lines[18]:
            five.append(tuple([ident,"-",rtsd0,"-",rtsd]))
            fiven+=1
# If 3' in TSD section, add to 3' truncated
        elif "3'" in lines[18] and "5'" not in lines[18]:
            threen+=1
            three.append(tuple([ident,ltsd0,"-",ltsd,"-"]))
# If both 5' and 3' in TSD section, add to normal ones
        elif "5'" not in lines[18] and "3'" not in lines[18] and lines[18]!="-":
            tsdn+=1
            normal.append(tuple([ident,lines[18]]))
        elif "5'" in lines[18] and "3'" in lines[18]:
            fntn+=1
            fiventhree.append(tuple([ident,"-"]))
normaltsds=[]


# Function that outputs intersection of two lists.
def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))
# Load real duplications from a separate file.
realduplication=[]
dupins=[]
with open('duplications.all.good.txt','r') as f:
    for line in f:
        realduplication.append(tuple(line.rstrip("\n").split(",")))
f.close()
realduplication=tuple(realduplication)
for i in realduplication:
    for j in i:
        dupins.append(j)
dupins=set(dupins)
dupfive=0
dupthree=0
dupfnt=0
duptsd=0
dupnotsd=0
# Count insertions in different categories.
#print(dupins)
for i in five:
    if i[0][0] in dupins:
        dupfive+=1
#        five.remove(i)
for i in tuple(dupins):
    for j in five:
        if j[0][0]==i:
            five.remove(j)
            break

for i in three:
    if i[0][0] in dupins:
        dupthree+=1
#        three.remove(i)
for i in tuple(dupins):
    for j in three:
        if j[0][0]==i:
            three.remove(j)
            break

for i in normal:
    if i[0][0] in dupins:
        duptsd+=1
#        normal.remove(i)
for i in tuple(dupins):
    for j in normal:
        if j[0][0]==i:
            normal.remove(j)
            break


for i in recombinators:
    if i[0][0] in dupins:
        dupnotsd+=1
#        recombinators.remove(i)
for i in tuple(dupins):
    for j in recombinators:
        if j[0][0]==i:
            recombinators.remove(j)
            break
for i in fiventhree:
    if i[0][0] in dupins:
        dupfnt+=1

for i in tuple(dupins):
    for j in fiventhree:
        if j[0][0]==i:
            fiventhree.remove(j)
            break


# Analyse all insertions with good TSDs recognisable on both ends. The script compares the TSDs against each other and groups them according to TSD recognised. It also creates a list of all TSDs from the good insertions to further compare the truncated ones with them and try to predict if any of the truncated ones might have been created by recombination/duplication of a good one. If a mutation occurs, both of the TSD variants are taken into account.

# Add all TSDs from normal inserts to a list, create a set.
#for i in normal:
#    print(i)
for i in normal:
    normaltsds.append(i[1])
normaltsds=list(set(normaltsds))
# Comparing normal to normal inserts (Proper TSDs)

# Put every TSD into a list.
for i in normaltsds:
    normaltsds[normaltsds.index(i)]=[normaltsds[normaltsds.index(i)]]

# Scan every normal insertion and put it into groups that have their TSDs.
for i in range(0,len(normaltsds)):
    for j in normal:
        if j[-1]==normaltsds[i][0]:
            normaltsds[i].append(j)
count=0

normaltonormal=[]


# Count the number of insertion involved in recombination.

for i in normaltsds:
    if len(i)>2:
        normaltonormal.append([])
        j=0

        for j in i[1:]:

            normaltonormal[-1].append(j[0][1].split(",")[0]+","+j[0][1].split(",")[3]+","+j[0][1].split(",")[5]+","+j[0][1].split(",")[6]+","+j[0][1].split(",")[7]+","+j[0][1].split(",")[4]+","+j[0][1].split(",")[3]+","+j[0][1].split(",")[18])



# Split results if mutation happened to account both TSD variants.
normaltsds=[]
for i in normal:
    if "/" in i[-1]:
        normaltsds.append(i[-1].split("/")[0])
        normaltsds.append(i[-1].split("/")[1])
    else:
        normaltsds.append(i[-1])
print("good to good")
for i in normaltonormal:
    print(i[0]+","+i[1])
    if len(i)>2:
        for j in i[2:]:
            print(i[0]+","+j)



print("--------------------")

# This section analyses all recombined insertions - ones that are not truncated but the TSDs don't match. The recognised flanking sequences are extended from 4-8bp and compared against each other, inserts with recognisable TSDs and 3' as well as 5' truncated insertions.

#print("Recombined")
retsds=[]
# Append all TSDs of recombined insertions into a list.
for i in recombinators:
    retsds.append((i[-4],i[-3],i[-2],i[-1]))
# Create variables for 5' to 3' and 3' to 5' recombination, 5' to 5' or 3' to 3' or (5' to 5' AND 3' to 3')
proper=0
half=0
segrecomdup=0
alreadyfound=[]
recomtorecomres=[]
halfhalf=[]
# Compare recombined insertion to other recombined insertions.

# Checking truncated insertions against each other for TSDs. Using flanks extending from 4 to 8bp and cross checking 5' to 3' and 3' to 5' as well as 5' to 5' and 3' to 3' ends.

# Go backwards on TSDs
for i in range(len(retsds) - 1, -1, -1):
    # Get a list of inserts minus the checked ones
    temp=retsds[:i]
    # Remove the ones that already recombined
    if len(alreadyfound)>0:
        for e in alreadyfound:
            if e in temp:
                temp.remove(e)
    # Go through rest inserts.
    for j in temp:
        # Check both extending flanks from 4-8bp for overlaps. Put the ones you find in a list.
        if any(elem in j[0]  for elem in retsds[i][3]) and any(elem in j[1]  for elem in retsds[i][2]):
            proper+=1
            alreadyfound.append(j)
            break
        # If just one extending flank appears in other insertion, put it in half recombined list.
        if any(elem in j[0]  for elem in retsds[i][3]) or any(elem in j[1]  for elem in retsds[i][2]):
            half+=1
            alreadyfound.append(j)
            halfhalf.append(j)
            break

# Check if 4-8bp from both flanks appear in other insertion in the same fashion.
for i in range(len(retsds) - 1, -1, -1):
    temp=retsds[:i]
    for j in temp:

        if any(elem in j[0]  for elem in retsds[i][2]) and any(elem in j[1]  for elem in retsds[i][3]):
            segrecomdup+=1
            break
recombpairs=[]

# Go through all the recombinating insertions again again. Create paris of insertions recombining 5' to 3' end and 3' to 5' end to find properly recombined insertions.
for i in range(len(recombinators)-1,-1,-1):
# Get a list of inserts minus the checked ones
    temp=recombinators[:i]
# Remove ones that have been found in previous iterations.
    if len(recombpairs)>0:
        for e in recombpairs:
            if e in temp:
                temp.remove(e)
# Scan inserts for 5' to 5' and 3' to 3' matches.
    for j in temp:

        if any(elem in j[4]  for elem in recombinators[i][2]):
            recombpairs.append([j[0][1],recombinators[i][0][1],max(intersection(j[4],recombinators[i][2]), key=len)])
        if any(elem in j[3]  for elem in recombinators[i][1]):
            recombpairs.append([j[0][1],recombinators[i][0][1],max(intersection(j[3],recombinators[i][1]), key=len)])
recombpairsres=[]
# Add results to list.

for i in recombpairs:
    recombpairsres.append(i[0].split(",")[0])
    recombpairsres.append(i[1].split(",")[0])

print("recom to recom")
for i in recombpairs:
    print(i[0].split(",")[0]+","+i[0].split(",")[3]+","+i[0].split(",")[5]+","+i[0].split(",")[6]+","+i[0].split(",")[7]+","+i[0].split(",")[4]+","+i[-1]+","+i[1].split(",")[0]+","+i[1].split(",")[3]+","+i[1].split(",")[5]+","+i[0].split(",")[6]+","+i[1].split(",")[7]+","+i[1].split(",")[4]+","+i[2]+",")
# The flanks of recombined insertions are compared to the list of TSDs recognised in the proper insertions to find if any of them have inserted due to duplication of another, existing insertion prior to recombination, or may be a result of a different event involving a known, proper insertion. The process involves scanning for a single mutation if it was observed to occur and outputting the longest match from the selected ones, since the 4-8bp window is compared from the recombined side.

recomfull=[]
recomtofl=[]
for i in range(len(retsds) - 1, -1, -1):
    longest=[]
    for j in normal:
# If there is mutation, compare each of the variants to the recombined TSDs.
        if "/" in j[-1]:
            if any(elem==j[-1].split("/")[0]  for elem in retsds[i][2]):# or any(elem==j[-1].split("/")[0]  for elem in retsds[i][3]):
                longest.append(j[-1].split("/")[0])
            if any(elem==j[-1].split("/")[1]  for elem in retsds[i][3]):# or any(elem==j[-1].split("/")[1]  for elem in retsds[i][2]) :
                longest.append(j[-1].split("/")[1])
# If there's no wooble, just do the comparison to the TSD.
        else:
            if any(elem==j[-1]  for elem in retsds[i][2]) or any(elem==j[-1]  for elem in retsds[i][3]):
                longest.append(j[-1])
# Get the longest match i.e. found ATCT in HERV1 vs HERV2, but found ATCTGT in HERV1 vs HERV3, so the pair HERV1-HERV3 is most likely correct.

    if len(longest)!=0:
        recomfull.append(max(longest, key=len))
        for nor in normal:
# Compare both TSDs if you see wooble. Assemble them in pairs. Put longest TSDs you find for a possible pair.
            if "/" in nor[-1]:
                if nor[-1].split("/")[0]==max(longest, key=len) or nor[-1].split("/")[1]==max(longest, key=len):
                    pair=nor[0][1].split(",")[0]+","+nor[0][1].split(",")[3]+","+nor[0][1].split(",")[5]+","+nor[0][1].split(",")[6]+","+nor[0][1].split(",")[4]+","+nor[0][1].split(",")[-1]
            else:
# Do the same for single, good TSDs. Assemble found inserts in paris.
                if nor[-1]==max(longest, key=len):

                    pair=nor[0][1].split(",")[0]+","+nor[0][1].split(",")[3]+","+nor[0][1].split(",")[5]+","+nor[0][1].split(",")[6]+","+nor[0][1].split(",")[4]+","+nor[0][1].split(",")[-1]
# Get the names of the insertions for each found TSDs to create complete pairs with all data on recombining inserts
        for z in recombinators:
            if z[-2]==retsds[i][-2] and z[-1]==retsds[i][-1]:
                recomtofl.append([z[0][1].split(",")[0]+","+z[0][1].split(",")[3]+","+z[0][1].split(",")[5]+","+z[0][1].split(",")[6]+","+z[0][1].split(",")[4]+","+max(longest, key=len),pair])

# 5'

# The 5' truncated insertions are compared to all of the other insertions in a similar fashion as the recombined inserts before. The 4-8bp possible matches following the 3' end of the insertion are taken into account, as well as the possibility of mutation. 5' inserts are checked against each ofther for duplications, 3' truncated inserts for recombinations as well as agains insertions with proper TSDs and ones with recognisable ends but non-matching TSDs.

# Create some lists for results.
fivefullres=[]
fivefull=[]

# Go through 5' truncated inserts, backwards.
for i in range(len(five) - 1, -1, -1):
    longest=[]
# Compare them with all known good TSDs.
    for j in normal:
# Use variants if mutation happens.
        if "/" in j[-1]:
            if any(elem==j[-1].split("/")[1]  for elem in five[i][-1]):
                longest.append(j[-1].split("/")[1])
        else:
            if any(elem==j[-1]  for elem in five[i][-1]):
                longest.append(j[-1])
# Choose the longest overlapping TSD combo.
    if len(longest)!=0:
        fivefull.append(max(longest, key=len))
# Get name of the insertion along with other data and put it into a pair. Account wooble.
        for zg in normal:
            if "/" in zg[-1] and (zg[-1].split("/")[1]==max(longest, key=len)):
                fivefullres.append([five[i][0][1].split(",")[0]+","+five[i][0][1].split(",")[3]+","+five[i][0][1].split(",")[5]+","+five[i][0][1].split(",")[6]+","+five[i][0][1].split(",")[4]+","+max(longest, key=len),zg[0][1].split(",")[0]+","+zg[0][1].split(",")[3]+","+zg[0][1].split(",")[5]+","+zg[0][1].split(",")[6]+","+zg[0][1].split(",")[4]+","+max(longest, key=len)])
# Get name of the insertion along with other data and put it into a pair.
            else:
                if zg[-1]==max(longest, key=len):
                    fivefullres.append([five[i][0][1].split(",")[0]+","+five[i][0][1].split(",")[3]+","+five[i][0][1].split(",")[5]+","+five[i][0][1].split(",")[6]+","+five[i][0][1].split(",")[4]+","+max(longest, key=len),zg[0][1].split(",")[0]+","+zg[0][1].split(",")[3]+","+zg[0][1].split(",")[5]+","+zg[0][1].split(",")[6]+","+zg[0][1].split(",")[4]+","+max(longest, key=len)])



# Compare 5' truncated insertions to 3' truncated insertions. 4-8bp possible TSDs are taken from each of the insertion and longest overlapping TSDs are selected.
fivethreeres=[]
fivethree=[]

# 5' insertions are compared to the recombined insertions (insertions with good ends but non-matching TSDs) below. Both of the ends of recombined insertion are compared against 3' end of the insert to find 4-8bp matches.

fiverecom=[]

# Go through all 5' trucated insertions, backwards.
for i in range(len(five) - 1, -1, -1):
    longest=[]
    counter=0
# Go through all recombined insertion TSDs.
    for j in retsds:
# Find longest matches on both ends
        if any(elem in j[-1]  for elem in five[i][-3]):
            longest+=intersection(j[-1], five[i][-3])

        counter+=1

# Choose the lontest match.
    if len(longest)!=0:
        for re in recombinators:
            if max(longest, key=len) in re[-1]:
                break
# Add the results to list with all the details for each pair of insertion.
        fiverecom.append([five[i][0][1].split(",")[0]+","+five[i][0][1].split(",")[3]+","+five[i][0][1].split(",")[5]+","+five[i][0][1].split(",")[6]+","+five[i][0][1].split(",")[4]+","+max(longest, key=len),re[0][1].split(",")[0]+","+re[0][1].split(",")[3]+","+re[0][1].split(",")[12]+","+re[0][1].split(",")[13]+","+re[0][1].split(",")[4]+","+max(longest, key=len)])


# The 3' truncated insertions are compared to all of the other insertions in a similar fashion as the recombined inserts before. The 4-8bp possible matches following the 5' end of the insertion are taken into account, as well as the possibility of mutation. 3' inserts are checked against each ofther for duplications, 5' truncated inserts for recombinations as well as agains insertions with proper TSDs and ones with recognisable ends but non-matching TSDs.

threefullres=[]
threefull=[]
# Go through 3' truncated insertions, backwards.
for i in range(len(three) - 1, -1, -1):
    longest=[]
# Go through normal ones.
    for j in normal:
# Get both TSD variants if wooble happens and compare with 3' truncated flank 4-8bp.
        if "/" in j[-1]:
            if any(elem==j[-1].split("/")[0]  for elem in three[i][-4]):

                longest.append(j[-1].split("/")[0])
# In case of no wooble, use the TSD and compare with 3' truncated flank 4-8bp.
        else:
            if any(elem==j[-1]  for elem in three[i][-4]):

                longest.append(j[-1])
# Get longest variant of recombination found.
    if len(longest)!=0:

        threefull.append(max(longest, key=len))
# Find details about each insertions, put pairs to a list.
        for zg in normal:
            if zg[-1]==max(longest, key=len):
                threefullres.append([three[i][0][1].split(",")[0]+","+three[i][0][1].split(",")[3]+","+three[i][0][1].split(",")[5]+","+three[i][0][1].split(",")[6]+","+three[i][0][1].split(",")[4]+","+max(longest, key=len),zg[0][1].split(",")[0]+","+zg[0][1].split(",")[3]+","+zg[0][1].split(",")[5]+","+zg[0][1].split(",")[6]+","+zg[0][1].split(",")[4]+","+max(longest, key=len)])

# Compare 3' truncated insertions to recombined insertions.
threerecom=[]
for i in range(len(three) - 1, -1, -1):
    longest=[]
    counter=0
# Get intersections of possible TSD combinations and choose the longest pair.
    for j in retsds:

        if any(elem in j[-2]  for elem in three[i][-4]):
            longest+=intersection(j[-2], three[i][-4])

        counter+=1
# Get the longest pair for recombinations. Append pairs to a list.
    if len(longest)!=0:
        for re in recombinators:
            if max(longest, key=len) in re[-1] or max(longest, key=len) in re[-2]:
                break
        threerecom.append([three[i][0][1].split(",")[0]+","+three[i][0][1].split(",")[3]+","+three[i][0][1].split(",")[5]+","+three[i][0][1].split(",")[6]+","+three[i][0][1].split(",")[4]+","+max(longest, key=len),re[0][1].split(",")[0]+","+re[0][1].split(",")[3]+","+re[0][1].split(",")[12]+","+re[0][1].split(",")[13]+","+re[0][1].split(",")[4]+","+max(longest, key=len)])

# Compare 3' truncated inserts against each other to find possible duplications:

threetothree=[]
for i in range(len(three) - 1, -1, -1):
    temp=three[:i]
    for j in temp:
# Accept results only if the TSD matches completely - possible duplication
        if any(elem in three[i][-2]  for elem in j[-4]):
            tsdl=[]
            for t in three[i][-2]:
                if t in j[-4]:
                    tsdl.append(t)
# Append results to a list.
            threetothree.append([three[i][0][1].split(",")[0]+","+three[i][0][1].split(",")[3]+","+three[i][0][1].split(",")[11]+","+three[i][0][1].split(",")[12]+","+three[i][0][1].split(",")[4]+","+max(tsdl,key=len),j[0][1].split(",")[0]+","+j[0][1].split(",")[3]+","+j[0][1].split(",")[11]+","+j[0][1].split(",")[12]+","+j[0][1].split(",")[4]+","+max(tsdl,key=len)])


# Compare 5' truncated inserts against each other to find possible duplications:

fivetofive=[]
for i in range(len(five) - 1, -1, -1):
    temp=five[:i]
    for j in temp:
# Accept results only if the TSD matches completely - possible duplication
        if any(elem in five[i][-1]  for elem in j[-3]):
            tsdl=[]
            for t in five[i][-1]:
                if t in j[-3]:
                    tsdl.append(t)
# Append results to a list.

            fivetofive.append([five[i][0][1].split(",")[0]+","+five[i][0][1].split(",")[3]+","+five[i][0][1].split(",")[11]+","+five[i][0][1].split(",")[12]+","+five[i][0][1].split(",")[4]+","+max(tsdl,key=len),j[0][1].split(",")[0]+","+j[0][1].split(",")[3]+","+j[0][1].split(",")[11]+","+j[0][1].split(",")[12]+","+j[0][1].split(",")[4]+","+max(tsdl,key=len)])



fivenames=[]
threenames=[]
allnames=[]
all5names=[]
all3names=[]


print("Recombined to FL")
for i in recomtofl:
    recombpairsres.append(i[0].split(",")[0])
    allnames.append(i[1].split(",")[0])
    for j in i[0:-1]:
        print(j,end=',')
    print(i[-1])
print("Recombined to 5'tr")
for i in fiverecom:
    fivenames.append(i[0].split(",")[0])
    recombpairsres.append(i[1].split(",")[0])
    for j in i[0:-1]:
        print(j,end=',')
    print(i[-1])
print("Recombined to 3'tr")
for i in threerecom:
    threenames.append(i[0].split(",")[0])
    recombpairsres.append(i[1].split(",")[0])
    for j in i[0:-1]:
        print(j,end=',')
    print(i[-1])
print("5'tr to 5'tr")
for i in fivetofive:
    fivenames.append(i[0].split(",")[0])
    fivenames.append(i[1].split(",")[0])
    for j in i[0:-1]:
       print(j,end=',')
    print(i[-1])
print("3'tr to 3'tr")
for i in threetothree:
    threenames.append(i[0].split(",")[0])
    threenames.append(i[1].split(",")[0])
    for j in i[0:-1]:
        print(j,end=',')
    print(i[-1])
print("5'tr to 3'tr")
for i in fivethreeres:
    fivenames.append(i[0].split(",")[0])
    threenames.append(i[1].split(",")[0])
    for j in i[0:-1]:
        print(j,end=',')
    print(i[-1])
print("5'tr to good")
for i in fivefullres:
    allnames.append(i[1].split(",")[0])
    fivenames.append(i[0].split(",")[0])
    for j in i[0:-1]:
        print(j,end=',')
    print(i[-1])
print("3'tr to good")
for i in threefullres:
    allnames.append(i[1].split(",")[0])
    threenames.append(i[0].split(",")[0])
    for j in i[0:-1]:
        print(j,end=',')
    print(i[-1])
# Create sets of results to remove insertions that have been recombining multiple times etc.
duplicationtsds=[]
with open('duplication_sum.csv','r') as f:
    for line in f:
        line=line.lstrip("\ufeff").rstrip("\n").split(",")
        if line[-2]=="-":
            line[-2]=""
        if line[-3]=="-":
            line[-3]=""
        if line[-1]=="t" or line[-1]=="m":
            line[-3]="/".join(singlemut(tuple([line[-3][-4:],line[-3][-5:],line[-3][-6:],line[-3][-7:],line[-3][-8:]])))
            line[-2]="/".join(singlemut(tuple([line[-2][:4],line[-2][:5],line[-2][:6],line[-2][:7],line[-2][:8]])))
        duplicationtsds.append(tuple(line))
f.close()

# Good to segdup

duplicationsfound=[]
for j in duplicationtsds:
    duplicationsfound.append([j])
    for i in normal:
        for z in j[1:]:
            z=z.split("/")
            if i[-1] in z:
                duplicationsfound[-1].append(i)
                break
print("dup to good")
for i in duplicationsfound:
    if len(i)>1:
        a=0
        print(i[0][0]+",,,,"+i[0][1]+"/"+i[0][2],end=',')
        for j in i[1:]:

            allnames.append(j[0][0])
            print(j[0][0]+","+j[0][1].split(",")[3]+","+j[0][1].split(",")[5]+","+j[0][1].split(",")[6]+","+j[0][1].split(",")[16]+","+j[0][1].split(",")[17]+","+j[0][1].split(",")[18])
            a+=1



print("")

# Five to segdup

duplicationsfound=[]
for j in duplicationtsds:
    duplicationsfound.append([j])
    for i in five:
        for z in j[1:]:
            z=z.split("/")
            if any(g in z for g in i[2]):
                duplicationsfound[-1].append(i)
                break


print("dup to 5tr")
for i in duplicationsfound:
    if len(i)>1:
        a=0
        print(i[0][0]+",,,,"+i[0][1]+"/"+i[0][2],end=',')

        for j in i[1:]:
            longest=[]
            if any(elem in i[0][1].split("/")  for elem in j[2]):
                longest+=intersection(i[0][1].split("/"), j[2])
            if any(elem in i[0][2].split("/") for elem in j[2]):
                longest+=intersection(i[0][2].split("/"), j[2])
            fivenames.append(j[0][0])
            print(j[0][0]+","+j[0][1].split(",")[3]+","+j[0][1].split(",")[5]+","+j[0][1].split(",")[6]+","+j[0][1].split(",")[16]+","+j[0][1].split(",")[17]+","+j[0][1].split(",")[18]+","+max(longest))
            a+=1



# three to segdup
duplicationsfound=[]
for j in duplicationtsds:
    duplicationsfound.append([j])
    for i in three:
        for z in j[1:]:
            z=z.split("/")
            if any(g in z for g in i[1]):
                duplicationsfound[-1].append(i)
                break
#print(duplicationsfound)
print("dup to 3tr")
for i in duplicationsfound:
    if len(i)>1:
        a=0
        print(i[0][0]+",,,,"+i[0][1]+"/"+i[0][2],end=',')
        for j in i[1:]:
            longest=[]
            if any(elem in i[0][1].split("/")  for elem in j[1]):
                longest+=intersection(i[0][1].split("/"), j[1])
            if any(elem in i[0][2].split("/") for elem in j[1]):
                longest+=intersection(i[0][2].split("/"), j[1])
            threenames.append(j[0][0])
            print(j[0][0]+","+j[0][1].split(",")[3]+","+j[0][1].split(",")[5]+","+j[0][1].split(",")[6]+","+j[0][1].split(",")[16]+","+j[0][1].split(",")[17]+","+j[0][1].split(",")[18]+","+max(longest))
            a+=1
# recombinators to segdup


duplicationsfound=[]
for j in duplicationtsds:
    duplicationsfound.append([j])
    for i in recombinators:
        for z in j[1:]:
            z=z.split("/")
            if any(g in z for g in i[1]) or any(g in z for g in i[2]):
                duplicationsfound[-1].append(i)
                break

print("dup to recom")
for i in duplicationsfound:
    if len(i)>1:
        a=0
        print(i[0][0]+",,,,"+i[0][1]+"/"+i[0][2],end=',')
        for j in i[1:]:
            longest=[]
            if any(elem in i[0][1].split("/")  for elem in j[1]):
                longest+=intersection(i[0][1].split("/"), j[1])
            if any(elem in i[0][2].split("/") for elem in j[1]):
                longest+=intersection(i[0][2].split("/"), j[1])
            if any(elem in i[0][1].split("/")  for elem in j[2]):
                longest+=intersection(i[0][1].split("/"), j[2])
            if any(elem in i[0][2].split("/") for elem in j[2]):
                longest+=intersection(i[0][2].split("/"), j[2])


            recombpairsres.append(j[0][0])
            print(j[0][0]+","+j[0][1].split(",")[3]+","+j[0][1].split(",")[5]+","+j[0][1].split(",")[6]+","+j[0][1].split(",")[16]+","+j[0][1].split(",")[17]+","+j[0][1].split(",")[18]+","+max(longest))
            a+=1

for i in normaltonormal:

    for j in i:
        allnames.append(j.split(",")[0])
recombpairsres=list(set(recombpairsres))

allnames=list(set(allnames))
fivenames=list(set(fivenames))
threenames=list(set(threenames))





# Print results in form of a table.

print("--------------")



print("1) Matching TSDs")
print("2) Non-matching TSDs")
print("3) 5' truncaded insertions")
print("4) 3' truncaded insertions")
print("5) 5' and 3' truncaded insertions")
print("--------------------------")
print("")
print("")
a1=tsdn-(len(allnames))-duptsd
a2=notsdn-len(recombpairsres)-dupnotsd
a3=fiven-len(fivenames)-dupfive
a4=threen-len(threenames)-dupthree
a5=fntn-dupfnt
b1=duptsd
b2=dupnotsd
b3=dupfive
b4=dupthree
b5=dupfnt
c1=len(allnames)
c2=len(recombpairsres)
c3=len(fivenames)
c4=len(threenames)
c5=0



print("")
print("                             |  1  |  2  |  3  |  4  |  5  |")
print("-----------------------------|-----|-----|-----|-----|-----|")
print("No evidence of recombination |                             |")
print("External flanks don't match  | "+str(a1)+" |"+str(a2)+" |"+str(a3)+" |"+str(a4)+" |"+str(a5)+" |")
print("                             |                             |")
print("-----------------------------|-----------------------------|")
print("External flanks match        |                             |")
print("Segmental duplications       | "+str(b1)+" |"+str(b2)+" |"+str(b3)+" |"+str(b4)+" |"+str(b5)+" |")
print("                             |                             |")
print("-----------------------------|-----------------------------|")
print("Insertions sharing TSDs      |                             |")
print("Ext flanks don't match       | "+str(c1)+" |"+str(c2)+" |"+str(c3)+" |"+str(c4)+" |"+str(c5)+" |")
print("Recombinations               |                             |")
print("-----------------------------|-----------------------------|")
j=0
