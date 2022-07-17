HERVK_pipeline - detecting duplications.

Tools used to detect groups of duplications involved with HERV-K(HML2) elements.

1) Step 1 - extraction of target sequences + 500bp of flanks.


 extractfromref_optimal.py - script, that extracts sequences from desired locations from .fasta reference file. It takes reference genome and a list of desired locations (toextrac.left.txt) as input, outputs fasta formatted sequences extracted from the reference on screen.

2) Step 2 - get the flanks

Script splitends.py gets 500bp of flanking sequences from each end (5' and 3') from fasta files separately.

3) Step 3 - BLAST detection.

Local BLAST is used to detect duplications separately for 5' flanks and 3' flanks, first making databases from each set of flanks, as follows:

makeblastdb -in flank.fasta -dbtype nucl
blastn -word_size 11 -evalue 0.001 -db flank.fasta -query flank.fasta -outfmt 6
-max_target_seqs 500000 -out output.txt

4) Step 4 - clean the blast results

cleanup.py is used to remove AvA - type duplications.

5) Step 5 - group results

group.py is used to create groups of same query-subject requence pairs.

5) Step 5 - check flanks

To remove some false-positives, selector.py is used for double-checking the results and including only ones that overlap the 10bp immmidiate flank.

6) Step 6 - check TSD.

To double-check the TSD homology between the selected candidates for duplications, comparetsd_segdupv2.py is used. It needs csv file with ins name in column 1, tsd status list for each insertion in column 19 (TRUNCATED for truncated insertions/- for no recognizable TSD), flanks in columnts 17 and 18 (5' and 3'; recommended 10bp).

7) Step 7 - connect 5' and 3' results.

Script connectends.py checks for duplications showing homology both on 5' and 3' flank.

8) Step 8 - create list of duplications.

Script getduplist.py creates a list of groups of duplications from previous output.

9) Step 9 - check remaining results and realign.

Script getoneside.py creates a list of results that show homology only on 5' or 3' flank, these can be manually checked in e.g. Geneious for some rare duplications that could be filtered out automatically and duplications of truncated/mutated elements, where tsd is difficult to recognize.
