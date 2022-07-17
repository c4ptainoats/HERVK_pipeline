HERVK_pipeline - Finding Human-specific Insertions

1) List of all insertions existing in Mhudiblu_PPA_v2 were provided by Hayley Free,
who used BLAST based search of PanPan reference against AY037928.fasta (1-968bp):

makeblastdb -in [reference genome assembly] -dbtype nucl -out [database file name]

blastn -db [database file name] -query [K113 LTR] -word_size 10 -out [blast output name] -outfmt 6

The data was loaded to MS excel, labelled plus orientation if subject end > subject start. Then removed all SVA elements, that aligned only with 1-339bp and 687-809bp of AY037928 LTR.
Eventually all remaining sequences ends were extracted (extractfromref_optimal.py) realigned and locations were provided.

2) Extracting elements - 3' end

extractfromref_optimal.py was used to extract 3' ends of each element: 100bp of 3' end and 100bp of 3' flank; from the Mhudiblu_PPA_v2 reference and all known 100bp 3' ends + 3' flanks for the human elements from GRCh38.

3) Ends comparisons

Extracted ends were compared using blast:

makeblastdb -in allhervk3.fasta -dbtype nucl
blastn -word_size 11 -evalue 0.001 -db allhervk3.fasta -query allPanP3.fasta -outfmt 6
-max_target_seqs 500000 -out Panp3.txt

allhervk3.fasta - database of HERV-K 3'human ends+flanks(100bp each)

allPanP3.fasta - database of HERV-K 3'PanPan ends+flanks(100bp each)

4) Collapse duplications

discard.py was used to leave only single BLAST result from duplications coming up in the p3 search. Input used was  Panp3.txt BLAST result and Human HERV-K duplications list, dups.csv

5) Results were filtered in excel  according to following criteria:

eval	<=0.00001
precent ident	>= 90
matchlen	>=120
q_e (query end)	>=119

6) Resulting matched Human elements were combined with non-reference Human HERV-K(HML2) elements to produce an alignment of HERV-K(HML2) human-specific insertions. In case of Full-length elements, only corresponding LTRs were extracted
