HERVK_pipeline - reference part

Tools used to extract the target sequences from the reference genome and analyse them. For detailed process steps and settings, please see each script and adjust input files/parameters internally, according to commets. See schem.png for further reference of each step.

1) Step 1 - extraction of target sequences.


 extractfromref_optimal.py - script, that extracts sequences from desired locations from .fasta reference file. It takes reference genome and a list of desired locations (toextrac.left.txt) as input, outputs fasta formatted sequences extracted from the reference on screen.

2) Step 2 - alignment of 20bp of insertion end + 20bp flanking sequence using Geneious, HERV-K113(AY037928) LTR used as a reference (AY037928.fasta; 1-968bp)

3) Step 3 - BLAST detection.

BLAST online was used to detect precise locations within 10kb of truncated insertions (default online blastn parameters: eval=0.05, word_size=11, max_target_seq=100, filter low complexity regions, mask for lookup table only, M/M score 2/-3, Gap cost 5/2)

4) Step 4 - BLAT reference check.

BLAT online was used to ascertain locations for each insertion in GRCh36/Hg18, GRCh37/Hg19, GRCh38/Hg38 Human reference genome versions.

5) Step 5 - Identifying TSDs

TSDs were listed for each insertion in csv format (flanks and TSD status in last 3 columns), incluting the possiblility of single nucleotide mutation - generated from initial csv using tsd_wooble.py + manual check.
