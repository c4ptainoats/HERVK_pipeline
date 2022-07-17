HERVK_pipeline - detecting duplications.

Tools used to detect groups of duplications involved with HERV-K(HML2) elements.

1) Step 1 - get duplication list


Refer to duplications directory and readme.


2) Step 2 - get TSD comparisons.

Script tsd.py analyses all possible combinations of TSDs, it takes 1bp mutations into account and input data in csv format, with insertions analysed in every line, 5' and 3' flanks in 15 and 16 columns respectively and TSD status in column 17 (5' TRUNCATED, 3' TRUNCATED, - for no recognizable TSD or the TSD sequence). Outputs all the numeric data for all the combinations in a table and outputs the combinations in csv format for further analysis.

3) Step 3 - optional - further selection.

You can analyse the resulting combinations and select the best combinations for all different types of TSDs. In the work, longest combinations between no recognizable TSD and insertions displaying good TSDs were considered.

4) Step 4 - Statistical analysis: Kolomogorov-Smirnov test

The results can be plotted against the chromosome on which every detected insertion is located. Such distribution can be compared to the general distribution of insertions using Kolomogorov-Smirnov goodness-of-fit test, using e.g. KSINV function from Real Statistics Resource Pack (https://real-statistics.com; file XRealStats-Mac.xlam comes from the website and installs the resource for mac version of Excel, tested with 365). The Kolmogorov-Smirnov D-statistic was calculated maximum vertical distance between the empirical cumulative distribution of observed donor insertions and the cumulative distribution of all HERV-K(HML2) loci :
D = max(F(x)-G(x)),
where F(x) is the distribution of selected inserts and G(x) is the distribution of all inserts.

5) Insertion rate in humans
