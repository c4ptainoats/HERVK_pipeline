HERVK_pipeline - Bayesian phylogeny

The   phylogeny   of   HERV-K(HML2)   elements   was   constructed   based   on   an alignment of all HERV-K(HML2) LTRs, including solo-LTR elements and full-length elements (where multiple  LTRs  from single  full-length element were  marked according to their position in the element).
The alignment was performed in Geneious   software   and   to   diminish   influence   of   random   mutations,   sites containing 95% of gaps were removed from the alignment. The alignment was loaded into BEAST v1.10.4, with following custom settings:

Bayesian skyline tree model (skyline.groupSize=15);
HKY with gamma distributed rate variation (G=4) substitution model   (empirical   base   frequencies);  
UPGMA   starting   tree;
uncorrelated, relaxed lognormal clock model;
250M MCMC chain length.

The trees were ran in BEAST using the following command:

 ./beast -beagle settings.xml

Three   chains   were   ran   and   combined   together   to   produce   satisfactory   ESS (>150, 70M burn-in, adjusted according to the trace), using logcombiner (part of BEAST   v1.10.4,   GUI).   


The   tree   annotation   of   individual   MCMC   chains   was performed in treeannotator (part of BEAST v1.10.4), using following command:

./treeannotator all_ltr.trees all_ltr.tree


Script tanglegram.R can be used to produce tanglegrams of corresponding fragments of trees, see comments inside.

Script stattrees.py can be used to analyse distribution of posterior probabilities from input trees and produce a statistic in csv format.
