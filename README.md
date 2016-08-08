Subtest
=========

A package for testing for differential genetic basis between two putative disease subtypes

An important consideration in genomic analysis is the question of whether disease heterogeneity arises from differences in genetic causality or differences in environment. Furthermore, accounting for genetically-driven disease heterogeneity can strengthen the ability to detect disease-causative variants. The question of whether putative subtypes of a disease phenotype have differential genetic basis can be addressed by analysis of SNP case-control data from GWAS and similar studies. This package implements a proposed method and generates supporting data.

We analyse disease heterogeneity as characterised either by a division of the case group into two subgroups or parametrised by a single quantitative variable. Our overall approach is to compute two sets of GWAS summary statistics, in the form of Z-scores:

*Z<sub>d<sub>*: characterising phenotypic heterogeneity (independent of controls)

*Z<sub>a<sub>*: comparing cases with controls

We then test for the presence of a set of SNPs which which show simultaneous evidence of inflation in both scores, corresponding to association both with case/control status and with within-case hetergeneity.

We accomplish this by fitting a bivariate mixture Gaussian model under two hypotheses (null and full). The difference in fit of the two models, and hence the evidence for differential genetic basis of subgroups, is assessed by means of an adapted likelihood-ratio test.

Broadly, the package performs four main functions: generation of Z scores, fitting of models, assessment of models by simulation of random subtypes, and analysis of single SNPs. Input should be a SnpMatrix object (package SnpStats) which has been rigorously QC'd and either pruned to minimal linkage disequilibrium between SNPs or had weights calculated by the LDAK algorithm ([http://dougspeed.com/ldak/]).

Output includes a p-value for the evidence of differential genetic basis in subgroups (under the null hypothesis of independence of causative basis of disease heterogeneity and disease causality) and details of fitted models. Fitted values give some indication of the genetic architecture of the disease; namely the approximate proportions of null SNPs and SNPs which are disease-causative without differentiating subtypes, and the (multivariate) distribution of effect sizes of causative SNPs.

Further information can be foundin our paper at http://biorxiv.org/content/early/2016/08/02/037713.
