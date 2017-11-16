# astral_pipeline
A pipeline for running pyRAD output through ASTRAL, a program for calculating a species tree under the multi-species coalescent model.  The pipeline was originally written to work with ASTRAL-II and has not yet been tested using more recent versions.  

The pipeline uses a naive binning approach to group short-read loci together at a user-defined clustering threshold prior to gene tree estimation in RAxML.  

This pipeline was used for analyses in a forthcoming paper: MR Bangs, MR Douglas, SM Mussmann, and ME Douglas (2017) Unraveling historical introgression and resolving phylogenetic discord within *Catostomus* (Pisces: Catostomidae).

The code for ASTRAL can be found here: https://github.com/smirarab/ASTRAL

The code for RAxML can be found here: https://github.com/stamatak/standard-RAxML

More information about naive binning can be found here: https://academic.oup.com/bioinformatics/article/29/18/2277/240280
