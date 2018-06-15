# Astral Pipeline
A pipeline for running pyRAD output through ASTRAL, a program for calculating a species tree under the multi-species coalescent model.  The pipeline was originally written to work with ASTRAL-II and has not yet been tested using more recent versions.  

The pipeline uses a naive binning approach to group short-read loci together at a user-defined clustering threshold prior to gene tree estimation in RAxML.  

The code for ASTRAL can be found here: https://github.com/smirarab/ASTRAL

The code for RAxML can be found here: https://github.com/stamatak/standard-RAxML

More information about naive binning can be found here: https://academic.oup.com/bioinformatics/article/29/18/2277/240280

## Installation:

A guide explaining how to install and use this pipeline is coming soon.  In short, you will have to modify the paths in perl script to reflect the location of your astral installation and possibly your RAxML installation (depending on how your system is configured).

## Using astral_pipeline

A user guide is coming soon.  Currently you can view user options for the pipeline by executing the perl script without any command line options.

## Citing astral_pipeline

If you use this method, please cite:

Bangs, Max R., Marlis R. Douglas, Steven M. Mussmann, and Michael E. Douglas. 2018. Unraveling historical introgression and resolving phylogenetic discord within *Catostomus* (Osteichthys: Catostomidae). BMC Evolutionary Biology 18:86. doi: https://doi.org/10.1186/s12862-018-1197-y
