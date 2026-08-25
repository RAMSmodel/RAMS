# RAMS
[![DOI](https://zenodo.org/badge/307815774.svg)](https://zenodo.org/badge/latestdoi/307815774)

RAMS model source code. To see full documentation, go to the [van den Heever group website](https://vandenheever.atmos.colostate.edu/vdhpage/rams/rams_docs.php).

In this fork, I have added BUGSRAD radiation as option 4 for ISWRTYP and ILWRTYP. A few notes about this implementation:

o  Only cloud droplets, pristine ice and snow are passed to the BUGSRAD radiation routines. This is ok for now as cloud droplets and ice crystals are the two most important hydrometeors for radiation

o  Cloud droplet mass and number are both passed. Only the sum of pristine ice and snow mass is passed.

o  No aerosol species are coupled to BUGSRAD radiation.

o  BUGSRAD contains separate albedo values for direct and diffuse radiation. I've set these to be equal. (see radiate/driver_read.f90, alndr, alndf, alvdr, alvdf)

o  BUGSRAD contains separate albedo values for shortwave and near IR. I've set these to be equal. (aln* - near IR, alv* - visible)

o  BUGSRAD has its own precision definitions (e.g. dbl_kind). They're still being used and aren't interfaced with the options to compile with double or single precision in the Makefile.

I've also added RTE+RRTMGP radiation (https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2019MS001621) as option 5

