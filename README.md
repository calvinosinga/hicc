# hicc
hicc stands for HI and Color Correlations, the project uses the Illustris TNG simulation to study the clustering bias in the statistics of atomic hydrogen and optical galaxy samples. This depends on the Pylians library from Francisco Villaescusa-Navarro.

## hiptl
This folder has the code that creates the fields for the atomic hydrogen prescriptions applied to the particle catalogue. It does so for both redshift space and in real space. Stored in 2048^3 grids.

TODO: add generalization for boxlength
changed the path to TNG data on deepthought2

## galaxy
This folder creates the fields for the subhalos separated into by optical color into blue and red populations, after removing unresolved subhalos. It does so for both redshift and real space. Stored in 2048^3 grids.

TODO: add generalization for boxlength
changed the path to TNG data on deepthought2

## power_spectrum
This folder calculates the various power spectra from the fields generated.

TODO: add generalization for boxlength
changed the path to TNG data on deepthought2
