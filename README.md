
10-October-2022

Scripts and datasets available for "Comparison of detection methods
and genome quality when quantifying nuclear mitochondrial insertions
in vertebrate genomes", revised submission to Frontiers in Genetics, October,
2022.

These files include all of the shell (.sh) and 'R' scripts, as well as
the datasets, used to make the figures and Supplementary Table 1 for the
manuscript.  These figures and table are created running the script:

```
rod_new_hit_aln_plots_20220930_s.sh  220930 yes
```

R Scripts were run using R version 4.2.1 (2022-06-23) on Macbook Pro
(Intel) under macOS Monterey (version 12.6).

R libraries required are:

library("RColorBrewer")
library("cowplot")
library("dplyr")
library("ggplot2")
library("ggtree")
library("ggtreeExtra")
library("patchwork")


