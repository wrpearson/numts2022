#!/bin/sh

## same as rod_new_hit_aln_plots_20220930.sh, except uses the R_scripts directory
##

## this version uses:
## plot_numts_hits_tree3_mat0.R
## plot_numts_hits_gsize_mat_g.R
## plot_numts_hits_alen_mat0.R

## which do not need to have the data combined
## their command line syntax is:
## plot_numts_hits_tree3_mat0.R tree.nwk new|rod_numts_huts_MAT_alen_ALG_2205xx[_m].tab
## for each of the MAT/ALG being plotted

## for plot_numts_hits_gsize_mat_g.R and plot_numts_hits_alen_mat0.R,
## a genome summary file (new_mtGenome_accs_class_size.tab,
## rod_mtGenome_accs_class_size.tab) is used

## the new plot programs look for all the taxa in the
## tree.nwk/class_size.tab files and add h_cnt=0 for missing results

## check for timestamp

BDIR=./R_scripts

if [[ -z $1 ]]; then
    echo timestamp required
    exit
fi

tstamp=$1

if [[ $2 ]]; then
    export PLOT_DOLAB=$2
    echo $PLOT_DOLAB
fi

new_mtgenome_class_file='new_mtGenome_accs_class_size_20220930.tab'

## plot distributions of numts vs tree for rod
## Fig. 1

echo "$BDIR/pplot_numts_hits_tree3_leg.R rod_mtGenomes.nwk rod_numts_hits_BLN_30_BLN_${tstamp}_m.tab ..."
$BDIR/pplot_numts_hits_tree3_leg.R rod_mtGenomes.nwk rod_numts_hits_BLN_30_BLN_${tstamp}_m.tab rod_numts_hits_BLNG_30_BLN_${tstamp}.tab rod_numts_hits_BLN_30_BLNT_${tstamp}_m.tab rod_numts_hits_BLNG_30_BLNT_${tstamp}.tab rod_numts_hits_MD10_10_TFX_${tstamp}_m.tab rod_numts_hits_MD40_10_TFX_${tstamp}_m.tab rod_numts_hits_MD40BG_10_TFXBLNT_${tstamp}.tab 

if [[ $PLOT_DOLAB == 'none' ]] ; then
    cp -p figs/rod_numts_hits_MD40BG_10_TFXBLNT_${tstamp}_tree3_leg_p_NL.pdf Frontiers_text/Figures/fig1_rod_tree3_${tstamp}_NL.pdf
fi

## plot distributions of numts vs genome size for new
## Fig. 5

echo "$BDIR/pplot_numts_hits_gsize_mat_g2c.R $new_mtgenome_class_file new_numts_hits_BLN_30_BLN_${tstamp}_m.tab ..."
$BDIR/pplot_numts_hits_gsize_mat_g2c.R $new_mtgenome_class_file new_numts_hits_BLN_30_BLN_${tstamp}_m.tab new_numts_hits_BLNG_30_BLN_${tstamp}.tab new_numts_hits_BLN_30_BLNT_${tstamp}_m.tab new_numts_hits_BLNG_30_BLNT_${tstamp}.tab new_numts_hits_MD10_10_TFX_${tstamp}_m.tab new_numts_hits_MD40_10_TFX_${tstamp}_m.tab new_numts_hits_MD40BG_10_TFXBLNT_${tstamp}.tab 

if [[ $PLOT_DOLAB == 'none' ]]; then
    cp -p figs/new_numts_hits_MD40BG_10_TFXBLNT_${tstamp}_gsize_mat_g2c_p_NL.pdf Frontiers_text/Figures/fig5_new_gsize_${tstamp}_NL.pdf
fi

## plot numt length vs genome size for new
## Suppl. Fig. 4

echo "$BDIR/pplot_numts_hits_gsize_alen_g2c.R $new_mtgenome_class_file new_numts_hits_BLN_30_BLN_${tstamp}_m.tab ..."
$BDIR/pplot_numts_hits_gsize_alen_g2c.R $new_mtgenome_class_file new_numts_hits_BLNG_30_BLNT_${tstamp}.tab new_numts_hits_MD40_10_TFX_${tstamp}_m.tab new_numts_hits_MD40BG_10_TFXBLNT_${tstamp}.tab

if [[ $PLOT_DOLAB == 'none' ]]; then
    cp -p figs/new_numts_hits_MD40BG_10_TFXBLNT_${tstamp}_gsize_alen_g2c_p_NL.pdf Frontiers_text/Figures/suppl_fig4_new_gsize_median${tstamp}_NL.pdf
fi

## plot distributions of numts vs alignment length
## Fig. 2

echo "$BDIR/pplot_numts_hits_alen_mat0.R $new_mtgenome_class_file new_numts_hits_BLN_30_BLN_${tstamp}_m.tab ..."
$BDIR/pplot_numts_hits_alen_mat0.R $new_mtgenome_class_file new_numts_hits_BLN_30_BLN_${tstamp}_m.tab new_numts_hits_BLNG_30_BLN_${tstamp}.tab new_numts_hits_BLN_30_BLNT_${tstamp}_m.tab new_numts_hits_BLNG_30_BLNT_${tstamp}.tab new_numts_hits_MD10_10_TFX_${tstamp}_m.tab new_numts_hits_MD40_10_TFX_${tstamp}_m.tab new_numts_hits_MD40BG_10_TFXBLNT_${tstamp}.tab 
if [[ $PLOT_DOLAB == 'none' ]]; then
    cp -p figs/new_numts_hits_MD40BG_10_TFXBLNT_${tstamp}_alen_mat0_NL.pdf Frontiers_text/Figures/suppl_fig1_new_alen_${tstamp}_NL.pdf
fi


echo "$BDIR/pplot_numts_hits_alen_mat0.R rod_mtGenome_accs_class_size.tab rod_numts_hits_BLN_30_BLN_${tstamp}_m.tab ..."
$BDIR/pplot_numts_hits_alen_mat0.R rod_mtGenome_accs_class_size.tab rod_numts_hits_BLN_30_BLN_${tstamp}_m.tab rod_numts_hits_BLNG_30_BLN_${tstamp}.tab rod_numts_hits_BLN_30_BLNT_${tstamp}_m.tab rod_numts_hits_BLNG_30_BLNT_${tstamp}.tab rod_numts_hits_MD10_10_TFX_${tstamp}_m.tab rod_numts_hits_MD40_10_TFX_${tstamp}_m.tab rod_numts_hits_MD40BG_10_TFXBLNT_${tstamp}.tab

if [[ $PLOT_DOLAB == 'none' ]]; then
    cp -p figs/rod_numts_hits_MD40BG_10_TFXBLNT_${tstamp}_alen_mat0_NL.pdf Frontiers_text/Figures/fig2_rod_alen_${tstamp}_NL.pdf
fi


## Fig. 4 plot numts per gene

echo "$BDIR/pplot_numts_hits_order1rs.R rod_mtGenome_accs_class_size.tab rod_numts_hits_MD40_10_TFX_${tstamp}.tab"
$BDIR/pplot_numts_hits_order1rs.R rod_mtGenome_accs_class_size.tab rod_numts_hits_MD40_10_TFX_${tstamp}.tab

if [[ $PLOT_DOLAB == 'none' ]]; then
    cp -p figs/rod_numts_hits_MD40_10_TFX_${tstamp}_ord1rs_p_NL.pdf Frontiers_text/Figures/fig4_rod_pergene_${tstamp}_NL.pdf
fi

## Fig. 3

echo "$BDIR/pplot_diff_percid_hist4.R new_numts_hits_MD40BLN_10_TFXBLNT_${tstamp}.tab rod_numts_hits_MD40BLN_10_TFXBLNT_${tstamp}.tab"
$BDIR/pplot_diff_percid_hist4.R new_numts_hits_MD40BLN_10_TFXBLNT_${tstamp}.tab rod_numts_hits_MD40BLN_10_TFXBLNT_${tstamp}.tab

if [[ $PLOT_DOLAB == 'none' ]]; then
    cp -p figs/rod_numts_hits_MD40BLN_10_TFXBLNT_${tstamp}_diff_alen_hist4_NL.pdf  Frontiers_text/Figures/fig3_rod_new_percid_${tstamp}_NL.pdf
fi

## suppl_Fig2

echo "(no bad rodents) $BDIR/pplot_numts_hits_alen_mat0x.R rod_mtGenome_accs_class_size.tab new_numts_hits_BLN_30_BLN_${tstamp}_m.tab ..."

export JOIN_PTS=''
export EXCLUDE_TAXA='bad_rodents.tab'
$BDIR/pplot_numts_hits_alen_mat0x.R rod_mtGenome_accs_class_size.tab rod_numts_hits_BLN_30_BLN_${tstamp}_m.tab rod_numts_hits_BLNG_30_BLN_${tstamp}.tab rod_numts_hits_BLN_30_BLNT_${tstamp}_m.tab rod_numts_hits_BLNG_30_BLNT_${tstamp}.tab rod_numts_hits_MD10_10_TFX_${tstamp}_m.tab rod_numts_hits_MD40_10_TFX_${tstamp}_m.tab rod_numts_hits_MD40BG_10_TFXBLNT_${tstamp}.tab 

EXCLUDE_TAXA=''

if [[ $PLOT_DOLAB == 'none' ]]; then
    cp -p figs/new_numts_hits_MD40BG_10_TFXBLNT_${tstamp}_alen_mat0x_NL.pdf Frontiers_text/Figures/suppl_fig2_good_rod_alen_${tstamp}_NL.pdf
fi

$BDIR/add_numts2table_hs.R supp_Table1_common_name_20220930.txt all_mtGenome_accs_class_size_20220930.tab rod_numts_hits_MD40BG_10_TFXBLNT_${tstamp}.tab new_numts_hits_MD40BG_10_TFXBLNT_${tstamp}.tab > supp_Table1_numt_cnt_len_hs.tab
if [[ $PLOT_DOLAB == 'none' ]]; then
    cp -p supp_Table1_numt_cnt_len_hs.tab Frontiers_text/Tables
fi

## suppl_Fig3 box plot across genome intervals

echo "$BDIR/pplot_numts_hits_mtgenome.R rod_numts_hits_BLNG_30_BLNT_${tstamp}_s.tab"
$BDIR/pplot_numts_hits_mtgenome.R rod_mtGenome_accs_class_size.tab rod_numts_hits_BLNG_30_BLNT_${tstamp}_s.tab

if [[ $PLOT_DOLAB == 'none' ]]; then
    cp -p figs/rod_numts_hits_BLNG_30_BLNT_${tstamp}_s_mtg_NL.pdf Frontiers_text/Figures/suppl_fig3_rod_segment_${tstamp}_NL.pdf
fi
