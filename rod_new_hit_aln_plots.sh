#!/bin/sh

## this version uses:
## pplot_numts_hits_tree3_mat0.R
## pplot_numts_hits_gsize_mat_g2c.R
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

if [[ -z $1 ]]; then
    echo timestamp required
    exit
fi

tstamp=$1

## to avoid plotting plot program/argument in plot, PLOT_DOLAB='none'

if [[ $2 ]]; then
    export PLOT_DOLAB=$2
    echo $PLOT_DOLAB
fi

## plot distributions of numts vs tree for rod
## Fig. 1

echo "pplot_numts_hits_tree3_leg.R rod_mtGenomes.nwk rod_numts_hits_BLN_30_BLN_${tstamp}_m.tab ..."
./pplot_numts_hits_tree3_leg.R rod_mtGenomes.nwk rod_numts_hits_BLN_30_BLN_${tstamp}_m.tab rod_numts_hits_BLNG_30_BLN_${tstamp}.tab rod_numts_hits_BLN_30_BLNT_${tstamp}_m.tab rod_numts_hits_BLNG_30_BLNT_${tstamp}.tab rod_numts_hits_MD10_10_TFX_${tstamp}_m.tab rod_numts_hits_MD40_10_TFX_${tstamp}_m.tab rod_numts_hits_MD40BG_10_TFXBLNT_${tstamp}.tab 

if [[ $PLOT_DOLAB == 'none' ]] ; then
    cp -p figs/rod_numts_hits_MD40BG_10_TFXBLNT_${tstamp}_tree3_leg_p_NL.pdf Frontiers_text/Figures/fig1_rod_tree3_${tstamp}_NL.pdf
fi

## plot distributions of numts vs genome size for new
## Fig. 5

echo "pplot_numts_hits_gsize_mat_g2c.R new_mtGenome_accs_class_size.tab new_numts_hits_BLN_30_BLN_${tstamp}_m.tab ..."
./pplot_numts_hits_gsize_mat_g2c.R new_mtGenome_accs_class_size.tab new_numts_hits_BLN_30_BLN_${tstamp}_m.tab new_numts_hits_BLNG_30_BLN_${tstamp}.tab new_numts_hits_BLN_30_BLNT_${tstamp}_m.tab new_numts_hits_BLNG_30_BLNT_${tstamp}.tab new_numts_hits_MD10_10_TFX_${tstamp}_m.tab new_numts_hits_MD40_10_TFX_${tstamp}_m.tab new_numts_hits_MD40BG_10_TFXBLNT_${tstamp}.tab 

if [[ $PLOT_DOLAB == 'none' ]]; then
    cp -p figs/new_numts_hits_MD40BG_10_TFXBLNT_${tstamp}_gsize_mat_g2c_p_NL.pdf Frontiers_text/Figures/fig5_new_gsize_${tstamp}_NL.pdf
fi

## plot distributions of numts vs alignment length
## Fig. 2

echo "./plot_numts_hits_alen_mat0.R rod_mtGenome_accs_class_size.tab rod_numts_hits_BLN_30_BLN_${tstamp}_m.tab ..."
./plot_numts_hits_alen_mat0.R rod_mtGenome_accs_class_size.tab rod_numts_hits_BLN_30_BLN_${tstamp}_m.tab rod_numts_hits_BLNG_30_BLN_${tstamp}.tab rod_numts_hits_BLN_30_BLNT_${tstamp}_m.tab rod_numts_hits_BLNG_30_BLNT_${tstamp}.tab rod_numts_hits_MD10_10_TFX_${tstamp}_m.tab rod_numts_hits_MD40_10_TFX_${tstamp}_m.tab rod_numts_hits_MD40BG_10_TFXBLNT_${tstamp}.tab

if [[ $PLOT_DOLAB == 'none' ]]; then
    cp -p figs/rod_numts_hits_MD40BG_10_TFXBLNT_${tstamp}_alen_mat0_NL.pdf Frontiers_text/Figures/fig2_rod_alen_${tstamp}_NL.pdf
fi

echo "./plot_numts_hits_alen_mat0.R new_mtGenome_accs_class_size.tab new_numts_hits_BLN_30_BLN_${tstamp}_m.tab ..."
./plot_numts_hits_alen_mat0.R new_mtGenome_accs_class_size.tab new_numts_hits_BLN_30_BLN_${tstamp}_m.tab new_numts_hits_BLNG_30_BLN_${tstamp}.tab new_numts_hits_BLN_30_BLNT_${tstamp}_m.tab new_numts_hits_BLNG_30_BLNT_${tstamp}.tab new_numts_hits_MD10_10_TFX_${tstamp}_m.tab new_numts_hits_MD40_10_TFX_${tstamp}_m.tab new_numts_hits_MD40BG_10_TFXBLNT_${tstamp}.tab 
if [[ $PLOT_DOLAB == 'none' ]]; then
    cp -p figs/new_numts_hits_MD40BG_10_TFXBLNT_${tstamp}_alen_mat0_NL.pdf Frontiers_text/Figures/suppl_fig1_new_alen_${tstamp}_NL.pdf
fi

## Fig. 3

echo "./plot_diff_percid_hist4.R new_numts_hits_MD40BLN_10_TFXBLNT_${tstamp}.tab rod_numts_hits_MD40BLN_10_TFXBLNT_${tstamp}.tab"
./plot_diff_percid_hist4.R new_numts_hits_MD40BLN_10_TFXBLNT_${tstamp}.tab rod_numts_hits_MD40BLN_10_TFXBLNT_${tstamp}.tab

if [[ $PLOT_DOLAB == 'none' ]]; then
    cp -p figs/rod_numts_hits_MD40BLN_10_TFXBLNT_${tstamp}_diff_alen_hist4_NL.pdf  Frontiers_text/Figures/fig3_rod_new_percid_${tstamp}_NL.pdf
fi

## Fig. 4 plot numts per gene

echo "./pplot_numts_hits_order1rs.R rod_mtGenome_accs_class_size.tab rod_numts_hits_MD40_10_TFX_${tstamp}.tab"
./pplot_numts_hits_order1rs.R rod_mtGenome_accs_class_size.tab rod_numts_hits_MD40_10_TFX_${tstamp}.tab

if [[ $PLOT_DOLAB == 'none' ]]; then
    cp -p figs/rod_numts_hits_MD40_10_TFX_${tstamp}_ord1rs_p_NL.pdf Frontiers_text/Figures/fig4_rod_pergene_${tstamp}_NL.pdf
fi

## suppl_Fig2

echo "./plot_numts_hits_alen_mat0x.R new_mtGenome_accs_class_size.tab new_numts_hits_BLN_30_BLN_${tstamp}_m.tab ..."

export JOIN_PTS=''
export EXCLUDE_TAXA='bad_rodents.tab'
./plot_numts_hits_alen_mat0x.R rod_mtGenome_accs_class_size.tab rod_numts_hits_BLN_30_BLN_${tstamp}_m.tab rod_numts_hits_BLNG_30_BLN_${tstamp}.tab rod_numts_hits_BLN_30_BLNT_${tstamp}_m.tab rod_numts_hits_BLNG_30_BLNT_${tstamp}.tab rod_numts_hits_MD10_10_TFX_${tstamp}_m.tab rod_numts_hits_MD40_10_TFX_${tstamp}_m.tab rod_numts_hits_MD40BG_10_TFXBLNT_${tstamp}.tab 

EXCLUDE_TAXA=''

if [[ $PLOT_DOLAB == 'none' ]]; then
    cp -p figs/rod_numts_hits_MD40BG_10_TFXBLNT_${tstamp}_alen_mat0x_NL.pdf Frontiers_text/Figures/suppl_fig2_rod_alen_${tstamp}_NL.pdf
fi


add_numts2table.R supp_Table_1.txt rod_numts_hits_MD40BG_10_TFXBLNT_${tstamp}.tab new_numts_hits_MD40BG_10_TFXBLNT_${tstamp}.tab > supp_Table1_numt_cnt_len.tab
if [[ $PLOT_DOLAB == 'none' ]]; then
    cp -p supp_Table1_numt_cnt_len.tab Frontiers_text/Tables
fi

