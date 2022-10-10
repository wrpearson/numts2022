#!/usr/bin/env Rscript --vanilla

################
## corr_gize_all.R [supp_Table1_numt_cnt_len.tab]
##
## this version makes a new dataframe with the summ of the hits from all the queries
## sorts plot based on genome size
## _g version prints genome size
##
## does aggregation on scoring matrix as well as taxon_id
##
## modified to insert zero's for taxa in new_mtGenome_accs_class_size.tab
##
################

p.name<-'corr_gsize_all.R'

args <- commandArgs(trailingOnly=TRUE)

plabel=paste(c(p.name,"\n",args),collapse=' ', sep=' ')

if (is.na(args[1])) {
  gt.file <- 'supp_Table1_numt_cnt_len.tab'
} else {
  gt.file <- args[1]
}

## get genome info
all_data <- read.table(file=gt.file,sep='\t',header=TRUE,stringsAsFactors=TRUE)
## dt.gt.data <- read.table(file=dt.tax.stats,sep='\t',header=TRUE,stringsAsFactors=TRUE)

head(all_data[,c('NCBI_Tax_ID','s_name','Nuclear_genome_size_.GB.','RefSeq.GenBank','numt_cnt')])

## head(summ_data)

zero_numt_taxa = c(7955, 8364, 31033, 35670)

all_data_nz <- all_data[! all_data$NCBI_Tax_ID %in% zero_numt_taxa,]

all_ref <- all_data[tolower(all_data$RefSeq.GenBank) == 'refseq',]

all_ref_nz <- all_data_nz[tolower(all_data_nz$RefSeq.GenBank) == 'refseq',]

cat(sprintf("all: %d\n",NROW(all_data)))
print(summary(all_data))
pearson_corr = cor.test(all_data$Nuclear_genome_size_.GB.,all_data$numt_cnt, method='pearson')
print(pearson_corr)

spearman_corr = cor.test(all_data$Nuclear_genome_size_.GB.,all_data$numt_cnt, method='spearman')
print(spearman_corr)

kendall_corr = cor.test(all_data$Nuclear_genome_size_.GB.,all_data$numt_cnt, method='kendall')
print(kendall_corr)

cat(sprintf("all_nz: %d\n",NROW(all_data_nz)))
print(summary(all_data_nz))
pearson_corr = cor.test(all_data_nz$Nuclear_genome_size_.GB.,all_data_nz$numt_cnt, method='pearson')
print(pearson_corr)

spearman_corr = cor.test(all_data_nz$Nuclear_genome_size_.GB.,all_data_nz$numt_cnt, method='spearman')
print(spearman_corr)

kendall_corr = cor.test(all_data_nz$Nuclear_genome_size_.GB.,all_data_nz$numt_cnt, method='kendall')
print(kendall_corr)

## refseq

cat(sprintf("refseq_all: %d\n",NROW(all_ref)))
print(summary(all_ref))

pearson_corr = cor.test(all_ref$Nuclear_genome_size_.GB.,all_ref$numt_cnt, method='pearson')
print(pearson_corr)

spearman_corr = cor.test(all_ref$Nuclear_genome_size_.GB.,all_ref$numt_cnt, method='spearman')
print(spearman_corr)

kendall_corr = cor.test(all_ref$Nuclear_genome_size_.GB.,all_ref$numt_cnt, method='kendall')
print(kendall_corr)

## refseq_nz

cat(sprintf("refseq_nz: %d\n",NROW(all_ref_nz)))
print(summary(all_ref_nz))

pearson_corr = cor.test(all_ref_nz$Nuclear_genome_size_.GB.,all_ref_nz$numt_cnt, method='pearson')
print(pearson_corr)

spearman_corr = cor.test(all_ref_nz$Nuclear_genome_size_.GB.,all_ref_nz$numt_cnt, method='spearman')
print(spearman_corr)

kendall_corr = cor.test(all_ref_nz$Nuclear_genome_size_.GB.,all_ref_nz$numt_cnt, method='kendall')
print(kendall_corr)
