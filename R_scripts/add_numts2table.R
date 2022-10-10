#!/usr/bin/env Rscript --vanilla

################
## add_numts2table.R rod_numts_hits_BLNG_20_BLNT_220614.tab new_numts_hits_BLNG_20_BLNT_220614.tab
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

library("dplyr")

p.name='add_numts2table.R'

args <- commandArgs(trailingOnly=TRUE)
cat(paste("#",p.name,paste(args,collapse=' ')),'\n')

plabel=paste(c(p.name,"\n",args),collapse=' ', sep=' ')

if (! is.na(args[1])) {
   ref.gb.stats=args[1]
} else {
    write(sprintf("%s supp_table1.txt  file required",p.name),stderr())
    quit(save='n')
}

if (! is.na(args[2])) {
   mtgenome.class.info=args[2]
} else {
    write(sprintf("%s all_mtGenome_class_size  file required",p.name),stderr())
    quit(save='n')
}

## get table info
gt.data <- read.table(file=ref.gb.stats,sep='\t',header=TRUE,stringsAsFactors=TRUE)
## dt.gt.data <- read.table(file=dt.tax.stats,sep='\t',header=TRUE,stringsAsFactors=TRUE)

sn.data <- read.table(file=mtgenome.class.info,header=TRUE,sep='\t')

hit_files = args[3:length(args)]
if (length(hit_files)<1) {
    write(sprintf("%s hits files required",p.name),stderr())
    quit(save='no')
}

data_dir = "data/"

if (!exists("read.numt.hits",mode='function')) source("./R_scripts/read_numt_hits.R")

save.cols=c('search_id','tag','canon_acc','common_name','qlen','s_acc','percid','alen','expect','taxon_id')

summ_data = NULL

## it makes more sense to aggregate first, then add matrix and algo

for (file in hit_files) {

    n_parts<-strsplit(file,'.', fixed=TRUE)
    t_name <- n_parts[[1]][1]

    file_parts <- strsplit(t_name,'_', fixed=TRUE)
    m_name <- file_parts[[1]][4]
    m_algo <- file_parts[[1]][6]
    aln_str <- file_parts [[1]][5]

    ## need to shift algo to matrix for BLN/BLNT

    message(paste("file=",file))
    message(paste(t_name, m_name, m_algo, aln_str))

    hit_data = read.numt.hits(data_dir, file, save.cols)

    t_data <- aggregate(cbind(numt_cnt=qlen) ~ taxon_id, data=hit_data, FUN=NROW)

    wt_data <- aggregate(cbind(numt_len=alen) ~ taxon_id, data=hit_data, FUN=sum)

    s5_data <- aggregate(alen ~ taxon_id, data=hit_data, FUN=fivenum)

    t_data = merge(t_data, wt_data[,c('taxon_id','numt_len')],by='taxon_id')

    med_mx_s5 = s5_data[,c('taxon_id')]
    med_mx_s5 = cbind(med_mx_s5,s5_data$alen[,3],s5_data$alen[,5])
    colnames(med_mx_s5) = c('taxon_id','med_len','max_len')

    ## print(head(med_mx_s5))

    t_data = merge(t_data, med_mx_s5, by='taxon_id')

    ## print("head(t_data1)")
    ## print(head(t_data))

    summ_data <- rbind(summ_data, t_data)
}

no_data <- subset(gt.data,!(NCBI_Tax_ID %in% summ_data$taxon_id))
message('no_data:')
message(no_data)

if (NROW(no_data) > 0) {
    zero_data <- data.frame(taxon_id=no_data$NCBI_Tax_ID, numt_cnt=0, numt_len=0, med_len=0, max_len=0)
    message("zero_data")
    message(zero_data)
    summ_data <- rbind(summ_data, zero_data)
}

summ_data_s <- summ_data[order(summ_data$taxon_id),]
summ_data_s <- summ_data_s[!duplicated(summ_data_s$taxon_id),c('taxon_id','numt_cnt','numt_len','med_len','max_len')]

## head(gt.data)
## head(summ_data_s)

new.gt.data <- merge(gt.data,summ_data_s,by.x='NCBI_Tax_ID',by.y='taxon_id')

## head(new.gt.data)
## head(sn.data)

new.gt.data <- merge(new.gt.data,sn.data[,c('taxon_id','s_name')],by.x='NCBI_Tax_ID',by.y='taxon_id')

## head(new.gt.data)
new_col_order = c('NCBI_Tax_ID','s_name','Species','Common_name','Mitochondrial_genome_accession','Nuclear_genome_accession','Nuclear_genome_size_.GB.','RefSeq.GenBank','dataset','numt_cnt','numt_len','med_len','max_len')

new.gt.data <- new.gt.data[order(new.gt.data$s_name),]

write.table(new.gt.data[,new_col_order],file='',quote=FALSE, sep='\t',row.names=FALSE, col.names=TRUE)
