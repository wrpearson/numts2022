#!/usr/bin/env Rscript --vanilla

################
## used to plot Figure 3

################
## plot_diff_percid_hist4.R new_numts_hits_MD40BLN_10_TFXBLNT_220603.tab
##
##
################

library("ggplot2")
library("cowplot")
library("RColorBrewer")

p.name<-'plot_diff_percid_hist4.R'

data_dir="data/"
pdf_dir="figs/"

## read in genbank/refseq labels

ref.gb.stats <- 'all_mtGenome_accs_class_size.tab'

gr_data <- read.table(ref.gb.stats,sep='\t',header=TRUE)

args <- commandArgs(trailingOnly=TRUE)

## print(args)

hit_files = args
if (length(hit_files)<1) {
    write(sprintf("%s hits files required",p.name),stderr())
    quit(save='no')
}

plabel = sprintf("%s\n%s",p.name, paste(args,collapse=' ', sep=' '))

if (!exists("read.numt.hits",mode='function')) source("read_numt_hits.R")

save.cols=c('search_id','tag','canon_acc','common_name','qlen','s_acc','percid','alen','expect','taxon_id','src')

summ_data = NULL

## it makes more sense to aggregate first, then add matrix and algo

for (file in hit_files) {

    n_parts<-strsplit(file,'.', fixed=TRUE)
    t_name <- n_parts[[1]][1]

    file_parts <- strsplit(t_name,'_', fixed=TRUE)

    data_name <- file_parts[[1]][1]
    m_name <- file_parts[[1]][4]
    m_algo <- file_parts[[1]][6]
    aln_str <- file_parts [[1]][5]

    ## need to shift algo to matrix for BLN/BLNT

    print(paste("file=",file))
    print(paste(t_name, data_name, m_name, m_algo, aln_str))

    t_data = read.numt.hits(data_dir, file, save.cols)
    t_data$d_name <- data_name

    ## print(head(t_data))

    m_name_l = strsplit(t_data$tag,'_',fixed=TRUE)
    ## print(head(m_name_l))

    m_name_sv = unlist(lapply(m_name_l,function(x) {x[length(x)]}))
    ## print(head(m_name_sv))

    t_data$f_algo = m_name_sv

    t_data$alen_f = t_data$alen

    ## print(fivenum(t_data$alen_f))

    t_data[t_data$f_algo == 'tfx',]$alen_f = t_data[t_data$f_algo == 'tfx',]$alen_f*3

    ## print(fivenum(t_data$alen_f))

    ## print("head(t_data)")
    ## print(head(t_data))

    summ_data <- rbind(summ_data, t_data)
}

## print(unique(summ_data$algo))

summ_data <- merge(summ_data, gr_data, by='taxon_id')

summ_data$src  <- factor(summ_data$src,levels=c('both','d1','d2'),ordered=TRUE)
## head(summ_data)
## summary(summ_data)

## aggregate(alen_f ~ s_name + src, data=summ_data, FUN=fivenum)

## aggregate(alen_f ~ src, data=summ_data, FUN=fivenum)

hs_summ_data <- summ_data[summ_data$taxon_id==9606,]
mm_summ_data <- summ_data[summ_data$taxon_id==10090,]
## head(mm_summ_data)

unique(summ_data$d_name)

rod_summ_data <- summ_data[summ_data$d_name=='rod' & summ_data$g_type=='refseq',]
head(rod_summ_data)

vert_summ_data <- summ_data[summ_data$d_name=='new' & summ_data$g_type=='refseq',]
head(vert_summ_data)

## head(h_summ_data)

s.colors <- brewer.pal(4,'Dark2') # 'Dark2', 'Set2', 'Paired'

sp.color <- scale_color_manual(values=c('both'=s.colors[1],'d1'=s.colors[2],'d2'=s.colors[3]), labels=c('both','TFASTX/MD40','BLASTN genomic (BLNGT)'))

pMM <- ggplot(mm_summ_data,aes(x=percid*100,color=src)) + geom_histogram(alpha=0.4,position='identity',bins=20,fill=I('transparent')) + scale_x_continuous(limits=c(20,100)) +  sp.color + theme_bw() + labs(title='A. mouse',x='percent identity')+ theme(plot.title=element_text(size=12),legend.position=c(0.33,0.85),legend.title=element_blank(), legend.text=element_text(size=9),legend.key.size=unit(0.35,'cm')) + scale_y_continuous(limits=c(0,50))

pHS <- ggplot(hs_summ_data,aes(x=percid*100,color=src)) + geom_histogram(alpha=0.2,position='identity',bins=20,fill=I('transparent')) + scale_x_continuous(limits=c(20,100)) +  sp.color + theme_bw() + labs(title='C. human', x='percent identity') + theme(plot.title=element_text(size=12), legend.position='none', legend.title=element_blank(), legend.text=element_text(size=9),legend.key.size=unit(0.35,'cm')) ## + scale_y_continuous(limits=c(0,200))

pV <- ggplot(vert_summ_data,aes(x=percid*100,color=src)) + geom_histogram(alpha=0.4,position='identity',bins=20,fill=I('transparent')) + scale_x_continuous(limits=c(20,100)) +  sp.color + theme_bw() + labs(title='D. vertebrates',x='percent identity')+ theme(plot.title=element_text(size=12),legend.position='none',legend.title=element_blank(), legend.text=element_text(size=9),legend.key.size=unit(0.35,'cm')) ##  + scale_y_continuous(limits=c(0,2000))

pR <- ggplot(rod_summ_data,aes(x=percid*100,color=src)) + geom_histogram(alpha=0.2,position='identity',bins=20,fill=I('transparent')) + scale_x_continuous(limits=c(20,100)) +  sp.color + theme_bw() + labs(title='B. rodents', x='percent identity') + theme(plot.title=element_text(size=12), legend.position='none', legend.title=element_blank(), legend.text=element_text(size=9),legend.key.size=unit(0.35,'cm')) + scale_y_continuous(limits=c(0,500))

pE <- ggplot() + labs(caption=plabel) + theme_void() + theme(plot.caption=element_text(size=8))

if (Sys.getenv(c("PLOT_DOLAB"))=='none') {
    pdf(file=paste0(pdf_dir, t_name, "_diff_alen_hist4_NL.pdf"),width=7, height=6)
    plot_grid(pMM,pR, pHS, pV,ncol=2,rel_heights=c(0.4,0.4))
} else {
    pdf(file=paste0(pdf_dir, t_name, "_diff_alen_hist4.pdf"),width=7, height=6.2)
    plot_grid(pMM,pR,pHS,pV,NULL, pE,ncol=2,rel_heights=c(0.4, 0.4, 0.04))
}
