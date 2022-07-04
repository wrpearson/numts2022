#!/usr/bin/env Rscript --vanilla

## plot_numts_hits_percid_mat0.R -- 0 versions of plot programs assume no self-matches
##

################
## plot_numts_hits_percid_mat.R -- plot box plots of all data by matrix and coverage -- .tab files in data/, results in figs/
##
################
##
## does a plot of numbers of numts as a function of percid, rather than alen.
##
################

library("ggplot2")
library("cowplot")
library("RColorBrewer")

p.name<-'plot_numts_hits_percid_mat0.R'

data_dir="data/"
pdf_dir="figs/"

args <- commandArgs(trailingOnly=TRUE)

plabel=sprintf("%s\n%s",p.name, paste(args,collapse=' ', sep=' '))

if (! is.na(args[1])) {
   ref.gb.stats=args[1]
} else {
    write(sprintf("%s new_mtGenome_accs_class_size.tab  file required",p.name),stderr())
    quit(save='n')
}

## get genome info
gt.data <- read.table(file=ref.gb.stats,sep='\t',header=TRUE,stringsAsFactors=TRUE)
## dt.gt.data <- read.table(file=dt.tax.stats,sep='\t',header=TRUE,stringsAsFactors=TRUE)

g_types = c('refseq','genbank')
g_types = factor(g_types, levels=g_types, ordered=TRUE)


hit_files = args[-1]
if (length(hit_files)<1) {
    write(sprintf("%s hits files required",p.name),stderr())
    quit(save='no')
}

data_dir = "data/"
pdf_dir = "figs/"

if (!exists("read.numt.hits",mode='function')) source("read_numt_hits.R")

save.cols=c('search_id','tag','canon_acc','common_name','qlen','s_acc','percid','alen','expect','taxon_id')

summ_data = NULL

## it makes more sense to aggregate first, then add matrix and algo

percids = c(0.25,0.70,0.85,0.95)

for (file in hit_files) {

    n_parts<-strsplit(file,'.', fixed=TRUE)
    t_name <- n_parts[[1]][1]

    file_parts <- strsplit(t_name,'_', fixed=TRUE)
    m_name0 <- file_parts[[1]][4]
    m_algo <- file_parts[[1]][6]
    aln_str <- file_parts [[1]][5]

    ## need to shift algo to matrix for BLN/BLNT

    print(paste("file=",file))
    print(paste(t_name, m_name0, m_algo, aln_str))

    hit_data = read.numt.hits(data_dir, file, save.cols)

    if (m_algo == 'BLNT') {
        m_name <- paste0(m_name0,'T')
    } else {
        m_name <- m_name0
    }

    for (percid in percids) {

        p_hit_data = hit_data[hit_data$percid >= percid,]
        print(head(p_hit_data))

        ## before aggregating (or perhaps during), we need calculate hit counts for longer alignment lengths

        t_data <- aggregate(cbind(h_cnt=qlen) ~ taxon_id, data=p_hit_data, FUN=NROW)
        print(head(t_data))

        t_data$matrix <- m_name
        t_data$algo <- m_algo
        t_data$percid <- percid

        ## check for zero-hit taxa
        no_data <- subset(gt.data,!(taxon_id %in% t_data$taxon_id))
        head(no_data)
        if (NROW(no_data) > 0) {
            zero_data <- data.frame(taxon_id=no_data$taxon_id, h_cnt=0, matrix=m_name, algo=m_algo, percid=percid)
            print("head(zero_data)")
            print(head(zero_data))
            t_data <- rbind(t_data, zero_data)
        }

        print("head(t_data):")
        print(head(t_data))

        summ_data <- rbind(summ_data, t_data)
    }
}

## is it refseq vs genbank

cat("head(summ_data)/gt.data:\n")
head(summ_data)
head(gt.data)

summ_data <- merge(summ_data, gt.data, by='taxon_id')

head(summ_data)

summ_data$percid_f <- sprintf("%.2f",summ_data$percid)

summ_data$percid_f <- factor(summ_data$percid_f,levels=sprintf("%.2f",percids),ordered=TRUE)

g_types = factor(c('refseq','genbank'),ordered=TRUE)
summ_data$g_type = factor(summ_data$g_type, levels=g_types, labels=g_types, ordered=TRUE)

a_types = factor(c('TFX','BLN'), ordered=TRUE)
summ_data$algo = factor(summ_data$algo, levels=a_types, labels=a_types,ordered=TRUE)

cat("before hc_meds aggregate\n")

hc_meds <- aggregate(h_cnt ~ g_type + matrix + percid_f, data=summ_data, median)
hc_meds$h_cnt = trunc(hc_meds$h_cnt)

hc_meds

if (nrow(summ_data[summ_data$h_cnt < 1,])>0) {
    summ_data[summ_data$h_cnt < 1,]$h_cnt = 0.1
    print(summ_data[summ_data$h_cnt < 1,])
}

g.colors <- brewer.pal(6,'Dark2') # 'Dark2', 'Set2', 'Paired'
g.shapes <- c(16,16,15,17)

a.shapes <- c(1,0)

gs.color <- scale_color_manual(values=g.colors)
gs.shape <- scale_shape_manual(values=g.shapes)
gs.fill <- scale_fill_manual(values=alpha(g.colors[1:2], 0.00))

as.color <- scale_color_manual(values=g.colors[1:2])
as.shape <- scale_shape_manual(values=a.shapes)
as.fill <- scale_fill_manual(values=g.colors[1:2])

y.limits = c(0.1,1000)

matrix.names=c('BLN','BLNT','MD10','MD20','MD40','BP62','BLNG','BLNGT')
summ_data$matrix = factor(summ_data$matrix, levels=matrix.names, labels=matrix.names, ordered=TRUE)
summary(summ_data)

plot.title=sprintf("(percid: %s%%)", unique(summ_data$percid_f))

pdf(file=paste0(pdf_dir, t_name, "_percid_mat0.pdf"),width=8.0, height=6.0)

pos.dodge = position_dodge(width=0.9)

pA <- ggplot(summ_data,aes(x=matrix,y=h_cnt, shape=g_type, fill=g_type)) +
  ## geom_violin(position=pos.dodge, width=0.4) + geom_boxplot(outlier.shape=NA, position=pos.dodge, width=0.4) +
  geom_point(aes(shape=g_type, color=g_type),size=2,alpha=0.4, position=position_jitterdodge(jitter.width=0.6,jitter.height=0.02)) +
  geom_boxplot(outlier.shape=NA) +
    ## geom_text(data=hc_meds,aes(x=matrix, group=g_type, y = h_cnt*1.1, label=h_cnt), vjust= 0, size=2.5, position=position_dodge(0.8)) +
    geom_text(data=hc_meds,aes(x=matrix, group=g_type, y = 1300, label=h_cnt), vjust= 0, size=3.0, position=position_dodge(0.8)) +
  scale_y_log10(name='numt count',breaks=c(0.1,1,10,100,1000),labels=c('0','1','10','100','1000')) +
  gs.color + gs.shape + gs.fill +
  theme_bw() + theme(legend.position="none") +
  annotate("text", x = 0.5, y = 2200, label = 'median:',hjust=0, vjust=0,size=3.0) +
  scale_x_discrete(name="Matrix") + labs(caption=plabel) + facet_wrap(~percid_f,nrow=2) # +ggtitle(paste("A. by matrix type",plot.title))

pA
