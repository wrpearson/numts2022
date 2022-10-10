#!/usr/bin/env Rscript --vanilla

## pplot_numts_alen_mat0.R -- 0 versions of plot programs assume no self-matches
##

################
## pplot_numts_alen_mat.R -- plot box plots of all data by matrix and coverage -- .tab files in data/, results in figs/
##
################
##
## does a plot of numbers of numts as a function of alen, rather than coverage.
## to compare blastn vs tfx, need to convert to "universal" (nt) alignment length
##
################

library("ggplot2")
library("cowplot")
library("RColorBrewer")

p.name<-'pplot_numts_hits_alen_mat0.R'

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

if (!exists("read.numt.hits",mode='function')) source("./R_scripts/read_numt_hits.R")

save.cols=c('search_id','tag','canon_acc','common_name','qlen','s_acc','percid','alen','expect','taxon_id')

summ_data = NULL

## it makes more sense to aggregate first, then add matrix and algo

a_lens = c(30,60,150,300)

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

    ## this would be much easier if:
    ## (1) we adjusted the dataframe$alen based on TFX
    if (m_algo == 'TFX') {
        hit_data$alen = hit_data$alen*3
    }

    ## (2) for combined TFX/BLNT results, we check for mtGenome
    if (m_algo == 'TFXBLNT') {
        hit_data[hit_data$common_name != 'mtGenome',]$alen <- hit_data[hit_data$common_name != 'mtGenome',]$alen * 3
    }

    for (a_len in a_lens) {

        a_hit_data = hit_data[hit_data$alen >= a_len,]

        ## cat(sprintf("%s %d %d / %d\n",m_algo, a_len, NROW(a_hit_data),NROW(hit_data)))
        ## before aggregating (or perhaps during), we need calculate hit counts for longer alignment lengths

        t_data <- aggregate(cbind(h_cnt=qlen) ~ taxon_id, data=a_hit_data, FUN=NROW)
        ## print(head(t_data))

        t_data$matrix <- m_name
        t_data$algo <- m_algo
        t_data$a_len <- a_len

        ## check for zero-hit taxa
        no_data <- subset(gt.data,!(taxon_id %in% t_data$taxon_id))
        ## head(no_data)
        if (NROW(no_data) > 0) {
            zero_data <- data.frame(taxon_id=no_data$taxon_id, h_cnt=0, matrix=m_name, algo=m_algo, a_len=a_len)
            print("head(zero_data)")
            print(head(zero_data))
            t_data <- rbind(t_data, zero_data)
        }

        ## print("head(t_data):")
        ## print(head(t_data))

        summ_data <- rbind(summ_data, t_data)
    }
}

## is it refseq vs genbank

## cat("head(summ_data)/gt.data:\n")
## head(summ_data)
## head(gt.data)

summ_data <- merge(summ_data, gt.data, by='taxon_id')

## head(summ_data)

summ_data$alen_f <- sprintf("%.0f",summ_data$a_len)
summ_data$alen_f <- factor(summ_data$alen_f,levels=a_lens,ordered=TRUE)

## head(summ_data)

alen_fs <- sprintf("%s alignment length: >= %.0f nt",c('A.','B.','C.','D.'), a_lens)
names(alen_fs) = a_lens

g_types = factor(c('refseq','genbank'),ordered=TRUE)
summ_data$g_type = factor(summ_data$g_type, levels=g_types, labels=g_types, ordered=TRUE)

a_types = factor(c('TFX','BLN','BLNT','TFXBLNT'), ordered=TRUE)
summ_data$algo = factor(summ_data$algo, levels=a_types, labels=a_types,ordered=TRUE)

## cat("before hc_meds aggregate\n")
hc_meds <- aggregate(h_cnt ~ g_type + matrix + alen_f, data=summ_data, median)
hc_meds$h_cnt = trunc(hc_meds$h_cnt)
## hc_meds

if (nrow(summ_data[summ_data$h_cnt < 1,])>0) {
    summ_data[summ_data$h_cnt < 1,]$h_cnt = 0.1
    ## print(summ_data[summ_data$h_cnt < 1,])
}

g.colors <- brewer.pal(6,'Dark2') # 'Dark2', 'Set2', 'Paired'
g.shapes <- c(16,16,15,17)

a.shapes <- c(1,0)

gs.color <- scale_color_manual(values=g.colors,labels=c("RefSeq","GenBank"))
gs.shape <- scale_shape_manual(values=g.shapes)
gs.fill <- scale_fill_manual(values=alpha(g.colors[1:2], 0.00),labels=c("RefSeq","GenBank"))

as.color <- scale_color_manual(values=g.colors[1:2])
as.shape <- scale_shape_manual(values=a.shapes)
as.fill <- scale_fill_manual(values=g.colors[1:2])

y.limits = c(0.1,10000)

matrix.names=c('BLN','BLNT','MD10','MD20','MD40','BP62','BLNG','BLNGT','MD40BG')
summ_data$matrix = factor(summ_data$matrix, levels=matrix.names, labels=matrix.names, ordered=TRUE)
m.shapes = c('BLN'=3,'BLNT'=4, 'MD10'=1, 'MD20'=6,'MD40'=16, 'BP62'=17,'BLNG'=0, 'BLNGT'=5,'MD40BG'=12)

m.shape <- scale_shape_manual(values=m.shapes,labels=matrix.names,guide='none')

## summary(summ_data)

pos.dodge = position_dodge(width=0.9)

mlab_coord = c(6000,3000)

pA <- ggplot(summ_data,aes(x=matrix,y=h_cnt, shape=matrix, fill=g_type)) +
  ## geom_violin(position=pos.dodge, width=0.4) + geom_boxplot(outlier.shape=NA, position=pos.dodge, width=0.4) +
  geom_point(aes(shape=matrix, color=g_type),size=2,alpha=0.7, position=position_jitterdodge(jitter.width=0.6,jitter.height=0.02)) +
  geom_boxplot(outlier.shape=NA) +
    ## geom_text(data=hc_meds,aes(x=matrix, group=g_type, y = h_cnt*1.1, label=h_cnt), vjust= 0, size=2.5, position=position_dodge(0.8)) +
    geom_text(data=hc_meds,aes(x=matrix, group=g_type, y = mlab_coord[2], label=h_cnt), vjust= 0, size=3.0, position=position_dodge(0.8)) +
  scale_y_log10(name='numt count',breaks=c(0.1,1,10,100,1000,10000),labels=c('0','1','10','100','1000','10000')) +
  gs.color + gs.fill + m.shape +
  theme_bw() + theme(legend.justification='top',legend.title=element_blank()) + theme(strip.text=element_text(hjust=0)) +
  annotate("text", x = 0.5, y = mlab_coord[1], label = 'median:',hjust=0, vjust=0,size=3.0) +
  scale_x_discrete(name="search strategy and scoring matrix") + facet_wrap(~alen_f,labeller=labeller(alen_f=alen_fs),nrow=2) # +ggtitle(paste("A. by matrix type",plot.title))


if (Sys.getenv(c("PLOT_DOLAB"))=='none') {
    pdf(file=paste0(pdf_dir, t_name, "_alen_mat0_NL.pdf"),width=10.0, height=6.0)
    pA
} else {
    pdf(file=paste0(pdf_dir, t_name, "_alen_mat0.pdf"),width=10.0, height=6.0)
    pA + labs(caption=plabel)
}
