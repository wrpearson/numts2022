#!/usr/bin/env Rscript --vanilla

####
## plot0.R assumes no self-matches -- does not correct

################
## plot_numts_hits_tree2_mat0.R  all_mito_geonmes_rr.nwk numts_summ_BLN-MD40_00_TFX_n.tab *
##
## this script reads the tree, then a set of hits files (with different scoring matrices/algorithms)
## aggregating the hits into counts
################
##
## modified to add zero's to missing entries based on tree organisms
##

library("ggplot2")
library("cowplot")
library("RColorBrewer")
library("ggnewscale")
library("ggtree")
library("ggtreeExtra")
## library("ape")

p.name<-'plot_numts_hits_tree2_mat0_wgt.R'

args <- commandArgs(trailingOnly=TRUE)

plabel=paste(c(p.name,"\n",args),collapse=' ', sep=' ')

if (! is.na(args[1])) {
   tree_file=args[1]
} else {
    write(sprintf("%s tree file required",p.name),stderr())
    quit(save='n')
}

## read the tree

print(paste("tree_file:",tree_file))
tr <- read.tree(tree_file)

ref.gb.stats <- 'mtGenome_accs.tab'
## get genome info
gt.data <- read.table(file=ref.gb.stats,sep='\t',header=TRUE,stringsAsFactors=TRUE)
## dt.gt.data <- read.table(file=dt.tax.stats,sep='\t',header=TRUE,stringsAsFactors=TRUE)
head(gt.data)

## label tree
g_types = c('refseq','genbank')
g_types = factor(g_types, levels=g_types, ordered=TRUE)

dd <- data.frame(label = tr$tip.label)
dd <- merge(dd, gt.data[,c('mtGenome_acc','s_name','taxon_id','g_type')], by.x='label', by.y='mtGenome_acc')
dd$s_name2 <- gsub("_",". ",dd$s_name)
## dd$name_acc <- sprintf("%s : %s",dd$label, dd$s_name2)

cat("head dd:\n")
head(dd)

dd$g_type = factor(dd$g_type, levels=g_types, ordered=TRUE)
## dd[dd$s_name=='H_sapiens',]$g_type=g_types[3]
## dd[dd$s_name=='M_muntjak',]$g_type=g_types[4]

##cat("dd\n")
## str(dd)
##dd[,c('label','s_name')]

################
## read in the data

data_dir = "data/"
pdf_dir = "figs/"

hit_files = args[-1]
if (length(hit_files)<1) {
    write(sprintf("%s hits files required",p.name),stderr())
    quit(save='no')
}

if (!exists("read.numt.hits",mode='function')) source("read_numt_hits.R")

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

    ## need to shift algo to matrix for BLN/BLNTe

    print(paste("file=",file))
    print(paste(t_name, m_name, m_algo, aln_str))

    hit_data = read.numt.hits(data_dir, file, save.cols)

    t_data <- aggregate(cbind(wh_cnt=alen) ~ taxon_id, data=hit_data, FUN=sum)

    if (m_algo == 'BLNT') {
        m_name <- paste0(m_name,'T')
    }

    t_data$matrix <- m_name
    t_data$algo <- m_algo
    t_data$l_aln <- aln_str

    no_data <- subset(dd,!(taxon_id %in% t_data$taxon_id))

    if (NROW(no_data) > 0) {
        zero_data <- data.frame(taxon_id=no_data$taxon_id, wh_cnt=0, matrix=m_name, algo=m_algo, l_aln=aln_str)
        print("zero_data")
        print(zero_data)
        t_data <- rbind(t_data, zero_data)
    }

    print("head(t_data)")
    print(head(t_data))

    summ_data <- rbind(summ_data, t_data)
}

head(summ_data)
summary(summ_data)

## removed for plot0.R
## all_data$h_cnt <- all_data$h_cnt - 1

m_types = c('BLN', 'BLNT', 'MD10','MD20','MD40','BP62','BLNG', 'BLNGT','MD40BG')
m.shapes = c('BLN'=3,'BLNT'=4, 'MD10'=1, 'MD20'=6,'MD40'=16, 'BP62'=17,'BLNG'=0, 'BLNGT'=5,'MD40BG'=12)

um_types <- unique(summ_data$matrix)
um_types

summ_data$matrix = factor(summ_data$matrix,levels=um_types,ordered=TRUE)

head(summ_data)
str(summ_data)

## plot.title=sprintf("(matrix: %s; coverage: %d%%)",m_name, unique(all_data$cov)*100)

## gt_data <- read.table('mtGenomes_taxIDs_20210330.f1-4',sep='\t',header=TRUE)
## str(g_types)
## print(c(g_types[3],g_types[4]))

summ_data <- merge(summ_data, gt.data[,c('taxon_id','g_type','mtGenome_acc','s_name')], by='taxon_id')

if (nrow(summ_data[summ_data$wh_cnt < 1,])>0) {
    summ_data[summ_data$wh_cnt < 1,]$wh_cnt = 0.1
    print(summ_data[summ_data$wh_cnt < 1,])
}

## cat("head(summ_data):\n")
## head(summ_data)

##print(paste("dim(gt.data):",dim(gt.data),"dim(summ_data):",dim(summ_data)))

## length(unique(all_data$taxon_id))

## print(all_data[all_data$s_name=='H_sapiens',])

summ_data$g_type = factor(summ_data$g_type, levels=g_types, labels=g_types, ordered=TRUE)
summ_data$matrix = factor(summ_data$matrix, levels=m_types, labels=m_types, ordered=TRUE)

## str(summ_data)

## colors for genome quality
br.gq_ref.colors <- brewer.pal(8,'GnBu')
gq.ref.colors <- br.gq_ref.colors[4:7]
br.gq_gb.colors <- brewer.pal(8,'OrRd')
gq.gb.colors <- br.gq_gb.colors[4:7]

g.colors <- brewer.pal(4,'Dark2') # 'Dark2', 'Set2', 'Paired'
g.shapes <- c(16,16,15,17)

s.color <- scale_color_manual(values=g.colors)
s.shape <- scale_shape_manual(values=g.shapes)

m.shape <- scale_shape_manual(values=m.shapes,labels=m_types)

pdf(file=paste0(pdf_dir, t_name, "_tree_m0_wgt.pdf"),width=8, height=5)

## cat("**tr$tip.label\n")
## str(tr)
## tr$tip.label

## cat("**dd\n")
## dd

pA <- ggtree(tr, layout='rectangular') %<+% dd
pA <- rotate(pA,40)

## cat("pA$data\n")
## print(pA$data[pA$data$isTip,c('y','label','s_name','g_type')],n=28)

## str(pA$data)

tip.dd <- pA$data[pA$data$isTip,c('y','label','s_name','s_name2')]
tip.dd$s_name = as.character(tip.dd$s_name)

## print(tip.dd,n=nrow(tip.dd))

## if xlim() is too low, then the tree will not plot.
## but xlim() for the rodents should be lower than xlim() for the vertebrates

pA.l <- pA + geom_tiplab(aes(label=s_name2, color=g_type),align=TRUE,size=4,offset=0.0,angle=90) + xlim(NA,0.35) +
        s.color + theme(legend.position="none") + # geom_text(aes(label=node)) +
        coord_flip()

## tip.dd[order(tip.dd$y, decreasing=FALSE),]

## cat("length(tip.dd)\n")
## str(tip.dd)
## length(tip.dd$y)

## cat("tip.dd levels\n")
## tip.dd[order(tip.dd$y,decreasing=FALSE),]$s_name

length(unique(summ_data$s_name))

summ_data$s_name <- factor(summ_data$s_name, levels=tip.dd[order(tip.dd$y,decreasing=FALSE),]$s_name, ordered=TRUE)

## summ_data$s_name <- factor(summ_data$s_name, levels=get_taxa_name(pA.l), ordered=TRUE)

## str(summ_data)

pB <- ggplot(summ_data,aes(x=s_name,y=wh_cnt)) +
   geom_jitter(aes(color=g_type,shape=matrix),size=2,alpha=0.8,width=0.1) +
    scale_y_log10(name='numt length')+
   ## scale_y_log10(name='numt count',breaks=c(0.1, 1, 10, 100, 1000,10000),labels=c('0','1','10','100','1000','10000'))+
    s.color + m.shape +
  theme_bw() +
  theme(axis.text.y = element_text(size=8,angle=90.0,hjust=0.5), axis.text.x=element_blank(),legend.position="none")+
  scale_x_discrete(name="")

pC <- ggplot() + labs(caption=plabel) + theme_void()

pgA <- plot_grid(NULL,pA.l,nrow=1,rel_widths=c(0.040,0.96))

plot_grid(pB, pgA, pC, nrow=3,rel_heights=c(0.4,0.6,0.1))
