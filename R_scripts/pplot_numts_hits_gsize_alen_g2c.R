#!/usr/bin/env Rscript --vanilla

################
## plot_numts_summ_gsize_mat_g2c.R numts_summ_BLN-MD40_00_TFX_n.tab
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
library("ggplot2")
## library("ggrepel")
library("RColorBrewer")
library("patchwork")

p.name<-'pplot_numts_hits_gsize_alen_g2c.R'

args <- commandArgs(trailingOnly=TRUE)

plabel=paste(c(p.name,"\n",args),collapse=' ', sep=' ')

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

for (file in hit_files) {

    n_parts<-strsplit(file,'.', fixed=TRUE)
    t_name <- n_parts[[1]][1]

    file_parts <- strsplit(t_name,'_', fixed=TRUE)
    m_name <- file_parts[[1]][4]
    m_algo <- file_parts[[1]][6]
    aln_str <- file_parts [[1]][5]

    ## need to shift algo to matrix for BLN/BLNT

    print(paste("file=",file))
    print(paste(t_name, m_name, m_algo, aln_str))

    hit_data = read.numt.hits(data_dir, file, save.cols)

    ## print("head(t_data1)")
    ## print(head(t_data))

    hit_data[hit_data$common_name != 'mtGenome',]$alen = hit_data[hit_data$common_name != 'mtGenome',]$alen * 3

    ## t_data <- aggregate(cbind(h_cnt=qlen) ~ taxon_id, data=hit_data, FUN=NROW)
    ## wt_data <- aggregate(cbind(wh_cnt=alen) ~ taxon_id, data=hit_data, FUN=sum)

    t_data <- aggregate(cbind(med=alen) ~ taxon_id, data=hit_data, FUN=median)

    ## t_data = merge(t_data, wt_data[,c('taxon_id','wh_cnt')],by='taxon_id')

    if (m_algo == 'BLNT') {
        m_name <- paste0(m_name,'T')
    }

    t_data$matrix <- m_name
    t_data$algo <- m_algo
    t_data$l_aln <- aln_str

    ## print("head(t_data)")
    ## print(head(t_data))

    no_data <- subset(gt.data,!(taxon_id %in% t_data$taxon_id))
    ## head(no_data)

    if (NROW(no_data) > 0) {
        zero_data <- data.frame(taxon_id=no_data$taxon_id, med=0, matrix=m_name, algo=m_algo, l_aln=aln_str)
        print("head(zero_data)")
        print(head(zero_data))
        t_data <- rbind(t_data, zero_data)
    }

    summ_data <- rbind(summ_data, t_data)
}

## matrix symbols: 0:square 3:+, 4:x, 1:circle, 6: upside down triangle, 16:filled circle, 17:filled triangle, 5:diamond, 12:square w/+
##                 BLNG     BLN  BLNT MD10      (MD20)                    MD40             (BP62)              BLNGT       MD40BG

m_types = c('MD40','BLNGT', 'MD40BG')
m.shapes = c('MD40'=16, 'BLNGT'=5,'MD40BG'=12)

## unique(summ_data$matrix)

summ_data$matrix = factor(summ_data$matrix,levels=m_types,ordered=TRUE)

summ_data <- merge(summ_data, gt.data[,c('taxon_id','g_type','mtGenome_acc','s_name','Genome_MB','Class')], by='taxon_id')

summ_data$gsize_f = sprintf("%.1fG",summ_data$Genome_MB/1000)

## cat("head(summ_data):\n")
## head(summ_data)

summ_data$g_type = factor(summ_data$g_type, levels=g_types, labels=g_types, ordered=TRUE)
summ_data$matrix = factor(summ_data$matrix, levels=m_types, labels=m_types, ordered=TRUE)

## str(summ_data)

## colors for genome quality
br.gq_ref.colors <- brewer.pal(8,'GnBu')
gq.ref.colors <- br.gq_ref.colors[4:7]
br.gq_gb.colors <- brewer.pal(8,'OrRd')
gq.gb.colors <- br.gq_gb.colors[4:7]


## genome classes
c.colors <- brewer.pal(7,'Dark2')
class.map <- c('mammal'=c.colors[1],'marsupial'=c.colors[2],'monotreme'=c.colors[3],'reptile'=c.colors[4],'amphibian'=c.colors[5],'bird'=c.colors[6],'fish'=c.colors[7])
class.names <- names(class.map)
summ_data$Class = factor(summ_data$Class, levels=class.names, labels=class.names, ordered=TRUE)

cl.colors <- scale_color_manual(values=class.map,labels=class.names)

g.colors <- brewer.pal(4,'Dark2') # 'Dark2', 'Set2', 'Paired'
g.shapes <- c(16,16,15,17)

g.s.color0 <- scale_color_manual(values=g.colors,guide='none')
g.s.color <- scale_color_manual(values=g.colors,labels=c('RefSeq','GenBank'))
g.s.shape <- scale_shape_manual(values=g.shapes)

m.shape <- scale_shape_manual(values=m.shapes,labels=m_types,guide=guide_legend(reverse=TRUE))
m.shape0 <- scale_shape_manual(values=m.shapes,guide='none')


## cat("**tr$tip.label\n")
## str(tr)
## tr$tip.label

## cat("**dd\n")
## dd

length(unique(summ_data$s_name))

## figure out levels based on genome size
## summ_data$s_name <- factor(summ_data$s_name, levels=[order(tip.dd$y,decreasing=FALSE),]$s_name, ordered=TRUE)
## summ_data$s_name <- factor(summ_data$s_name, levels=get_taxa_name(pA.l), ordered=TRUE)

summ_data$s_name2 <- gsub("_",". ",summ_data$s_name)
## summ_data$s_name2 <- factor(summ_data$s_name2, levels=order(summ_data$Genome_MB),ordered=TRUE)

summ_data_o <- summ_data %>% arrange(Genome_MB) %>% mutate(s_name2o = factor(s_name2, levels=unique(s_name2), ordered=TRUE))

## cat("head/tail(summ_data_o)\n")
## head(summ_data_o)
## tail(summ_data_o)

u_sname2o = data.frame(s_name2o=unique(summ_data_o$s_name2o))

if (NROW(u_sname2o) %% 2 == 1) {
    u_sname2o$y_off = c(rep(c(0,1),(NROW(u_sname2o))/2),0)
} else {
    u_sname2o$y_off = c(rep(c(0,1),(NROW(u_sname2o))/2))
}

u_sname2o <- merge(u_sname2o, summ_data_o[,c('s_name2o','gsize_f')], by='s_name2o')

u_sname2o <- unique(u_sname2o)

## str(summ_data_o)

class_order = gt.data[order(gt.data$Genome_MB),]$Class

summ.d.class = class.map[class_order]
## summ.d.class

ymax=500

cat("low-data\n")
summ_data_o[summ_data_o$taxon_id %in% c(7955,8364,35670,31033),]

pB <-   ggplot(summ_data_o, aes(x=s_name2o,y=med)) +
   geom_jitter(aes(color=g_type,shape=matrix),size=2,alpha=0.8,width=0.2) +
   geom_text(data=u_sname2o, aes(x=s_name2o,y=ymax-25-25*y_off,label=gsize_f),size=2.5) +
   scale_y_continuous(name='length (nt)')+ g.s.color + m.shape +
  theme_bw() +
##   theme(axis.text.y = element_text(size=8,angle=90.0,hjust=0.5),axis.text.x = element_blank(),legend.title=element_blank()) +
  theme(axis.text.y = element_text(size=8,angle=90.0,hjust=0.5),axis.text.x = element_text(size=8,angle=45.0,vjust=1, hjust=1,color=summ.d.class),legend.title=element_blank())+
  ggtitle("median numt length by genome size") + scale_x_discrete(name="")

## cat("pB\n")
## ggplot_build(pB)

pC <- ggplot() + labs(caption=plabel) + theme_void()

if (Sys.getenv(c("PLOT_DOLAB"))=='none') {
    pdf(file=paste0(pdf_dir, t_name, "_gsize_alen_g2c_p_NL.pdf"),width=8, height=4.0)
    pB
} else {
    pdf(file=paste0(pdf_dir, t_name, "_gsize_alen_g2c_p.pdf"),width=8, height=4.2)
    pB / pC + plot_layout(heights = c(1,.025))
}

