#!/usr/bin/env Rscript --vanilla

################
## plot_numts_order1rs.R -- plot boxplots for each mitogene with lines connecting species
##
##
################

library("dplyr")
library("ggplot2")
library("patchwork")
library("RColorBrewer")

## this plot orders from left to right based on gene size.
## the rs version scales the counts by the gene length

args <- commandArgs(trailingOnly=TRUE)

## this plot orders from left to right based on gene size.
## the rs version scales the counts by the gene length

p.name<-'pplot_numts_hits_order2rs.R'
plabel=paste(c(p.name,args),collapse=' ', sep=' ')

data_dir='data/'
pdf_dir='figs/'

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

hit_files

if (!exists("read.numt.hits",mode='function')) source("./R_scripts/read_numt_hits.R")
save.cols=c('search_id','tag','canon_acc','common_name','qlen','s_acc','percid','alen','expect','taxon_id')

summ_data = NULL

## read in the lengths of the human proteins
hum_prots <- read.table('human_mito_prot.tab',header=TRUE, sep='\t')

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

    t_data <- aggregate(cbind(h_cnt=qlen) ~ taxon_id + common_name, data=hit_data, FUN=NROW)

    if (m_algo == 'BLNT') {
        m_name <- paste0(m_name,'T')
    }

    t_data$matrix <- m_name
    t_data$algo <- m_algo
    t_data$l_aln <- aln_str

    no_data <- subset(gt.data, !(taxon_id %in% t_data$taxon_id))

    if (NROW(no_data) > 0) {
        zero_data <- data.frame(taxon_id=no_data$taxon_id, h_cnt=0, matrix=m_name, algo=m_algo, l_aln=aln_str, common_name='XXXX')
        print("zero_data")
        print(zero_data)
        t_data <- rbind(t_data, zero_data)
    }

    ## print("head(t_data)")
    ## print(head(t_data))

    summ_data <- rbind(summ_data, t_data)
}

## head(summ_data)
## summary(summ_data)

if (nrow(summ_data[summ_data$h_cnt < 1,])>0) {
    summ_data[summ_data$h_cnt < 1,]$h_cnt = 0.1
}

plot.title=sprintf("(matrix: %s)",m_name)

## is it refseq vs genbank

g_types = c('refseq','genbank','H. sapiens','M. muntjak')
g_types = factor(g_types, levels=g_types, ordered=TRUE)

summ_data <- merge(summ_data, gt.data[,c('taxon_id','g_type','s_name')], by='taxon_id')

print(summ_data[summ_data$s_name=='H_sapiens',])

summ_data$g_type = factor(summ_data$g_type, levels=g_types, labels=g_types, ordered=TRUE)
summ_data[summ_data$s_name=='H_sapiens',]$g_type=g_types[3]
summ_data[summ_data$s_name=='M_muntjak',]$g_type=g_types[4]

hum_prots_srt <- hum_prots[order(hum_prots$qlen),]

hum_ave = mean(hum_prots$qlen)
hum_prots$q_scale = hum_prots$qlen/hum_ave
## hum_prots

## head(summ_data)

summ_data <- merge(summ_data, hum_prots[,c('common_name','q_scale')],by='common_name')
summ_data$h_cnt_s <- summ_data$h_cnt / summ_data$q_scale

## head(summ_data)

summ_data$common_name <- factor(summ_data$common_name,levels=hum_prots_srt$common_name, ordered=TRUE)
## summ_data$common_name <- factor(summ_data$common_name, ordered=TRUE)

hc_meds <- aggregate(h_cnt ~ common_name + g_type, data=summ_data, median)
## print("hc_meds: ")
## hc_meds

hc_s_meds <- aggregate(h_cnt_s ~ common_name + g_type, data=summ_data, median)
## print("hc_s_meds: ")
## hc_s_meds

hc_meds2 = reshape(hc_meds, idvar="common_name", timevar='g_type',direction='wide')
hc_meds2 <- merge(hc_meds2, hum_prots[,c('common_name','qlen')],by='common_name')
## hc_meds2

hc_s_meds2 = reshape(hc_s_meds, idvar="common_name", timevar='g_type',direction='wide')


## hc_meds2$h_cnt.refseq
## hc_meds2$h_cnt.genbank

hc_meds2$display = paste0(hc_meds2$h_cnt.refseq,'/', hc_meds2$h_cnt.genbank)
## hc_meds2

hc_s_meds2$display = sprintf("%.0f/%.0f",hc_s_meds2$h_cnt_s.refseq,hc_s_meds2$h_cnt_s.genbank)
## hc_s_meds2

## afterwards, it should be the factor  (not used for order2
## summ_data[summ_data$s_name=='H_sapiens',]$g_type=g_types[3]
## summ_data[summ_data$s_name=='M_muntjak',]$g_type=g_types[4]

## print(summ_data[summ_data$s_name=='H_sapiens',])
## print(summ_data[is.na(summ_data$g_type),])
## summary(summ_data)

g.colors <- brewer.pal(4,'Dark2') # 'Dark2', 'Set2', 'Paired'
g.shapes <- c(16,16,15,2)
s.color  <- scale_color_manual(values=g.colors)
sb.color <- scale_color_manual(values=g.colors,guide='none')
s.shape  <- scale_shape_manual(values=g.shapes)
sb.shape <- scale_shape_manual(values=g.shapes,guide='none')
s.fill <- scale_fill_manual(values=alpha(g.colors,0.00), labels=g_types)

y.limits = c(0.1,10000)

## uncomment after changing g_types
fivenum(summ_data$h_cnt)
gene_li <- with(summ_data, tapply(h_cnt, common_name, fivenum))
gene_fr <-  data.frame(matrix(unlist(gene_li),ncol=5, byrow=TRUE))
colnames(gene_fr) = c('min','q1','med','q3','max')
gene_fr$common_name <- rownames(gene_li)

## uncomment after changing g_types
## fivenum(summ_data$h_cnt_s)
gene_li_s <- with(summ_data, tapply(h_cnt_s, common_name, fivenum))
gene_fr_s <-  data.frame(matrix(unlist(gene_li_s),ncol=5, byrow=TRUE))
colnames(gene_fr_s) = c('min','q1','med','q3','max')
gene_fr_s$common_name <- rownames(gene_li_s)

## head(gene_fr)
## rownames(gene_fr)=rownames(gene_li)
## head(gene_fr)

## print(paste(t_name, "by gene summary"))
## summary(gene_fr)
print("quantiles of h_cnt medians")
quantile(gene_fr$med,probs=seq(0,1,0.1))

print("quantiles of h_cnts medians")
quantile(gene_fr_s$med,probs=seq(0,1,0.1))

taxon_li <- with(summ_data,tapply(h_cnt, s_name,      fivenum))
taxon_fr <- data.frame(matrix(unlist(taxon_li),ncol=5, byrow=TRUE),row.names=rownames(taxon_li))
colnames(taxon_fr) = c('min','q1','med','q3','max')
## print("h_cnt by taxon summary")
## summary(taxon_fr)
print("quantiles of medians")
quantile(taxon_fr$med,probs=seq(0,1,0.1))
## taxon_fr
## with(summ_data,tapply(h_cnt, s_name,      fivenum))

taxon_li_s <- with(summ_data,tapply(h_cnt_s, s_name, fivenum))
taxon_fr_s <- data.frame(matrix(unlist(taxon_li_s),ncol=5, byrow=TRUE),row.names=rownames(taxon_li_s))
colnames(taxon_fr_s) = c('min','q1','med','q3','max')
## print("h_cnt_s by taxon summary")
## summary(taxon_fr_s)
print("h_cnt_s quantiles of medians")
quantile(taxon_fr_s$med,probs=seq(0,1,0.1))
## taxon_fr
## with(summ_data,tapply(h_cnt, s_name,      fivenum))

## summ_data_s <-summ_data %>% arrange(common_name, s_name)

pA <- ggplot(summ_data,aes(x=common_name,y=h_cnt,fill=g_type)) +
   geom_boxplot(outlier.shape=NA)+
   geom_line(aes(color=g_type ),alpha=0.4) + ## ,position=position_jitterdodge(seed=123,jitter.width=0.8,jitter.height=0.01)) +
   geom_point(aes(shape=g_type, color=g_type),size=2,alpha=0.4) + ##  position=position_jitterdodge(seed=123,jitter.width=0.8,jitter.height=0.01)) +
    ##   geom_text(data=hc_meds,aes(x=common_name, y = h_cnt*1.1, label=h_cnt), vjust= 0, size=3) +
       geom_text(data=hc_meds2,aes(x=common_name, y = 0.7*y.limits[2] , label=display), vjust= 0, size=4) +
       geom_text(data=hc_meds2,aes(x=common_name, y = 0.15 , label=qlen), vjust= 0, size=4) +
  scale_y_log10(name='numt count',limits=y.limits,breaks=c(0.1,1,10,100,1000,10000), labels=c('0','1','10','100','1000','10000')) + s.color + s.shape +
  theme_bw() +
  theme(legend.title=element_blank(), legend.text=element_text(size=12), axis.text=element_text(size=12),axis.title=element_text(size=14)) +
  annotate("text", x = 14.25, y = 0.7 * y.limits[2], label = 'median:\nrefseq/genbank',vjust=0.5, hjust=0) +
  annotate("text", x = 14.25, y = 0.15, label = 'gene length (aa)',vjust=0, hjust=0) +
  coord_cartesian(xlim = c(0.75, 13.25), clip = "off") +
  scale_x_discrete(name='mitochondrial gene') + ggtitle(paste("A. numts by gene",plot.title))

pB <- ggplot(summ_data,aes(x=common_name,y=h_cnt_s,fill=g_type)) +
   geom_boxplot(outlier.shape=NA)+
   geom_line(aes(group=interaction(s_name,g_type), color=g_type),alpha=0.4) + ## ,position=position_jitterdodge(seed=123,jitter.width=0.8,jitter.height=0.01)) +
   geom_point(aes(group=interaction(s_name,g_type), color=g_type),size=2,alpha=0.4) + ##  position=position_jitterdodge(seed=123,jitter.width=0.8,jitter.height=0.01)) +
    ##   geom_text(data=hc_meds,aes(x=common_name, y = h_cnt*1.1, label=h_cnt), vjust= 0, size=3) +
       geom_text(data=hc_s_meds2,aes(x=common_name, y = 0.7*y.limits[2] , label=display), vjust= 0, size=4) +
  scale_y_log10(name='scaled numt count',limits=y.limits,breaks=c(0.1,1,10,100,1000,10000), labels=c('0','1','10','100','1000','10000')) + s.color + s.shape +
  theme_bw() +
  theme(legend.title=element_blank(), legend.position='none', axis.text=element_text(size=12),axis.title=element_text(size=14)) +
  annotate("text", x = 14.25, y = 0.7 * y.limits[2], label = 'median:\nrefseq/genbank',vjust=0.5, hjust=0) +
  coord_cartesian(xlim = c(0.75, 13.25), clip = "off") +
  scale_x_discrete(name='mitochondrial gene') + ggtitle(paste("B. numts scaled by average gene length",plot.title))

pC <- ggplot() + labs(caption=plabel) + theme_void()

## pdf(file=paste0(pdf_dir, t_name, "_ord1rs.pdf"),width=10, height=7)

if (Sys.getenv(c("PLOT_DOLAB"))=='none') {
    pdf(file=paste0(pdf_dir, t_name, "_ord2rs_p_NL.pdf"),width=10, height=7)
    pA / pB
} else {
    pdf(file=paste0(pdf_dir, t_name, "_ord2rs_p.pdf"),width=10, height=7.3)
    pA / pB / pC + plot_layout(heights=c(1,1,0.025))
}


