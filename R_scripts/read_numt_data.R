
read.numt.data <- function (dir, file, fields) {

    print(paste0(dir,file))
    data <-read.table(paste0(dir,file),header=T,sep='\t',row.names=NULL)

#    print(colnames(data))
#    print(fields)
#    print(head(data))

    ## extracted shared set of data
    t.data <- data[,fields]

    ## create f_name from tag

    ## if (unique(data$algo)=='TFASTX') {
    ##    m_str <- paste0('_',unique(data$matrix),'m_tfx')
    ##    t.data$f_name<-sub(m_str,"",t.data$tag)
    ## }
    ## if (unique(tolower(data$algo))=='blastn') {
    ##    m_str <- '_blastn_blastnm_blastn'
    ##    t.data$f_name<-sub(m_str,"",t.data$tag)
    ## }
    ## else {
    ##    t.data$f_name=data$name
    ## }

    ## print(dim(t.data))
#    print(head(t.data))
    t.data
}
