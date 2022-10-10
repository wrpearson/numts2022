
read.numt.hits <- function (dir, file, fields) {

##    print(paste0(dir,file))
    data <-read.table(paste0(dir,file),header=TRUE,sep='\t',row.names=NULL)

#    print(colnames(data))
#    print(fields)
#    print(head(data))

    ## extracted shared set of data
    t.data <- data[,fields]

##    print(head(t.data))
    t.data
}
