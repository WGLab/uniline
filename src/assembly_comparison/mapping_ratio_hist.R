require('ggplot2')
args <- commandArgs(TRUE)
input <- as.character(args[1])
output <- paste(input,".pdf",sep="")
x = read.delim(input,as.is=TRUE)
x$Mapping_ratio = as.numeric(sub("%","",x$Mapping_ratio))
x$CatByLength = cut(breaks=c(0,1e4,1e5,6e5,1e6,1e8),x$Total_length)
pdf(output)
print(ggplot(x,aes(Mapping_ratio))+geom_histogram(aes(group = CatByLength, fill = CatByLength),position="dodge")+xlab("Mapping ratio(%)"))
print(ggplot(x,aes(Mapping_ratio))+geom_histogram(aes(group = CatByLength, fill = CatByLength),position="dodge")+xlab("Mapping ratio(%)")+xlim(c(0,0.99)))
print(ggplot(x,aes(Mapping_ratio))+geom_histogram(aes(group = CatByLength, fill = CatByLength),position="dodge")+xlab("Mapping ratio(%)")+xlim(c(0,0.75)))
print(ggplot(x,aes(Mapping_ratio))+geom_histogram(aes(group = CatByLength, fill = CatByLength),position="dodge")+xlab("Mapping ratio(%)")+xlim(c(0,0.5)))
dev.off()
