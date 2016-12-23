library("ggplot2")
library("reshape2")
#library("Cairo")

args<-commandArgs(T)

d=read.table(args[1],header=T)
d=d[,-c(2,3)]
#lines_color=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#8B1A1A","#CC6666","#000000", "#9999CC", "#66CC99","#A020F0")
agcd<-melt(d,id.vars="Reads")
png(args[2])
#p=ggplot(d)+geom_line(aes(x=Reads,y=Pathogenic))+geom_line(aes(x=Reads,y=X12D))+geom_line(aes(x=Reads,y=X12V))+ geom_line(aes(x=Reads,y=X13D))
p=ggplot(agcd,aes(x=Reads,y=value,colour=variable))+geom_line()
#print(p)
#dev.off()
ggsave(file=args[2],p,width=4.5,height=4.5,type="cairo-png")
