library("ggplot2")
library("Cairo")
#library("gridExtra")


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
	library(grid)
	# Make a list from the ... arguments and plotlist
	plots <- c(list(...), plotlist)
	numPlots = length(plots)
	# If layout is NULL, then use 'cols' to determine layout
	if (is.null(layout)) {
	# Make the panel
	# ncol: Number of columns of plots
	# nrow: Number of rows needed, calculated from # of cols
		layout <- matrix(seq(1, cols * ceiling(numPlots/cols*1.5)),	ncol = cols, nrow = ceiling(numPlots/cols))
	}

	if (numPlots==1) {
		print(plots[[1]])
	} else {
	# Set up the page
		grid.newpage()
		pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

		# Make each plot, in the correct location
		for (i in 1:numPlots) {
		# Get the i,j matrix positions of the regions that contain this subplot
			matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
			print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
			layout.pos.col = matchidx$col))
		}
	}
}

d=read.table("../input.xls",header=T)
rownames(d)=d[,2]
dd=d[,-c(1,3:20)]
label_for_barplot=names(dd)[2:14]
label_for_piechart=names(dd)[3:13]
color_for_bar=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#8B1A1A","#CC6666","#000000", "#9999CC", "#66CC99","#A020F0")
color_for_pie=c("#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#8B1A1A","#CC6666","#000000", "#9999CC", "#66CC99","#A020F0")
#color_for_barplot=c()
#color_for_piechart=c()
for(i in 1:nrow(dd))
{
	file=paste(dd[i,1],".png",sep="")
	png(file,width=1000,height=400)
	value_for_barplot=c(dd[i,2])

	value_for_piechart=c()
	for (j in 3:13){
		value_for_barplot = c(value_for_barplot,dd[i,j])
		value_for_piechart = c(value_for_piechart,dd[i,j])	
	}
	value_for_barplot = c(value_for_barplot,dd[i,14])
	value_for_barplot_final = data.frame(site=label_for_barplot,ratio=value_for_barplot)
	#draw barplot
	#cat(value_for_barplot_final$ratio)
	barplot = ggplot(value_for_barplot_final) + geom_bar(aes(x=site,y=ratio,fill=site),stat="identity") + xlab("") + scale_fill_manual(values=color_for_bar) + theme(axis.text.x=element_text(angle=270,size=10))

	#raw piechart
	percent_str <- paste(round(value_for_piechart/sum(value_for_piechart) * 100,2), "%", sep="")
	vfc=data.frame(Percentage=round(value_for_piechart/sum(value_for_piechart)*100,2), Type=label_for_piechart,percent=percent_str)
	#pie + scale_fill_manual(values = colours,labels = labels)
	piechart <- ggplot(vfc) +  geom_bar(aes(x=factor(1),y=Percentage,fill=Type),stat="identity")+coord_polar(theta="y") + xlab('Percentage') + ylab('')+ scale_fill_manual(values=color_for_pie)
	#grid.arrange(barplot, piechart, ncol=2, nrow=1, widths=c(2,1), heights=1)
	#print(piechart)
	#print(barplot)
	multiplot(barplot, piechart, cols=2)
	dev.off()
}
