# ks: this is a set of graphics functions designed for irt model

# input : data.panel() or a list() of data panels
# 			first column of the data.panel corresponds to the x-axis
#			remaining columns describe items  
#
# call: default call to the plot functions is: plotitems(mydatapanel)

# to do: a) validation methods 
# 		 b) color modules; now choose one of: colset = c("gray","black","red","blue")
# 		 c) write a short vignette for plotitems()
# -----------------------------------------------------------------
plotItemsGrd <- function(...){
	UseMethod("plotItemsGrd")
}
plotItemsGrd.default <- function(x,...){
	args <- list(...)
	itemplot(x,...)
}
plotItemsGrd.list <- function(x,...){
	itempanelplot(x,...)
}
# -----------------------------------------------------------------
itempanelplot <- function(datalist,ncol=2,newpage = TRUE,vp=NULL,coltype="gray",xscale=NULL,itemnm=NULL){
	
	nplots <- length(datalist)
	if(is.null(itemnm)){
		if( is.null(names(datalist)) )	itemnm <-  paste("item",1:nplots,sep=" ")
		else itemnm <-  names(datalist)
	}

	if(iseven(nplots)) nrow <- nplots%/%ncol
	else nrow <- nplots%/%ncol + 1

	if(nplots==1) ncol <- 1 
	
	if (newpage) grid.newpage()	# new page
 	if (!is.null(vp)) pushViewport(vp) # view port
#	pushViewport(viewport(clip="on"))
	temp.vp <- viewport(w=0.7, h=.7,layout = grid.layout(nrow, ncol))
	pushViewport(temp.vp)

  for (i in 1:nplots) {
    col <- (i - 1)%%ncol + 1
    row <- (i - 1)%/%ncol + 1
    # We go to each entry of the table
    panel.vp <- viewport(layout.pos.row = row, layout.pos.col = col)
	itemcellplot(data=datalist[[i]], 
               axis.left = (col == 1), 
               axis.left.label = isodd(row),
               axis.right = (col == ncol || i == nplots), 
               axis.right.label = iseven(row),
               axis.bottom = (row == nrow), 
               axis.bottom.label = isodd(col),
               axis.top = (row == 1), 
               axis.top.label = iseven(col),
			   coltype =coltype,
               vp = panel.vp,
			   itemnm=itemnm[i],
			   xscale=xscale,
			   newpage=FALSE) 
	}

	popViewport()
	if (!is.null(vp)) popViewport()
}
# -----------------------------------------------------------------
itemplot <- function(data,xscale=NULL,itemnm="item",newpage = FALSE,vp=NULL,coltype="gray"){
	if (newpage) grid.newpage()	# new page
 	if (!is.null(vp)) pushViewport(vp) # view port
	
	itemcellplot(data=data, 
               axis.left = TRUE, 
               axis.left.label = TRUE,
               axis.bottom = TRUE, 
               axis.bottom.label = TRUE,
			   coltype =coltype,
               vp = viewport(w=0.7, h=.7),
			   itemnm=itemnm,
			   xscale=xscale,
			   newpage=FALSE) 
  
  if (!is.null(vp))  popViewport()
}
# -----------------------------------------------------------------
itemcellplot <- function (
		  data=data.frame(x=runif(10),y=runif(10)), 		  
		  itemnm="item",
          xscale = NULL, 
          axis.left = TRUE,   axis.left.label = TRUE,
          axis.right = FALSE, axis.right.label = TRUE, 
          axis.bottom = TRUE, axis.bottom.label = TRUE, 
          axis.top = FALSE,   axis.top.label = TRUE,
          vp = NULL, coltype = "gray",
		  newpage = TRUE
		)
{
#	newpage <- TRUE
	if (newpage) grid.newpage()	# new page
	cltp <- setcolplot(coltype)	
	
	# xy-axis settings; default is range(x); otherwise input xscale with 5 x-ticks and 5 grill lines (equaly spaced)
	yscale = c(0,1) 
#	xscale = c(0,1)
	if(is.null(xscale)) xscale <- range(data[[1]]) 
	scaleprt <- seq(xscale[1],xscale[2],by=diff(xscale)/5)
	scaleprt <- round(scaleprt,2)	
	if (!is.null(vp)) pushViewport(vp)

	# divide a Viewport into Up and Bottom: 
	vpMain <- viewport(layout = grid.layout(2, 1,heights=unit(c(1,1),c("lines","null") )) )
	pushViewport(vpMain)
	
	# Upper viewport; title panel + (conditional axis.top=TRUE) x-axis 
	vpUp <- viewport(layout.pos.col=1, layout.pos.row=1,xscale=xscale  )
	pushViewport(vpUp)
	grid.rect(gp=gpar(fill=cltp$rectUpfill,lwd=1,col=cltp$rectUpcol))
	grid.text(itemnm,gp=gpar(col=cltp$textUpcol,fontface="bold",cex=0.8)) # whitesmoke

	# axis font setting
	axfont <- gpar(lwd=0.6,fontsize=9,lineend="round",col=cltp$rectUpcol,lineheight=0.7)	
  	if (axis.top) grid.xaxis(at=scaleprt,main = FALSE, label = axis.top.label,gp=axfont)
	popViewport() # Up

  # Bottom viewport: plot + (conditional axis.*=TRUE) xy-axis
	vpBottom <- viewport(layout.pos.col=1, layout.pos.row=2,xscale=xscale,yscale=yscale)
	pushViewport(vpBottom)
	grid.rect(gp=gpar(fill=cltp$rectBtmfill,lwd=1,col=cltp$rectBtmcol)) #snow

	if (axis.left) grid.yaxis(label = axis.left.label,gp=axfont)
	if (axis.right) grid.yaxis(main = FALSE, label = axis.right.label,gp=axfont)
	if (axis.bottom) grid.xaxis(at=scaleprt,label = axis.bottom.label,gp=axfont)
	pushViewport(dataViewport(name="plotRegion",xscale=xscale,yscale=yscale,clip="on"))
#	grid.grill(h = unit(seq(0.2, 0.8, 0.2), "npc"), gp=gpar(col = "grey",lwd=1))

	pushViewport(dataViewport(xscale=xscale,yscale=yscale,clip="on"))
 	grid.grill(h = unit(seq(0.2, 0.8, 0.2), "npc"),v=unit(scaleprt,"native"),default.units = "native", gp=gpar(col = "grey80",lwd=0.3))

  # Plot data
	nc <- ncol(data); colset <- cltp$itemcurves
	for(i in 2:nc){
		switch(colset,
				gray = cl <- paste("gray",1+(i-2)*15,sep=""),
				black = cl <- "black",
				red = cl <- paste("red",(i-1),sep=""),
				orange = ifelse(nc > 5, cl <- "orange3",cl <- paste("orange",(i-1),sep="") ),			
				blue = cl <- "deepskyblue3")
		grid.lines(x=unit(data[,1],"native"),y=unit(data[,i],"native"),gp=gpar(lwd=1.1,col=cl))	
	}
 
	popViewport(4) # close 3 views: Data, Buttom and Main

	if (!is.null(vp)) popViewport()
invisible(list(vpUp = vpUp, vpBottom = vpBottom))
}
# -----------------------------------------------------------------
# color settings
setcolplot <- function(colset=c("gray","black","red","blue")){
	gr <- list(
		rectUpfill= "lightskyblue3",
		rectUpcol = "lightskyblue4",
		textUpcol = "whitesmoke",
		rectBtmfill = "grey98",
		rectBtmcol = "lightskyblue4",
		itemcurves = "gray"
	)
	bl <- list(
		rectUpfill= "deepskyblue",
		rectUpcol = "deepskyblue2",
		textUpcol = "white",
		rectBtmfill = "grey98",
		rectBtmcol = "deepskyblue2",
		itemcurves = "blue"
	)
	bw <- list(
		rectUpfill= "gray10",
		rectUpcol = "gray1",
		textUpcol = "whitesmoke",
		rectBtmfill = "grey98",
		rectBtmcol = "gray1",
		itemcurves = "black"
	)
	rd <- list(
		rectUpfill= "salmon4",
		rectUpcol = "salmon3",
		textUpcol = "white",
		rectBtmfill = "grey98",
		rectBtmcol = "salmon4",
		itemcurves = "orange"
	)
	type=match.arg(colset)
	switch(type,
		gray = ans <- gr,
		black = ans <- bw,
		red = ans <- rd,
		blue = ans <- bl
		)
return(ans)	
}
# -----------------------------------------------------------------
# --------------------------------------------------------------
# this plot is is not attach to the plotitems class. It serves as an info plot
# during optimization
plotItemsGrdInfo <- function(Ftab,rnames,plotInfoSet,xscale=c(-2.5,2.5),showinfo=FALSE){
	assignList(plotInfoSet)
	lnInfo <- length(plotInfoSet)
	nmInfo <- names(plotInfoSet)

	timenow <- proc.time()[3]
	if( ( (timenow-plotInfoSet$starttime)%%10 )< 0.3){
		if(showinfo){
			pushViewport(viewport(h=0.95,w=0.95))
			vpMain <- viewport(layout = grid.layout(2, 1,heights=unit(c(1+lnInfo,1),c("lines","null") )) )
			pushViewport(vpMain)
			vpUp <- viewport(layout.pos.col=1, layout.pos.row=1)
			pushViewport(vpUp)
			grid.rect(gp=gpar(fill="gray99",lwd=.5,col="gray90"))

			gpn <- gpar(col="black",fontface="bold",just="left",cex=0.8)
			gpninfo <- gpar(col="red4",cex=0.8)
			for(i in 2:lnInfo) grid.text(paste(nmInfo[i],"=",plotInfoSet[[i]]),just="left",x=unit(3, "lines"),y=unit(i, "lines"),gp=gpninfo) # whitesmoke
			grid.text("IRT Info panel",y=unit(lnInfo+1, "lines"),gp=gpn) # whitesmoke

			popViewport() # Up

		  # Bottom viewport: plot + (conditional axis.*=TRUE) xy-axis
			vpBottom <- viewport(layout.pos.col=1, layout.pos.row=2)
			pushViewport(vpBottom)
			plotItemsGrd(Ftab,itemnm=rnames,xscale=xscale,coltype="black",newpage=FALSE)
			popViewport(3)
		}else{
			plotItemsGrd(Ftab,itemnm=rnames,xscale=xscale,coltype="black",newpage=FALSE)			
		}
	}
}

