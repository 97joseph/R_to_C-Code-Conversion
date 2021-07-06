# basic error functions for HA
# also routines to plot range charts, etc

ColumnRangeExtension<-function(y){
	# compute the range extension implied by a column of data, sorted in time
	lp=which(y>0)
	lpv=lp[1]:lp[length(lp)]
	sum(y[lpv]==0,na.rm=TRUE)
}

AshRangeExtension<-function(y){
	lp=which(y>0)
	lpv=lp[1]:lp[length(lp)]
	sum(y[lpv]==0,na.rm=TRUE)+sum(is.na(y[lpv]))
}
AshRangeExtensionStrong<-function(y){
	lp=which(y>0)
	lpv=lp[length(lp)]-lp[1]
}

# this is the passing error function, y1 is the column of the lower object, y2 is the column of the upper object
# y1type =0, appearence, 1-FAD, 2-LAD,  y2type is same info

PassingError<-function(y1,y1type,y2,y2type)
{
	if(y1type==0)
	{
		y1pos=which(y1==1)
	}
	else
	{
		if(y1type==1)
		{
			y1pos=min(which(y1==1),na.rm=TRUE)
		}
		else
		{
			y1pos=max(which(y1==1),na.rm=TRUE)
		}
	}
	if(y2type==0)
	{
		y2pos=which(y2==1)
	}
	else
	{
		if(y2type==1)
		{
			y2pos=min(which(y2==1),na.rm=TRUE)
		}
		else
		{
			y2pos=max(which(y2==1),na.rm=TRUE)
		}
	}
	
	if(y1pos>y2pos)
	{
		return(y1pos-y2pos)
	}
	else
	{
		return(0)
	}
	
}

PassingError2<-function(y1,y1type,y2,y2type)
{
	if(y1type==0)						# appearance, find list of all
	{
		y1pos=which(y1==1)
	}
	else
	{
		if(y1type==1)									#  Fad, find the lowest
		{
			y1pos=min(which(y1==1),na.rm=TRUE)
		}
		else											# Lad, find the highest
		{
			y1pos=max(which(y1==1),na.rm=TRUE)
		}
	}
	if(y2type==0)
	{
		y2pos=which(y2==1)
	}
	else
	{
		if(y2type==1)
		{
			y2pos=min(which(y2==1),na.rm=TRUE)
		}
		else
		{
			y2pos=max(which(y2==1),na.rm=TRUE)
		}
	}
	
	if(y1pos>y2pos)															# error count when y1pos is higher than y2pos
	{
		temp=sum((which(y1==1)-y2pos)*(which(y1==1)>y2pos))
		temp=temp+sum(-1*(y1pos-which(y2==1))*(which(y2==1)<y2pos))
		return(temp)
	}
	else
	{
		return(0)
	}
	
}




# Net Range Extension calculation
NetRangeExtension<-function(y){
	
	sum(apply(y,2,ColumnRangeExtension))
	#rx=parallel(sum(apply(y,2,ColumnRangeExtension)))
	#collect(rx)
	#return(rx)
}

# Net Range Extension calculation
NetRangeExtension2<-function(y){
	
	#sum(apply(y,2,ColumnRangeExtension))
	rx<-parallel(sum(apply(y,2,ColumnRangeExtension)))
	a=collect(rx)
	b=a[[1]]
	return(b)
}




# isotope error function, fdistance based on y values only!
fdistance<-function(x)
{
	# column 2 is the y value at the x value in column 1
	# this is a distance based slotting error calculation
	#nd=sum((diff(x[,2])^2+diff(x[,1])^2)^0.5)
	nd=sum((diff(x[,2])^2)^0.5)

		return(nd)
}

# isotope error function, based on transitions

ftransitions<-function(x)
{
	# column 2 is the y value at the x value in column 1
	y=diff(x[,2])
	y=sign(y)
	yalt=y[-1]-y[-length(y)]
	reversals=sum(yalt!=0)
	return(reversals)
}

ftransitions2<-function(x)
{
	# column 2 is the y value at the x value in column 1
	# assumes a single column of data
	y=diff(x[!is.na(x)])
	y=sign(y)
	yalt=diff(y)
	reversals=sum(abs(yalt)/2)
	return(reversals)
}


# plot a range chart, also calculates FAD/LAD listing in column format

PlotRangeChart<-function(j,doplot=TRUE){
	# j$pen, j$initpen,j$d, j$history
	pt=dim(j$d)
	print("size of data set")
	print(pt)
	trange=apply(j$d[,3:pt[2]],2,ColumnFADLAD)
	trange=matrix(unlist(trange),ncol=2,byrow=TRUE)
    y=1:(pt[2]-2)
	y=cbind(y,trange)
	if(doplot){
		quartz()
		plot(c(0,(pt[2]-2)),c(0,(pt[1]+1)),type="n",xlab="Taxa Number",ylab="Position in Composite")
		points(y[,1],y[,2],col="black")
		points(y[,1],y[,3],col="black")
		for(k in 1:(pt[2]-3)){
		   lines(c(y[k,1],y[k,1]),c(y[k,2],y[k,3]))	
		}
	}
	
    return(y)	
}


#------------------------------------------------------------------------
# plot a range chart, also calculates FAD/LAD listing in column format


PlotRangeChartPS2<-function(j_ch,species_list,doplot=TRUE,overplot=FALSE,hbar=0.66,tlabels=FALSE,barOrder="None",MarkSections=TRUE){
	# j$pen, j$initpen,j$d, j$history
	# if multiple plotting- alter hbar, try hbar=0.4
	# note that this version assumes taxa start in column 5,  so score is column 1, section is column 2, height=3	
	pt=dim(j_ch)
	print("size of data set")
	print(pt)
	# get taxa range
	toffset=4;				# offset for start of taxa, Horizon anneal 3 uses offset of 4
	trange=apply(j_ch[,species_list+toffset],2,ColumnFADLAD)
	#print(trange)
	trange=matrix(unlist(trange),ncol=2,byrow=TRUE)
    y=species_list
    npts=length(species_list)
	y=cbind(y,trange)
	yoff=min(y[,2])
	#y[,2]=y[,2]-yoff
	#y[,3]=y[,3]-yoff
	
	#labelactive=(class(tlabels)=="data.frame")
	cat(tlabels)
	tlabels_val=colnames(j_ch[,4+species_list])
	
	if(barOrder=='Fad')
	{
		forder=order(y[,2])
		y=y[forder,]
		tlabels_val=tlabels_val[forder]
	}
	if(barOrder=="Lad")
	{
	  forder=order(y[,3])
	  y=y[forder,]
	  tlabels_val=tlabels_val[forder]
	}  
	

	
	
	if(doplot){
		if(overplot==FALSE)
		{
	#	quartz()
		plot(c(0,npts+1),c(-(pt[1])/5,(pt[1]+1)),type="n",xlab="Taxa Number",ylab="Position in Composite")
		  for(k in 1:(npts)){
		  	 lines(c(k,k+hbar),c(y[k,2],y[k,2]),col="black")
	         lines(c(k,k+hbar),c(y[k,3],y[k,3]),col="black")
		     lines(c(k,k),c(y[k,2],y[k,3]))	
		     lines(c(k,k)+hbar,c(y[k,2],y[k,3])) 
		     if(tlabels)
		     {
		     	text(k+0.5,y[k,2],c(tlabels_val[k],"   "),adj=c(1,0),cex=0.6,srt=90,font=3) 
		     	#cat(tlabels$TaxaName[y[k,1]],"\n")
		     }
		  }
		}
		else
		{	
		  hoffset=0.45
		  
		  for(k in 1:(npts)){
		     lines(c(y[k,1],y[k,1]+hbar)+hoffset,c(y[k,2],y[k,2]),col="black",lty=2)
	         lines(c(y[k,1],y[k,1]+hbar)+hoffset,c(y[k,3],y[k,3]),col="black",lty=2)
		     lines(c(y[k,1],y[k,1])+hoffset,c(y[k,2],y[k,3]),col="black",lty=2)	
		     lines(c(y[k,1],y[k,1])+hbar+hoffset,c(y[k,2],y[k,3]),col="black",lty=2)	
		  }
	
		}
		
	}
	
	# now go through the taxa, and plot the sections they were found in
	if(MarkSections)
	{  
	for(k in 1:npts)
	{
	    myfinds=which(j_ch[,toffset+y[k,1]]==1)
	    myfinds=myfinds[!is.na(myfinds)]
	    mySecs=j_ch[myfinds,2]
	    for(k2 in 1:length(myfinds))
	    {
	      mytext=sprintf("%i",mySecs[k2])
	      text(k+0.5,myfinds[k2],mytext,adj=c(1,0),cex=0.5,col='red')
	    }
	  }
  }
	
	
	# use labels from column names in j_ch
	
	if(tlabels)
	{
		y2=data.frame(c9code=y[,1],FAD=y[,2],LAD=y[,3],SpeciesName=tlabels_val)
	}
	else
	{
		y2=data.frame(c9code=y[,1],FAD=y[,2],LAD=y[,3],SpeciesName=rep(NA,length(species_list)))
	}
    return(y2)	
}





# plot a range chart based off a C9 style input data set

PlotRangeChartC9<-function(c9sol,species_list,doplot=TRUE,overplot=FALSE,hbar=0.66,tlabels=c("na"),barOrder="None"){
	# c9sol- solution in c9 format
	maxtaxa=max(c9sol[,1])
	sort_ord=order(c9sol[,1]+(maxtaxa+10)*c9sol[,2])
	c9sol=c9sol[sort_ord,]
	c9length=dim(c9sol)[1]/2
	y=cbind(c9sol[1:c9length,1],c9sol[c9sol[,2]==1,3],c9sol[c9sol[,2]==2,3])
	y=y[species_list,]
	npts=length(species_list)
	if(barOrder=='Fad')
	{
		forder=order(y[,2])
		y=y[forder,]
	}
	pt=dim(y)
	labelactive=(class(tlabels)=="data.frame")
	if(doplot){
		if(overplot==FALSE)
		{
	#	quartz()
		plot(c(0,npts+1),c(-(max(c9sol[,3]))/5,(max(c9sol[,3])+1)),type="n",xlab="Taxa Number",ylab="Position in Composite")
		  for(k in 1:(npts)){
		  	 lines(c(k,k+hbar),c(y[k,2],y[k,2]),col="black")
	         lines(c(k,k+hbar),c(y[k,3],y[k,3]),col="black")
		     lines(c(k,k),c(y[k,2],y[k,3]))	
		     lines(c(k,k)+hbar,c(y[k,2],y[k,3])) 
		     if(labelactive)
		     {
		     	text(k+0.5,y[k,2],c(tlabels$TaxaName[y[k,1]],"   "),adj=c(1,0),cex=0.6,srt=90,font=3) 
		     }
		  }
		}
		else
		{	
		  hoffset=0.45
		  
		  for(k in 1:(npts)){
		     lines(c(y[k,1],y[k,1]+hbar)+hoffset,c(y[k,2],y[k,2]),col="black",lty=2)
	         lines(c(y[k,1],y[k,1]+hbar)+hoffset,c(y[k,3],y[k,3]),col="black",lty=2)
		     lines(c(y[k,1],y[k,1])+hoffset,c(y[k,2],y[k,3]),col="black",lty=2)	
		     lines(c(y[k,1],y[k,1])+hbar+hoffset,c(y[k,2],y[k,3]),col="black",lty=2)	
		  }
	
		}
		
	}
	if(labelactive)
	{
		y=cbind(y,tlabels$TaxaName[y[,1]])
		colnames(y)<-c("c9Code","FAD","LAD","SpeciesName")
	}
    return(y)	
}

# plot two c9 solutions for comparison
ComparisonPlotC9<-function(c9sol1,c9sol2,species_list,doplot=TRUE,overplot=FALSE,hbar=0.66,tlabels=c("na"),addcounts=FALSE,taxacounts='na')
{
	#addcounts=TRUE, plots the number of sections taxa are seen in with the label
	#taxacounts should be the output of InformativeTaxa()- second column is the site counts
	
    maxtaxa=max(c9sol1[,1])
    sort_ord=order(c9sol1[,3])
	c9sol1=c9sol1[sort_ord,]
	c9sol1[,3]=1:dim(c9sol1)[1]                        # ordinal numbering
	sort_ord=order(c9sol1[,1]+(maxtaxa+10)*c9sol1[,2])
	c9sol1=c9sol1[sort_ord,]							#ordered by taxa number
	c9length1=dim(c9sol1)[1]/2
	
	sort_ord=order(c9sol2[,3])
	c9sol2=c9sol2[sort_ord,]
	c9sol2[,3]=1:dim(c9sol2)[1]								#ordinal numbering
	sort_ord=order(c9sol2[,1]+(maxtaxa+10)*c9sol2[,2])		# ordered by taxa number
	c9sol2=c9sol2[sort_ord,]
	c9length2=dim(c9sol2)[1]/2
	if(c9length1!=c9length2)
	{
		cat("Unequal length data")
		return(-1)
	}
	y1=cbind(c9sol1[1:c9length1,1],c9sol1[c9sol1[,2]==1,3],c9sol1[c9sol1[,2]==2,3])
	y1=y1[species_list,]
	y2=cbind(c9sol2[1:c9length1,1],c9sol2[c9sol2[,2]==1,3],c9sol2[c9sol2[,2]==2,3])
	y2=y2[species_list,]
	npts=length(species_list)
	
	
	if(TRUE)
	{
		forder=order(y1[,2])
		y1=y1[forder,]
		y2=y2[forder,]
	}
	pt=dim(y1)
	labelactive=(class(tlabels)=="data.frame")
	if(doplot){
		if(overplot==FALSE)
		{
			hbar=0.4
			hoffset=0.45
	#	quartz()
		plot(c(0,npts+1),c(-(max(c9sol1[,3]))/5,(max(c9sol1[,3])+1)),type="n",xlab="Taxa Number",ylab="Position in Composite")
		  for(k in 1:(npts)){
		  	 lines(c(k,k+hbar),c(y1[k,2],y1[k,2]),col="black")
	         lines(c(k,k+hbar),c(y1[k,3],y1[k,3]),col="black")
		     lines(c(k,k),c(y1[k,2],y1[k,3]))	
		     lines(c(k,k)+hbar,c(y1[k,2],y1[k,3])) 
		     
		     lines(c(k,k+hbar)+hoffset,c(y2[k,2],y2[k,2]),col='red')
	         lines(c(k,k+hbar)+hoffset,c(y2[k,3],y2[k,3]),col="red")
		     lines(c(k,k)+hoffset,c(y2[k,2],y2[k,3]),col="red",lty=4)	
		     lines(c(k,k)+hbar+hoffset,c(y2[k,2],y2[k,3]),col='red',lty=4) 
		     
		     if(labelactive)
		     {
		     	if(addcounts)
		     	{
		     	  text(k+0.5,min(c(y1[k,2],y2[k,2])),paste(tlabels$TaxaName[y1[k,1]]," ",as.character(taxacounts[y1[k,1],2])," -"),adj=c(1,0),cex=0.6,srt=90,font=3) 
		     	}
		     	else
		     	{
		     	  text(k+0.5,min(c(y1[k,2],y2[k,2])),paste(tlabels$TaxaName[y1[k,1]],"-"),adj=c(1,0),cex=0.6,srt=90,font=3)	
		     	}
		     }
		  }
		}
		else
		{	
		  hoffset=0.45
		  
		  for(k in 1:(npts)){
		     lines(c(k,k+hbar),c(y1[k,2],y1[k,2]),col="black")
	         lines(c(k,k+hbar),c(y1[k,3],y1[k,3]),col="black")
		     lines(c(k,k),c(y1[k,2],y1[k,3]))	
		     lines(c(k,k)+hbar,c(y1[k,2],y1[k,3])) 
		     
		     lines(c(k,k+hbar),c(y2[k,2],y2[k,2]),col="black")
	         lines(c(k,k+hbar),c(y2[k,3],y2[k,3]),col="black")
		     lines(c(k,k),c(y2[k,2],y2[k,3]))	
		     lines(c(k,k)+hbar,c(y2[k,2],y2[k,3])) 
		  }
	
		}
		
	}
	if(labelactive)
	{
		y1=cbind(y1,tlabels$TaxaName[y1[,1]])
		colnames(y1)<-c("c9Code","FAD","LAD","SpeciesName")
	}
    return(y1)	
}
	
   


# get a summary of taxa range using the CH4 format, which has score, section, horizon number+horizon
# height in the first four columns
# xval and yval are points used to extrapolate the composite positions into time.
GetRangeSummary<-function(ch,xval,yval)
{
	pt=dim(ch)
	print("size of data set")
	print(pt)
	trange=apply(ch[,5:pt[2]],2,ColumnFADLAD)
	trange=matrix(unlist(trange),ncol=2,byrow=TRUE)
	sval1=ch[trange[,1],1]
	sval2=ch[trange[,2],1]
    cat(head(trange))
	y=1:(pt[2]-4)
	y=cbind(y,trange,sval1,sval2)
	if(length(xval)>2)
	{
		tval1=approxExtrap(xval,yval,trange[,1])
		tval2=approxExtrap(xval,yval,trange[,2])
		y=cbind(y,tval1$y,tval2$y)
	}
	return(y)
}

# plot a Range Chart with Taxa Labels

PlotRangeChartLabels<-function(j,TaxDic,doplot=TRUE,FullName=FALSE,FADOrder=TRUE){
	# j$pen, j$initpen,j$d, j$history
	# TaxDic is the taxanomic dictionary -$TNumber, $Abr,$Name
	pt=dim(j$d)
	stTaxa=5
	trange=apply(j$d[,stTaxa:pt[2]],2,ColumnFADLAD)
	trange=matrix(unlist(trange),ncol=2,byrow=TRUE)
    y=1:(pt[2]-stTaxa+1)
	y=cbind(y,trange)
	if(FADOrder){
	  tord=order(y[,2])
	  y=y[tord,]	
	}
	if(doplot){
		#quartz()
		if(FullName){
		    plot(c(0,(pt[2]-2)),c(-100,(pt[1]+1)),type="n",ylab="Ordinal Position",xlab="")
		  }else{
		  	  plot(c(0,(pt[2]-2)),c(-10,(pt[1]+1)),type="n",ylab="Ordinal Position",xlab="")
		  	}
		points(1:(pt[2]-4),y[,2],col="black")
		points(1:(pt[2]-4),y[,3],col="black")
		for(k in 1:(pt[2]-4)){
		   lines(c(k,k),c(y[k,2],y[k,3]))	
		   if(FullName){
		   	 text(k,y[k,2]-78,labels=TaxDic$Name[y[k,1]],cex=0.7,srt=90,font=3)
		   }else{
		     text(k,y[k,2]-9,labels=TaxDic$Abr[y[k,1]],cex=0.7,srt=90,font=3)
		   }
		}
	}
	
    return(y)	
}



#columnwise listing of FADLAD values

ColumnFADLAD<-function(y){
	lp=which(y==1)
	if(length(lp)==1){
		x=c(lp,lp)
	}
	else{
		x=c(lp[1],lp[length(lp)])
	}
	return(x)
}

# routine for showing history of HA solution
PlotHAHistory<-function(y){
	quartz()
	pmin=min(y$history[,2])
	pmax=max(y$history[,3])
	plot(c(1,length(y$history[,1])),c(pmin,pmax))
	points(y$history[,3],col='red')
	points(y$history[,2],col='green')
}

#compute a CONOP style output listing from the D3 form of the HA solution

CH2C9Sol<-function(j){
	# note j is the solution output from HA-
	# j$pen, j$initpen,j$d, j$history
	# call PlotRangeChart(j,doplot=FALSE) to get the columns, TAXA#, FAD position, LAD position
	y=PlotRangeChart(j,doplot=FALSE)
	ypt=dim(y)
	z=cbind(y[,1],rep(1,ypt[1]),y[,2])
	z2=cbind(y[,1],rep(2,ypt[1]),y[,3])
	znet=rbind(z,z2)
	#print(znet)
	zord=order(znet[,3])
	y=cbind(znet[zord,1:2],1:(2*ypt[1]))
	
}
#compute a CONOP style output listing from the D3 form of the HA solution
# this version uses PlotRangeChartPS2 to compute the FADs/LADs

CH2C9Sol2<-function(j_ch,species_list,tlabels=c("na")){	
	# call PlotRangeChartPS2(j,doplot=FALSE) to get the columns, TAXA#, FAD position, LAD position
	y=PlotRangeChartPS2(j_ch,species_list,doplot=FALSE,overplot=FALSE,hbar=0.66,tlabels,barOrder="None")
	ypt=dim(y)
	z=cbind(y[,1],rep(1,ypt[1]),y[,2])
	z2=cbind(y[,1],rep(2,ypt[1]),y[,3])
	znet=rbind(z,z2)
	#print(znet)
	zord=order(znet[,3])
	y=cbind(znet[zord,1:2],1:(2*ypt[1]))	
}



# plot section range charts

PlotSectionRange<-function(y,ordinal=TRUE,sectNames=c(NA,NA),overplot=FALSE){
	par(pch='_')
	# set up range
	
	if(overplot)
	{
	  hbar=0.4
	}
	else
	{
		hbar=0.66
	}  
	nsections=max(y$d[,2])
	npts=dim(y$d)[1]
	sectmax=rep(0,nsections)
	sectmin=rep(1e9,nsections)
	labelsActive=(class(sectNames)=="data.frame")
	lspace=0
    
	if(ordinal)
	{
		if(labelsActive)
		{
		  plot(c(0,max(y$d[,2])+1),c(-(dim(y$d)[1]/5),dim(y$d)[1]+2),type="n",xlab="Section Number",ylab="Position in Composite")
		}
		else
		{
		  plot(c(0,max(y$d[,2])+1),c(0,dim(y$d)[1]+2),type="n",xlab="Section Number",ylab="Position in Composite")	
		}  
		dpt=dim(y$d)
		dord=1:dpt[1]
		for(i in 1:npts)
		{
		  lines(c(y$d[i,2],y$d[i,2]+hbar),c(dord[i],dord[i]))
		  if(dord[i]<sectmin[y$d[i,2]]){
		     sectmin[y$d[i,2]]=dord[i]
		     }
		  if(dord[i]>sectmax[y$d[i,2]]){
		     sectmax[y$d[i,2]]=dord[i]
		     }
		}
		for(i in 1:nsections)
		{
			lines(c(i,i),c(sectmin[i],sectmax[i]))
			lines(c(i,i)+hbar,c(sectmin[i],sectmax[i]))
			if(labelsActive)
		    {
		      text(i+0.5,sectmin[i],paste(sectNames$SectionName[i]," -"),adj=c(1,0),cex=0.6,srt=90,font=3) 
		    }
			
		}
		
	}
	else{	
		 if(labelsActive)
		 {
		   plot(c(0,max(y$d[,2])+1),c(-0.2,max(y$d[,1]+0.1)),type="n",xlab="Section Number",ylab="Position in Composite")
		 }
		 else
		 {
		   	plot(c(0,max(y$d[,2])+1),c(0,max(y$d[,1]+0.1)),type="n",xlab="Section Number",ylab="Position in Composite")
		 }
		 for(i in 1:npts)
		 {
	      lines(c(y$d[i,2],y$d[i,2]+hbar),c(y$d[i,1],y$d[i,1]))
	      if(y$d[i,1]<sectmin[y$d[i,2]]){
		     sectmin[y$d[i,2]]=y$d[i,1]
		     }
		  if(y$d[i,1]>sectmax[y$d[i,2]]){
		     sectmax[y$d[i,2]]=y$d[i,1]
		     }
	    }
	    for(i in 1:nsections)
		{
			lines(c(i,i),c(sectmin[i],sectmax[i]))
			lines(c(i,i)+hbar,c(sectmin[i],sectmax[i]))
			if(labelsActive)
		    {
		      text(i+0.5,sectmin[i],paste(sectNames$SectionName[i]," -"),adj=c(1,0),cex=0.6,srt=90,font=3) 
		    }
		}
	  } 
	  
	  
}

# plot section range charts
# this version shows the number of stratigraphically informative
# taxa supporting the positioning of the horizon

PlotSectionRange2<-function(y,ordinal=TRUE,sectNames=c(NA,NA),overplot=FALSE){
	par(pch='_')
	# set up range
	
	if(overplot)
	{
	  hbar=0.25
	}
	else
	{
		hbar=0.5
	}  
	nsections=max(y$d[,2])
	npts=dim(y$d)[1]
	sectmax=rep(0,nsections)
	sectmin=rep(1e9,nsections)
	labelsActive=(class(sectNames)=="data.frame")
	lspace=0
    
	if(ordinal)
	{
		if(labelsActive)
		{
		  plot(c(0,max(y$d[,2])+1),c(-(dim(y$d)[1]/5),dim(y$d)[1]+2),type="n",xlab="Section Number",ylab="Position in Composite")
		}
		else
		{
		  plot(c(0,max(y$d[,2])+1),c(0,dim(y$d)[1]+2),type="n",xlab="Section Number",ylab="Position in Composite")	
		}  
		dpt=dim(y$d)
		dord=1:dpt[1]
		for(i in 1:npts)
		{
		  lines(c(y$d[i,2],y$d[i,2]+hbar),c(dord[i],dord[i]))
		  stext=sprintf("%i",as.integer(sum(y$d[i,-(1:4)]>0,na.rm=TRUE)))
		  text(y$d[i,2]+hbar*1.3, dord[i],stext,cex=1,col="red")
		  if(dord[i]<sectmin[y$d[i,2]]){
		     sectmin[y$d[i,2]]=dord[i]
		     }
		  if(dord[i]>sectmax[y$d[i,2]]){
		     sectmax[y$d[i,2]]=dord[i]
		     }
		}
		for(i in 1:nsections)
		{
			lines(c(i,i),c(sectmin[i],sectmax[i]))
			lines(c(i,i)+hbar,c(sectmin[i],sectmax[i]))
			if(labelsActive)
		    {
		      text(i+0.5,sectmin[i],paste(sectNames$SectionName[i]," -"),adj=c(1,0),cex=1,srt=90,font=3) 
		    }
			
		}
		
	}
	else{	
		 if(labelsActive)
		 {
		   plot(c(0,max(y$d[,2])+1),c(-0.2,max(y$d[,1]+0.1)),type="n",xlab="Section Number",ylab="Position in Composite")
		 }
		 else
		 {
		   	plot(c(0,max(y$d[,2])+1),c(0,max(y$d[,1]+0.1)),type="n",xlab="Section Number",ylab="Position in Composite")
		 }
		 for(i in 1:npts)
		 {
	      lines(c(y$d[i,2],y$d[i,2]+hbar),c(y$d[i,1],y$d[i,1]))
	      if(y$d[i,1]<sectmin[y$d[i,2]]){
		     sectmin[y$d[i,2]]=y$d[i,1]
		     }
		  if(y$d[i,1]>sectmax[y$d[i,2]]){
		     sectmax[y$d[i,2]]=y$d[i,1]
		     }
	    }
	    for(i in 1:nsections)
		{
			lines(c(i,i),c(sectmin[i],sectmax[i]))
			lines(c(i,i)+hbar,c(sectmin[i],sectmax[i]))
			if(labelsActive)
		    {
		      text(i+0.5,sectmin[i],paste(sectNames$SectionName[i]," -"),adj=c(1,0),cex=1,srt=90,font=3) 
		    }
		}
	  } 
	  
	  
}



#--------------------------------
# adding a second section for comparison
PlotSectionRange2nd<-function(y,ordinal=TRUE){
	par(pch='_')
	# set up range
	
	
	  hbar=0.4
	  hoffset=0.45
	
	 
	nsections=max(y$d[,2])
	npts=dim(y$d)[1]
	sectmax=rep(0,nsections)
	sectmin=rep(1e9,nsections)
	labelsActive=FALSE
	lspace=0
    
	if(ordinal)
	{
		 
		dpt=dim(y$d)
		dord=1:dpt[1]
		for(i in 1:npts)
		{
		  lines(c(y$d[i,2],y$d[i,2]+hbar)+hoffset,c(dord[i],dord[i]),col='blue')
		  if(dord[i]<sectmin[y$d[i,2]]){
		     sectmin[y$d[i,2]]=dord[i]
		     }
		  if(dord[i]>sectmax[y$d[i,2]]){
		     sectmax[y$d[i,2]]=dord[i]
		     }
		}
		for(i in 1:nsections)
		{
			lines(c(i,i)+hoffset,c(sectmin[i],sectmax[i]),col='blue',lty=4)
			lines(c(i,i)+hbar+hoffset,c(sectmin[i],sectmax[i]),col='blue',lty=4)
		
		}
		
	}
	else{	
		 
		 for(i in 1:npts)
		 {
	      lines(c(y$d[i,2],y$d[i,2]+hbar)+hoffset,c(y$d[i,1],y$d[i,1]),col='blue')
	      if(y$d[i,1]<sectmin[y$d[i,2]]){
		     sectmin[y$d[i,2]]=y$d[i,1]
		     }
		  if(y$d[i,1]>sectmax[y$d[i,2]]){
		     sectmax[y$d[i,2]]=y$d[i,1]
		     }
	    }
	    for(i in 1:nsections)
		{
			lines(c(i,i)+hoffset,c(sectmin[i],sectmax[i]),col='blue',lty=4)
			lines(c(i,i)+hbar+hoffset,c(sectmin[i],sectmax[i]),col='blue',lty=4)
			
		}
	  } 
}



# this is a tool to overplot several different range charts on one graph

PlotSectionRangeRed<-function(y,ordinal=TRUE){
	par(pch='_')
	par(col="red")
	if(ordinal)
	{
		dpt=dim(y$d)
		dord=1:dpt[1]
		 points(y$d[,2]+0.5,dord,xlab='Section Number',ylab='Position in Composite (Ordinal)')
	}
	else{
	  points(y$d[,2]+0.5,y$d[,1],xlab='Section Number',ylab='Position in Composite')
	  }
}






# plot range extension per species by section
# visual depiction of where the error is,  input y is the output structure from HA
# note: assumes score is column 1, section is column 2, 3 horizon number, 5 is horizon height, followed by all taxa


ShowRangeExtensions<-function(y_ch,TaxDic=0,UseTaxDic=FALSE,SecDic=0,UseSecDic=FALSE,DoPlot=TRUE,boxscale=20){
	# get dimensions of grand data matrix
	dimdat=dim(y_ch)
	# figure out number of sections and number of taxa from size of data matrix
	nsections=max(y_ch[,2])
	ntaxa=dimdat[2]-4
	rexmat=matrix(NA,nrow=nsections,ncol=ntaxa)
	duration=rep(0,ntaxa)
	for(i in 1:ntaxa){
	  lp=which(y_ch[,4+i]>0)
	  lpv=lp[1]:lp[length(lp)]
	  duration[i]=lp[length(lp)]-lp[1]+1
	  for(j in lpv){
	  	  if(!is.na(y_ch[j,4+i])){
	  	  	if(is.na(rexmat[y_ch[j,2],i]))
	  	  	{
	  	  		 rexmat[y_ch[j,2],i]=0
	  	  	}	
	    	if(y_ch[j,i+4]==0){
	  	    	rexmat[y_ch[j,2],i]=rexmat[y_ch[j,2],i]+1
	  	    }else{
	  		    rexmat[y_ch[j,2],i]=rexmat[y_ch[j,2],i]+0
	  	    }
	  	}
	  	
	  }
	}
	# rescale durations
	if(DoPlot==FALSE)
	{
		return(rexmat)
	}
	
	
	# now plot rexmat
	
	#boxscale=20
	plot(c(0,ntaxa+1),c(0,nsections+1),type='n',xlab='Taxa Number',ylab='Section Number')
	for(i in 1:ntaxa){
		for(j in 1:nsections){
			if(is.na(rexmat[j,i])){
				points(i,j,pch='-')
			}
			else{
				expa=rexmat[j,i]/boxscale +0.25
				points(i,j,pch=0,cex=expa,col='red')
			}
		}
		lines(c(i,i),c(nsections+0.05,(duration[i]/max(duration))+nsections+0.05))
	}
	if(UseTaxDic){
		text(1:ntaxa,0.2,labels=TaxDic$Abr,cex=0.7,srt=90,font=3)
	}
    if(UseSecDic){
    	text(0.1,1:nsections,labels=SecDic$Abr,cex=0.7)
    }
    
	return(rexmat)
}

# plot range extension per species by section
# visual depiction of where the error is,  input y is the output structure from HA
# note: assumes score is column 1, section is column 2, 3 horizon number, 5 is horizon height, followed by all taxa
# this version is meant to use the full name of the taxa and sections-uses more screen

ShowRangeExtensions2<-function(y_ch,TaxDic=0,UseTaxDic=FALSE,SecDic=0,UseSecDic=FALSE,DoPlot=TRUE,boxscale=20){
	# get dimensions of grand data matrix
	dimdat=dim(y_ch)
	# figure out number of sections and number of taxa from size of data matrix
	nsections=max(y_ch[,2])
	ntaxa=dimdat[2]-4
	rexmat=matrix(NA,nrow=nsections,ncol=ntaxa)
	duration=rep(0,ntaxa)
	for(i in 1:ntaxa){
	  lp=which(y_ch[,4+i]>0)
	  lpv=lp[1]:lp[length(lp)]
	  duration[i]=lp[length(lp)]-lp[1]+1
	  for(j in lpv){
	  	  if(!is.na(y_ch[j,4+i])){
	  	  	if(is.na(rexmat[y_ch[j,2],i]))
	  	  	{
	  	  		 rexmat[y_ch[j,2],i]=0
	  	  	}	
	    	if(y_ch[j,i+4]==0){
	  	    	rexmat[y_ch[j,2],i]=rexmat[y_ch[j,2],i]+1
	  	    }else{
	  		    rexmat[y_ch[j,2],i]=rexmat[y_ch[j,2],i]+0
	  	    }
	  	}
	  	
	  }
	}
	# rescale durations
	if(DoPlot==FALSE)
	{
		return(rexmat)
	}
	
	
	# now plot rexmat
	
	#boxscale=20
	plot(c(-8,ntaxa+1),c(-6,nsections+1),type='n',xlab='Taxa Number',ylab='Section Number')
	for(i in 1:ntaxa){
		for(j in 1:nsections){
			if(is.na(rexmat[j,i])){
				points(i,j,pch='-')
			}
			else{
				expa=rexmat[j,i]/boxscale +0.25
				points(i,j,pch=0,cex=expa,col='red')
			}
		}
		lines(c(i,i),c(nsections+0.05,(duration[i]/max(duration))+nsections+0.05))
	}
	if(UseTaxDic){
		text(1:ntaxa,-5.8,labels=TaxDic$Abr,cex=0.8,srt=90,font=3)
	}
    if(UseSecDic){
    	text(-6.5,1:nsections,labels=SecDic$Abr,cex=0.8)
    }
    
	return(rexmat)
}


# now a variant on ShowRangeExtension to show range extensions over time...

ShowRangeExtensionsTime<-function(yin,TaxDic=0,UseTaxDic=FALSE,SecDic=0,UseSecDic=FALSE,timebins=20,DoPlot=TRUE,boxscale=20){
	# note use of timebins to indicate the number of temporal bins to use
	# get dimensions of grand data matrix
	y=list(d=yin)
	dimdat=dim(y$d)
	# figure out number of sections and number of taxa from size of data matrix
	nsections=max(y$d[,2])
	ntaxa=dimdat[2]-3
	nhorizons=dimdat[1]
	rexmat=matrix(NA,nrow=timebins,ncol=ntaxa)
   
	duration=rep(0,ntaxa)
	for(i in 1:ntaxa){
	  lp=which(y$d[,3+i]>0)
	  lpv=lp[1]:lp[length(lp)]
	  duration[i]=lp[length(lp)]-lp[1]+1
	  for(j in lpv){
	  	  # figure out current bin
	  	  cbin=ceiling((j/nhorizons)*timebins)
	  	  if(!is.na(y$d[j,3+i])){
	  	  	if(is.na(rexmat[cbin,i]))
	  	  	{
	  	  		 rexmat[cbin,i]=0
	  	  	}	
	    	if(y$d[j,i+3]==0){
	  	    	rexmat[cbin,i]=rexmat[cbin,i]+1
	  	    }else{
	  		    rexmat[cbin,i]=rexmat[cbin,i]+0
	  	    }
	  	}
	  	
	  }
	}
	# just return rexmat if plot not desired
    if(DoPlot==FALSE)
    {
    	return(rexmat)
    }

	# rescale durations
	
	# now plot rexmat
	
	#boxscale=20- now a defaulted input parameter
	plot(c(0,ntaxa+1),c(0,timebins+1),type='n',xlab='Taxa Number',ylab='Time Bin')
	for(i in 1:ntaxa){
		for(j in 1:timebins){
			if(is.na(rexmat[j,i])){
				points(i,j,pch='-')
			}
			else{
				expa=rexmat[j,i]/boxscale +0.25
				points(i,j,pch=0,cex=expa,col='red')
			}
		}
		lines(c(i,i),c(timebins+0.05,(duration[i]/max(duration))+timebins+0.05))
	}
	if(UseTaxDic){
		text(1:ntaxa,0.2,labels=TaxDic$Abr,cex=0.7,srt=90,font=3)
	}
    if(UseSecDic){
    	text(0.1,1:nsections,labels=SecDic$Abr,cex=0.7)
    }
    
	return(rexmat)
}

# comparison plot of range extensions, shows the differences in range extensions between solutions
# by taxa, and either time or section number- note you must use ShowRangeExtensions, or ShowRangeExtentionsTime
# to generate the rexmat matrices are then compared using ShowDifRangeExtension

ShowDifRangeExtensions<-function(y1,y2,TaxDic=0,UseTaxDic=FALSE,SecDic=0,UseSecDic=FALSE, boxscale=10){
	# note use of timebins to indicate the number of temporal bins to use
	# get dimensions of grand data matrix
	
	#find dimensions
	y1dim=dim(y1)
	ntaxa=y1dim[2]
	timebins=y1dim[1]
	
	rexmat=y2-y1		#change is from 1 to 2,  or final minus initial
	
	# now plot rexmat
	
	#boxscale=10
	plot(c(0,ntaxa+1),c(0,timebins+1),type='n',xlab='Taxa Number',ylab='Time Bin')
	for(i in 1:ntaxa){
		for(j in 1:timebins){
			if(is.na(rexmat[j,i])){
				points(i,j,pch='-')
			}
			else{
				if(rexmat[j,i]>0)
				{
				  expa=rexmat[j,i]/boxscale +0.25
				  points(i,j,pch=0,cex=expa,col='red')
				}
				else
				{
				  expa=abs(rexmat[j,i]/boxscale) +0.25
				  points(i,j,pch=1,cex=expa,col='green')

				}
			}
		}
		#lines(c(i,i),c(timebins+0.05,(duration[i]/max(duration))+timebins+0.05))
	}
	if(UseTaxDic){
		text(1:ntaxa,0.2,labels=TaxDic$Abr,cex=0.7,srt=90,font=3)
	}
    if(UseSecDic){
    	text(0.1,1:nsections,labels=SecDic$Abr,cex=0.7)
    }
    
	return(rexmat)
}





# routine to show the slotting of isotope data
# input is a j structure output out of of HA-
# column 1 is the relative position in the composite
# column 2 is the variable to be slotted
# column 3 is the section number

ShowIsotopeSlotting<-function(j){
	nsections<-max(j$d[,3])
	par(mfcol=c(1,nsections+1))
	cpts=1:length(j$d[,1])
	plot(j$d[,2],cpts,xlab="",ylab="",type="n",main="Composite") 
	lines(j$d[,2],cpts)
	points(j$d[,2],cpts)
	for(i in 1:nsections){
		plot(j$d[,2],cpts,xlab="",ylab="",type="n")
	    lines(j$d[j$d[,3]==i,2],cpts[j$d[,3]==i])
	    points(j$d[j$d[,3]==i,2],cpts[j$d[,3]==i])
	}
}

# compare stratigraphic solutions in C9 format
#

CompareC9Solutions<-function(sol1,sol2){
# sol1, sol2 are 3 column matrices.  Col 1 is the event number, Col2 is the type and Col3 is the 
# ordinal position in the solution
  ord1=order(sol1[,1]+10000*sol1[,2])
  sol1=sol1[ord1,]
  ord2=order(sol2[,1]+10000*sol2[,2])
  sol2=sol2[ord2,]
  cval=cor.test(sol1[,3],sol2[,3])
  cval2=cor.test(sol1[,3],sol2[,3],method="kendall")
  print(cval2)
  cval3=kStress(sol1[,3],sol2[,3])
  return(list(c1=cval,c2=cval2,stress=cval3))

	
}

# k-stress calculation
kStress<-function(x,y)
{
	n=length(x)
	dx=matrix(0,nrow=n,ncol=n)
	dy=dx
	for(i in 1:(n-1))
	{
		for(j in (i+1):n)
		{
			dx[i,j]=abs(x[i]-x[j])
			dx[j,i]=dx[i,j]
			dy[i,j]=abs(y[i]-y[j])
			dy[j,i]=dy[i,j]
		}
	}
	z=100*sum( (dx-dy)^2)/(sum(dx^2)*sum(dy^2))^0.5
}


# plot range charts of two sections

DualPlotRangeChart<-function(j,j2,doplot=TRUE){
	# j$pen, j$initpen,j$d, j$history
	pt=dim(j$d)
	print("size of data set")
	print(pt)
	trange=apply(j$d[,3:pt[2]],2,ColumnFADLAD)
	trange=matrix(unlist(trange),ncol=2,byrow=TRUE)
    y=1:(pt[2]-2)
	y=cbind(y,trange)
	trange2=apply(j2$d[,3:pt[2]],2,ColumnFADLAD)
	trange2=matrix(unlist(trange2),ncol=2,byrow=TRUE)
    y2=1:(pt[2]-2)
	y2=cbind(y2,trange2)
	print(y)
	print(y2)
	if(doplot){
		#quartz()
		par(cex=2,cex.axis=0.5,cex.lab=0.5)
		plot(c(0,(pt[2]-2)),c(0,(pt[1]+1)),type="n",xlab="Taxa Number",ylab="Position in Composite")
		points(y[,1],y[,2],col="black")
		points(y[,1],y[,3],col="black")
		points(y2[,1]+0.25,y2[,2],col="red")
		points(y2[,1]+0.25,y2[,3],col="red")
		for(k in 1:(pt[2]-2)){
		   lines(c(y[k,1],y[k,1]),c(y[k,2],y[k,3]),col="black")	
		   lines(c(y2[k,1]+0.25,y2[k,1]+0.25),c(y2[k,2],y2[k,3]),col="red")	
		}
	}
	
    return(list(y=y,y2=y2,t=trange,t2=trange2))	
}


# plot section range charts

DualPlotSectionRange<-function(y,y2,ordinal=TRUE){
	# 
	par(pch='_',cex=2,cex.axis=0.5)
	dpt=dim(y$d)
	print(max(y$d[,2]))
	plot(c(1,max(y$d[,2])+1),c(0,dpt[1]+2),type="n",xlab='Section Number',ylab='Position in Composite (Ordinal)')
	if(ordinal)
	{
		
		dord=1:dpt[1]
		 
		 points(y$d[,2],dord)
		 points(y2$d[,2]+0.25,dord,col="red")
	}
	else{
	  points(y$d[,2],y$d[,1])
	  points(y2$d[,2]+0.25,y2$d[,1],col="red")
	  
	  }
}


# alteration to handle data columns with no 1 values in them

ColumnRangeExtension2<-function(y){
	# compute the range extension implied by a column of data, sorted in time
	lp=which(y>0)
	if(length(lp)==0)
	{return(0)}
	lpv=lp[1]:lp[length(lp)]
	sum(y[lpv]==0,na.rm=TRUE)
}
# Net Range Extension calculation using ColumnRangeExtension2- allowing for empty columns
NetRangeExtension2<-function(y){
	
	sum(apply(y,2,ColumnRangeExtension2))
	#rx=parallel(sum(apply(y,2,ColumnRangeExtension2)))
	#collect(rx)
	#return(rx)
}

# determining which taxa are stratigraphically informative
# given a capture history, with a penalty structure, figure out how many sections each taxa
# in the penalty structure are in
InformativeTaxa<-function(ch,taxaList)
{
	# ch is the capture history, column 1 is the score, column 2 is the section number
	# taxaList is the list of taxa that should be used
	
	taxaNumbers=1:length(taxaList)
	sectcnt=rep(0,length(taxaList))
	for(i in taxaNumbers)
	{
		sectv=unique(ch[ (ch[,taxaList[i]+4]==1)&(!is.na(ch[,taxaList[i]+4]==1)),2])
		sectcnt[i]=length(sectv)
	}
	y=cbind(taxaNumbers,sectcnt)
	return(y)
}



SummarizeErrorHA<-function(d3,pen_str)
{
	#  Summarize the HA error for the data set d3 and the penalty function pen_str
	# penalty_spec=list(n_biostrat=n_biostrat,biostrat=biostrat,n_pmag=n_pmag,pmag=pmag,n_dates=n_dates,dates=dates,n_ashes=n_ashes,n_continuous=n_continuous,continuous=continuous)
	
	
	
  # assumes that the data has position scores in the first column, reorder the data to match	
  

  d3ord=order(d3[,1])
  d3a=d3[d3ord,]
  d3a[,1]=(d3a[,1]-min(d3a[,1]))/(max(d3a[,1])-min(d3a[,1]))
  
  
pcv3=0 
if(pen_str$n_biostrat)
{  print(sprintf("Net taxa range extension %f",NetRangeExtension(d3a[,2+pen_str$biostrat])))
	pcv3=NetRangeExtension(d3a[,2+pen_str$biostrat]) }
	
pcvm=0
if(pen_str$n_pmag>0)
 {
   for(i in pen_str$pmag)
   {
      pcvm=+ftransitions2(d3a[,i+2])*30
   }    
 }	
 pcv3=pcv3+pcvm
 print(sprintf("Pmag Penalty: %f",pcvm))
 print(sprintf("Pmag Reversals: %f",pcvm/30))
 acv=0
if(pen_str$n_ashes>0)
  {
  		for(i in pen_str$ashes[,1])
  		{
  		    acv=acv+AshRangeExtensionStrong(d3a[,i+2])*1000
  		    # note AshRangeExtension and AshRangeExtensionStrong are slightly different
  		}    
  }
  pcv3=pcv3+acv
print(sprintf("Ash term: %i",acv))
print(sprintf("Ash levels %i",acv/1000))

        ccv=0;
		if(pen_str$n_continuous>0)
  		{
  			for(i in pen_str$continuous[,1])
  			{
  		    	ccv=ccv+ftransitions2(d3a[,i+2])
  			}    
  		}
print(sprintf("Continuous variable term: %f ",ccv))
pcv3=pcv3+ccv
		pcv=0
		if(pen_str$n_dates>0)
  		{
  			for(i in 1:pen_str$n_dates)
  			{
  		    	pcv=pcv+PassingError(d3a[,pen_str$dates[i,1]+2],pen_str$dates[i,2],d3a[,pen_str$dates[i,3]+2],pen_str$dates[i,4])*pen_str$dates[i,5]
  			}    
  		}
print(sprintf("Passing Error: %i",pcv))
pcv3=pcv3+pcv
 print(sprintf("Best Total Penalty: %f ", pcv3))
 }