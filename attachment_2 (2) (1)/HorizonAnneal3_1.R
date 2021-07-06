# Version 3 of Horizon Annealing in R
# allows for a variable number of geological variables- chemostratigraphy, sequence stratigraphy, paleomagnetic

# modified 5/14/2018
# call this version 3_1 meaning 3.1


HorizonAnneal3_1<-function(d3,param=list(nouter=400,ninner=100,temperature=5,cooling=0.9),pen_str)
{
	#  simulated annnealing based on biostratigraphy+chemostratigraphy
	# d3[,1] - score
	# d3[,2] -section codes
  # d3[,3]- horizon number
  # d3[,4]- horizon height
	# d3[,5+1:lastcolumn] - geological variables- of all types
  
  
  # defines location of data relative to the first column of the matrix
  data_offset=4
	
	# param is the list of parameter values controlling the annealing
	#
	# pen_str is the penalty structure specifying how the penalty is to be calculated
	# penalty_spec=list(n_biostrat=n_biostrat,biostrat=biostrat,n_pmag=n_pmag,pmag=pmag,n_dates=n_dates,dates=dates,n_ashes=n_ashes,n_continuous=n_continuous,continuous=continuous)
	
	
	
   ptm <- proc.time()

  # assumes that the data has position scores in the first column, reorder the data to match	
  

  d3ord=order(d3[,1])
  d3a=d3[d3ord,]
  d3a[,1]=(d3a[,1]-min(d3a[,1]))/(max(d3a[,1])-min(d3a[,1]))
  gaprange=mean(diff(d3a[,1]))						# size of typical gap	
 
  movetrack=c(0,0,0,0,0)
  movetry=c(0,0,0,0,0)
  cat("Initial Penalty Calculation","\n")
 					# compute the biostratigraphic range extension  
  print(pen_str$biostrat+data_offset)
  if(pen_str$n_biostrat>0)
  {
    cv=NetRangeExtension(d3a[,(data_offset+pen_str$biostrat)])
  }
  else
  {
  	cv=0
  }
  cat('biostrat error was ',cv,'\n')
  				#compute the paleomagnetic or binary range extension
  if(pen_str$n_pmag==1)
  {
  		    cv=cv+ftransitions2(d3a[,pen_str$pmag+data_offset])*30    
  }	
  else
  {
  	  if(pen_str$n_pmag>1)
  	  {
  	  	cv=cv+sum(apply(d3a[,pen_str$pmag+data_offset],2,ftransitions2))*30
  	  }
  }
  				# compute the AshRangeExtensionPenalty
  if(pen_str$n_ashes==1)
  {
  		    cv=cv+AshRangeExtension(d3a[,pen_str$ashes[,1]+data_offset])*10  
  		    cat('value of penalty string on ashes',pen_str$ashes[,1]+2,"\n")
  		    cat("ash was ",AshRangeExtension(d3a[,pen_str$ashes[,1]+2]), "\n" )  
  }
  else
  {
  	   if(pen_str$n_ashes>1)
  	   {
  	   	  cv=cv+sum(apply(d3a[,pen_str$ashes[,1]+data_offset],2,AshRangeExtension))*10  
  	   }
  }
  				# compute the continuous variable penalty
  
  if(pen_str$n_continuous==1)
  {
  		    cv=cv+ftransitions2(d3a[,pen_str$continuous[,1]+data_offset])*pen_str$continuous[,2]   
  }
  else
  {
  	if(pen_str$n_continuous>1)
  	{
  		    cat("in calculation of continuous, n>1")
  			cv=cv+sum(apply(d3a[,pen_str$continuous[,1]+data_offset],2,ftransitions2))*pen_str$continuous[1,2]
  			cat(cv)
  	}
  }
  
  			#compute the "Passing Penalty"
  
  if(pen_str$n_dates>0)
  {
  		for(i in 1:pen_str$n_dates)
  		{
  		    cv=cv+PassingError2(d3a[,pen_str$dates[i,1]+data_offset],pen_str$dates[i,2],d3a[,pen_str$dates[i,3]+data_offset],pen_str$dates[i,4])*pen_str$dates[i,5]
  		    cat(i,pen_str$dates[i,1]+2,pen_str$dates[i,2],pen_str$dates[i,3]+2,pen_str$dates[i,4])
  		    cat(i,'Passing Error', PassingError2(d3a[,pen_str$dates[i,1]+2],pen_str$dates[i,2],d3a[,pen_str$dates[i,3]+2],pen_str$dates[i,4])*pen_str$dates[i,5],"\n")
  		}    
  }	   
  
  bestcv=cv
  cat("Starting penalty")
  print(bestcv)
  bestd3=d3a
  temperature=param$temperature
  nsections=max(d3a[,2])
  history=matrix(0,ncol=3,nrow=param$nouter)
  nhorizons=dim(d3a)[1]

  
for(j in 1:param$nouter)
{
  for(i in 1:param$ninner)
  {
  	 	pd3=d3a
		# propose change
		psec=floor(runif(1,min=1.000001,max=nsections+0.99999))      # pick the section
		pmove=runif(1)
		#cat(j,i,"\n")
		#cat(psec,pmove,"\n")
		if(pmove<0.2)             #shift up or down
		{							
			shval=runif(1,min=-0.1,max=0.1)
			ps=(pd3[,2]==psec)
			pd3[ps,1]=pd3[ps,1]+shval
			movetype=1
			movetry[1]=movetry[1]+1
			
		}											#end if
		else
		{
		if(pmove<0.4)                      #expand/contract
		{
			shval=runif(1,min=-0.05,max=0.05)+1;
			ps=(pd3[,2]==psec)
			pd3min=mean(pd3[ps,1])
			pd3[ps,1]=(pd3[ps,1]-pd3min)*shval+pd3min;
			movetype=2
			movetry[2]=movetry[2]+1
			
		}												#end if
		else									#insert or remove gap
		{
		if(pmove<0.6)	
		  {
			ps=(pd3[,2]==psec)
			while(sum(ps)<3)
			{
				psec=floor(runif(1,min=1.000001,max=nsections+0.99999))
				ps=(pd3[,2]==psec)
			}
		#	cat("Section ",psec,"\n")
			breakpoint=floor(runif(1,min=2,max=sum(ps)-0.001))
			w=which(ps)	
			gap=runif(1,min=-0.59,max=5)*(pd3[w[breakpoint+1],1]-pd3[w[breakpoint],1])
			w=w[-(1:breakpoint)]
			pd3[w,1]=pd3[w,1]+gap
			movetype=3
			movetry[3]=movetry[3]+1
		
		  }
		  else 
		  {	 
		  	if(pmove<0.8)								#insert dogleg  (was at 0.6)
		  	{
		  		shval=runif(1,min=-0.1,max=0.1)+1;
				ps=(pd3[,2]==psec)
			    while(sum(ps)<3)
			    {
				   psec=floor(runif(1,min=1.000001,max=nsections+0.99999))
				   ps=(pd3[,2]==psec)
			    }
				w=which(ps)
				breakpt=floor(runif(1,min=2,max=sum(ps)-0.001))
				gapval=diff(pd3[w,1])
				upchoice=runif(1,min=0,max=1)
				if(upchoice>0.5)
				{
					gapval[(breakpt+1):sum(ps)-1]=gapval[(breakpt+1):sum(ps)-1]*shval;
				}
				end
				{
					gapval[1:(breakpt)]=gapval[1:(breakpt)]*shval;
				}
				newval=cumsum(c(pd3[w[1],1],gapval))
				pd3[w,1]=newval
				movetype=4
				movetry[4]=movetry[4]+1
		  	} 				# end else at dogleg
		  	else				# insert shuffling routine
		  	{
		  		target=floor(runif(1,min=1.000001,max=nhorizons+.99))
		  		dmove=floor(runif(1,min=0.01,max=1.99))
		  		nmove=ceiling(abs(rnorm(1,0,4)))
		  		nmove=1
		  		movetype=5
		  		movetry[5]=movetry[5]+1

		  		if(dmove==0)					# shuffle up
		  		{
		  			startsection=pd3[target,2]
		  			while(nmove>0)
		  			{
		  				if(target==nhorizons)
		  				{
		  					nmove=0
		  				}
		  				else
		  				{
		  				#	cat("target ",target, pd3[target,2], pd3[target,1],"\n")
		  				#	cat("target+1",target+1,pd3[target+1,2],pd3[target+1,1],"\n")
		  					if(pd3[target+1,2]==startsection)
		  					{
		  						nmove=0
		  					}
		  					else
		  					{
		  						temp=pd3[target+1,1]
		  						pd3[target+1,1]=pd3[target,1]
		  						pd3[target,1]=temp
		  						target=target+1
		  						nmove=nmove-1
		  					}
		  				}     #end else on target
		  			}			#end while on nmove
		  			#absess=sum(!(pd3$DelC13[pd3$Section==1]==bestd3$DelC13[bestd3$Section==1]))
		  			#cat("after ripple",absess,"\n")
		  			#absess=sum(!(pd3$DelC13[pd3$Section==2]==bestd3$DelC13[bestd3$Section==2]))
		  			#cat("after ripple",absess,"\n")
		  			#absess=sum(!(pd3$DelC13[pd3$Section==3]==bestd3$DelC13[bestd3$Section==3]))
		  			#cat("after ripple",absess,"\n")
		  			#absess=sum(!(pd3$DelC13[pd3$Section==4]==bestd3$DelC13[bestd3$Section==4]))
		  			#cat("after ripple",absess,"\n")
		  			#absess=sum(!(pd3$DelC13[pd3$Section==5]==bestd3$DelC13[bestd3$Section==5]))
		  			#cat("after ripple",absess,"\n")
		  			
		  			
		  		}			# end of if on dmove
		  		else								#shuffle down
		  		{
		  			while(nmove>0)
		  			{
		  				if(target==1)
		  				{
		  					nmove=0
		  				}
		  				else
		  				{
		  					if(pd3[target-1,2]==pd3[target,2])
		  					{
		  						target=target-1
		  					}
		  					else
		  					{
		  						temp=pd3[(target-1),1]
		  						pd3[(target-1),1]=pd3[target,1]
		  						pd3[(target),1]=temp
		  						target=target-1
		  						nmove=nmove-1
		  					}
		  				}     #end else on target
		  			}			#end while on nmove
		  		}   # end of else based on dmove
		  	}
		  }
		  
		}												# end else at gap
		}												# end else at expand/contract
										# order the proposed solution
										
										
		#if(movetype==5)
	    #{
		#absess=sum(!(pd3$DelC13[pd3$Section==1]==bestd3$DelC13[bestd3$Section==1]))
		#cat("prior sort",absess,"\n")
		 # 			absess=sum(!(pd3$DelC13[pd3$Section==2]==bestd3$DelC13[bestd3$Section==2]))
		  #			cat("prior sort",absess,"\n")
		 # 			absess=sum(!(pd3$DelC13[pd3$Section==3]==bestd3$DelC13[bestd3$Section==3]))
		 ## 			cat("prior sort",absess,"\n")
		 # 			absess=sum(!(pd3$DelC13[pd3$Section==4]==bestd3$DelC13[bestd3$Section==4]))
		 # 			cat("prior sort",absess,"\n")
		 # 			absess=sum(!(pd3$DelC13[pd3$Section==5]==bestd3$DelC13[bestd3$Section==5]))
		 # 			cat("prior sort",absess,"\n")
		 # 			}								
		pd3ord=order(pd3[,1])
		pd3=pd3[pd3ord,]
	    #if(movetype==5)
	    #{
		#absess=sum(!(pd3$DelC13[pd3$Section==1]==bestd3$DelC13[bestd3$Section==1]))
		#cat("after sort",absess,"\n")
		#  			absess=sum(!(pd3$DelC13[pd3$Section==2]==bestd3$DelC13[bestd3$Section==2]))
		#  			cat("after sort",absess,"\n")
		# 			absess=sum(!(pd3$DelC13[pd3$Section==3]==bestd3$DelC13[bestd3$Section==3]))
		#  			cat("after sort",absess,"\n")
		#  			absess=sum(!(pd3$DelC13[pd3$Section==4]==bestd3$DelC13[bestd3$Section==4]))
		#  			cat("after sort",absess,"\n")
		#  			absess=sum(!(pd3$DelC13[pd3$Section==5]==bestd3$DelC13[bestd3$Section==5]))
		#  			cat("after sort",absess,"\n")
		#  			
		#  			}
#		abstotal=sum(!(pd3$DelC13[pd3$Section==1]==bestd3$DelC13[bestd3$Section==1]))+sum(!(pd3$DelC13#[pd3$Section==2]==bestd3$DelC13[bestd3$Section==2]))+sum(!(pd3$DelC13[pd3$Section==3]==bestd3$DelC13[bestd3$Section==3]))+sum(!(pd3$DelC13[pd3$Section==4]==bestd3$DelC13[bestd3$Section==4]))+sum(!(pd3$DelC13[pd3$Section==5]==bestd3$DelC13[bestd3$Section==5])) 			
#		if(abstotal>0)
#		{
#			cat(pd3ord)
#			cat(unique(pd3ord))
#			return(-1)
#		}
		
		
										# find the error for the proposed solution, starting with biostrat data
		if(pen_str$n_biostrat>0)
  		{
    		pcv3=NetRangeExtension(pd3[,data_offset+pen_str$biostrat])  		}
  		else
  		{
  			pcv3=0
  		}
 
         				#compute the paleomagnetic or binary range extension
  		if(pen_str$n_pmag==1)
  		{
  		    pcv3=pcv3+ftransitions2(pd3[,pen_str$pmag+data_offset])*30    
  		}	
  		else
  		{
  	  		if(pen_str$n_pmag>1)
  	  		{
  	  			pcv3=pcv3+sum(apply(pd3[,pen_str$pmag+data_offset],2,ftransitions2))*30
  	  		}
  		}

					# compute the AshRangeExtensionPenalty
  		if(pen_str$n_ashes==1)
  		{
  		    pcv3=pcv3+AshRangeExtension(pd3[,pen_str$ashes[,1]+data_offset])*10    
  		}
  		else
  		{
  	   		if(pen_str$n_ashes>1)
  	   		{
  	   	  		pcv3=pcv3+sum(apply(pd3[,pen_str$ashes[,1]+data_offset],2,AshRangeExtension))*10  
  	   		}
  		}
     				# compute the continuous variable range extension    
		if(pen_str$n_continuous>0)
  		{
  			for(i in pen_str$continuous[,1])
  			{
  		    	pcv3=pcv3+ftransitions2(pd3[,i+data_offset])*pen_str$continuous[1,2]
  			}    
  		}
  					# compute the passing error penalty
		if(pen_str$n_dates>0)
  		{
  			for(i in 1:pen_str$n_dates)
  			{
  		    	pcv3=pcv3+PassingError2(pd3[,pen_str$dates[i,1]+data_offset],pen_str$dates[i,2],pd3[,pen_str$dates[i,3]+data_offset],pen_str$dates[i,4])*pen_str$dates[i,5]
  			}    
  		}

									# now decide whether to accept the proposed solution	
		if(pcv3<=bestcv)						# if the proposed is better than the best ever seen, accept it
		{
			if(pcv3<bestcv)
			  { print(bestcv)
			  	movetrack[movetype]=movetrack[movetype]+1 }
			else
			  { cat("swop","\t")}
			bestcv=pcv3
			bestd3=pd3
			d3a=pd3
			cv=pcv3
			
		}else
		{
			pch=runif(1)
			if(pch<exp(-(pcv3-cv)/temperature))				# use the Bolzman factor to move up sometimes
			{
			 d3a=pd3
			 cv=pcv3
			}
			else
			{
								# don't accept the change
			}
		}	
							# force rescale
		#d3a[,1]=(d3a[,1]-min(d3a[,1]))/(max(d3a[,1])-min(d3a[,1]))
		d3a[,1]=(1:length(d3a[,1]))/length(d3a[,1])
  }
temperature=temperature*param$cooling
history[j,1]=temperature
history[j,2]=bestcv
history[j,3]=pcv3
cat("\nN outer: ",j, "T: ",temperature,"Best pen: ",bestcv,pmove,psec,"Recent prop pen: ", pcv3,"\n")
}
ptime=proc.time() - ptm
y=list(pen=bestcv,initpen=cv, d=bestd3,history=history)

print(ptime)
print(sprintf("Best Total Penalty: %f ", bestcv))
if(pen_str$n_biostrat)
{  print(sprintf("Net taxa range extension %f",NetRangeExtension(bestd3[,data_offset+pen_str$biostrat]))) }
pcv3=0
if(pen_str$n_pmag>0)
 {
   for(i in pen_str$pmag)
   {
      pcv3=pcv3+ftransitions2(bestd3[,i+data_offset])*30
   }    
 }	
 print(sprintf("Pmag Penalty: %f",pcv3))
 print(sprintf("Pmag Reversals: %f",pcv3/30))
 acv=0
if(pen_str$n_ashes>0)
  {
  		for(i in pen_str$ashes[,1])
  		{
  		    acv=acv+AshRangeExtension(bestd3[,i+data_offset])*10
  		}    
  }
print(sprintf("Ash term: %i",acv))
print(sprintf("Ash levels %i",acv/10))

        ccv=0;
		if(pen_str$n_continuous>0)
  		{
  			for(i in pen_str$continuous[,1])
  			{
  		    	ccv=ccv+ftransitions2(bestd3[,i+data_offset])
  			}    
  		}
print(sprintf("Continuous variable term: %f ",ccv))

		pcv=0
		if(pen_str$n_dates>0)
  		{
  			for(i in 1:pen_str$n_dates)
  			{
  		    	pcv=pcv+PassingError2(bestd3[,pen_str$dates[i,1]+data_offset],pen_str$dates[i,2],bestd3[,pen_str$dates[i,3]+data_offset],pen_str$dates[i,4])*pen_str$dates[i,5]
  			}    
  		}
print(sprintf("Passing Error: %i",pcv))
movetrack=1000*movetrack/movetry
print(sprintf("Move track values %f %f %f %f %f",movetrack[1],movetrack[2],movetrack[3],movetrack[4],movetrack[5]))
cat(gaprange)
#PlotHAHistory(y)
return(y)
}
