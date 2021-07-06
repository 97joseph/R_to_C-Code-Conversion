n_biostrat=62
biostrat=1:62
# for biostrat data, or Taxa FADs,LADs, the biostrat variable is the numbers of the columns with taxa

n_pmag=1
pmag=63

# pmag is a list of the column(s) with paleomagnetic signals, or really any binary data, NA values are not counted

n_dates=3
dates=matrix(c(64,0,65,0,1000,65,0,66,0,1000,66,0,67,0,1000),nrow=3,byrow=TRUE)

# each row of the dates matrix is a set of data to be entered into the passing penalty
# the first entry on each row is the column of the lower variale
# second entry on a row is the data type 0- singular date,  1- FAD, 2-LAD
# third and fourth entries on each row are the column and type of the second variable
# fifth value on each row is the weight

n_ashes=2
ashes=matrix(c(68,100,69,100),nrow=2,byrow=TRUE)
n_continuous=2
continuous=matrix(c(70,5,71,5),nrow=2,byrow=TRUE)
penalty_spec=list(n_biostrat=n_biostrat,biostrat=biostrat,n_pmag=n_pmag,pmag=pmag,n_dates=n_dates,dates=dates,n_ashes=n_ashes,ashes=ashes,n_continuous=n_continuous,continuous=continuous)

# penalty_spec is the list format for the penalty structure in the form of HA meant to handle multiple penalties
# the above example is set up to handle the Riley system with this penalty structure

