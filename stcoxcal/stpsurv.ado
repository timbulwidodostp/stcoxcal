*************************************************************************
* Stpsurv.
*
* Pseudo observations for the survival function.
*
* Author: Erik Parner, Department of Biostatistics, University of Aarhus.
*         Email: parner@biostat.au.dk
* 
* 14August2010.
*************************************************************************

program stpsurv
	version 11
	st_is 2 analysis
	
	syntax [if] [in] , At(numlist>0 min=1 ascending) [Generate(string) Failure ]

	* Weights information.
	local w: char _dta[st_wv]
	local fw `"`w'"'
	
	* Mark the sample.
	marksample touse
	quietly replace `touse' = 0 if _st==0
	
	* Check if some event time are left truncated.
	quietly sum _t0 if(`touse')
	local max=r(max)
	if(`max'>0) {
		display as error "The function does not support left truncated event times."
		exit
	}
	
	* Defining and checking tmax.
	quietly sum _t if(`touse')
	local max=r(max)
	foreach time of numlist `at' {
		if(`time'>`max') {
			display as error "The time point `time' is greater than the largest event time."
			exit
		}
	}
	numlist "`at'"
	local tmax=r(numlist)
	
	* Start the progress indication.
	nois _dots 0, title(Computing pseudo observations (progress dots indicate percent completed))

	tempname w
	if(`"`fw'"'!="") {
		gen `w'=`fw'
	}
	else {
		gen `w'=1
	}
	local name="pseudo"
	if(`"`generate'"'!="") {
		local name="`generate'"
	}
	local type="survival"
	if(`"`failure'"'!="") {
		local type="failure"
	}
	mata: psurv("_t","_d","`w'","`touse'","`tmax'","`name'","`type'")
end  


version  11
mata:
void  psurv(string  scalar  tname, string scalar dname, string scalar wname, string scalar tousename,string tmaxname,string scalar var,string scalar type)
{
	t=st_data(.,  tname)
	d=st_data(.,  dname)
	w=st_data(.,  wname)
	touse=st_data(., tousename)
	tmax=strtoreal(tokens(tmaxname))
	ntmax=length(tmax)
	
	// Collapsing the data.
	ts=select(t,(touse:==1))
	ts=uniqrows(ts)
	ts=sort(ts,1)
	ds=J(length(ts),1,0)
	cs=J(length(ts),1,0)
	for (i=1; i<=length(ts); i++) {
		for (j=1; j<=length(t); j++) {
			if(ts[i]==t[j] & d[j]==1 & touse[j]==1) {
				ds[i]=ds[i]+w[j]
			}
			if(ts[i]==t[j] & d[j]==0 & touse[j]==1) {
				cs[i]=cs[i]+w[j]
			}
		}
	}

	// The total KM.
	noobs=colsum(cs+ds)
	Ytot=J(length(ts),1,0)
	Stot=J(length(ts),1,0)
	Ytot[1]=noobs
	Stot[1]=1-ds[1]/Ytot[1]
	for (i=2; i<=length(ts); i++) {
		Ytot[i]=Ytot[i-1]-ds[i-1]-cs[i-1]
		Stot[i]=Stot[i-1]*(1-ds[i]/Ytot[i])
	}
	Stottmax=J(1,ntmax,0)
	for(i=1; i<=ntmax; i++) {
		Stottmax[i]=min(select(Stot,(ts:<=tmax[i])))
		if(Stottmax[i]==.) Stottmax[i]=1
	}

	// The pseudo values.
	pseudods=J(length(ts),ntmax,0)
	pseudocs=J(length(ts),ntmax,0)
	S=J(length(ts),1,0)
	obsno=0
	lastobsno=0
	for(i=1; i<=length(ts); i++) {
		for(j=0; j<=1; j++) {
			if( (j==0 & cs[i]>0) | (j==1 & ds[i]>0) ) {
				if(j==0) cs[i]=cs[i]-1
				else     ds[i]=ds[i]-1			
				Y=Ytot
				for(k=1; k<=i; k++) {
					Y[k]=Y[k]-1
				}
				S[1]=1-ds[1]/Y[1]
				for(k=2; k<=length(ts); k++) {
					S[k]=S[k-1]*(1-ds[k]/Y[k])
				}
				Stmax=J(ntmax,1,0)
				for(k=1; k<=ntmax; k++) {
					Stmax[k]=min(select(S,(ts:<=tmax[k])))
					if(Stmax[k]==.) Stmax[k]=1
					if(j==0) pseudocs[i,k]=noobs*Stottmax[k]-(noobs-1)*Stmax[k]
					else     pseudods[i,k]=noobs*Stottmax[k]-(noobs-1)*Stmax[k]
				}
				if(j==0) cs[i]=cs[i]+1 
				else     ds[i]=ds[i]+1
				
				// Progress indication.
				lastobsno=obsno
				if(j==0) obsno=obsno+cs[i]
				else obsno=obsno+ds[i]
				if( floor( (obsno*100)/noobs ) > floor( (lastobsno*100)/noobs ) ) {
					number=floor( (lastobsno*100)/noobs)
					cc=floor( (obsno*100)/noobs ) - floor( (lastobsno*100)/noobs )
					for (k=1; k<=cc; k++) {
						stata("noisily _dots "+strofreal(number+k)+" 0")
					}
				}
			}
		}
	}	
	
	// Expanding the pseudo values.
	pseudo=J(length(t),ntmax,0)
	_fillmissing(pseudo)
	for (i=1; i<=length(t); i++) {
		for (j=1; j<=length(ts); j++) {
			if(touse[i]==1 & t[i]==ts[j]) {
				for (k=1; k<=ntmax; k++) {
					if(d[i]==0) pseudo[i,k]=pseudocs[j,k]
					else        pseudo[i,k]=pseudods[j,k]
				}
			}
		}
	}
	
	// Returning the pseudo values.
	if(type=="failure") {
		pseudo=1:-pseudo
	}
	if(ntmax==1) {
		ids=st_addvar("double",var)
		st_varlabel(ids,"Pseudo values at time "+strofreal(tmax[1])+".")
		st_store(.,ids,pseudo)
		stata(`"disp as result "Generated pseudo variable: "'+var+`"""')
	}
	else {
		for (k=1; k<=ntmax; k++) {
			name=var+strofreal(k)
			ids=st_addvar("double",name)
			st_varlabel(ids,"Pseudo values at time "+strofreal(tmax[k])+".")
			st_store(.,ids,pseudo[,k])
		}
		stata(`"disp as result "Generated pseudo variables: "'+var+strofreal(1)+"-"+var+strofreal(ntmax)+`"""')
	}
}
end
