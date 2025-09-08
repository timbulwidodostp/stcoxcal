*************************************************************************
* Stpci.
*
* Pseudo observations for the cumulative incidence function under 
* competing risks.
*
* Author: Erik Parner, Department of Biostatistics, University of Aarhus.
*         Email: parner@biostat.au.dk
* 
* 14August2010.
*************************************************************************

program stpci
	version 11
	st_is 2 analysis
	
	syntax varlist(min=1 max=1 numeric) [if] [in] , At(numlist>0 min=1 ascending) [Generate(string)]

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
	
	* Checking the variable containing competing risk events.
	tempvar compvar
	quietly gen `compvar'=`varlist'
	capture confirm numeric variable `compvar'
	if(_rc) {
	 	display as error "The competing risk variable `compet' is not numeric."
	 	exit
	}
	quietly count if( `compvar'==1 & _d==1 )
	local count=r(N)
	if(`count'>0) {
		display as error "The competing risk variable marks events as being competing risk."
		exit
	} 
	
	* Start the progress indication.
	nois _dots 0, title(Computing pseudo observations (progress dots indicate percent completed))

	tempname w d
	if(`"`fw'"'!="") {
		gen `w'=`fw'
	}
	else {
		gen `w'=1
	}
	gen `d'=0*(_d==0 & `compvar'!=1)+1*(_d==1)+2*(_d==0 & `compvar'==1)
	local name="pseudo"
	if(`"`generate'"'!="") {
		local name="`generate'"
	} 
	mata: pci("_t","`d'","`w'","`touse'","`tmax'","`name'")
end  



version  11
mata:
void  pci(string  scalar  tname,  string  scalar  dname, string scalar wname, string scalar tousename,string tmaxname,string scalar var)
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
	cs=J(length(ts),1,0)
	ds1=J(length(ts),1,0)
	ds2=J(length(ts),1,0)
	for (i=1; i<=length(ts); i++) {
		for (j=1; j<=length(t); j++) {
			if(ts[i]==t[j] & d[j]==0 & touse[j]==1) {
				cs[i]=cs[i]+w[j]
			}		
			if(ts[i]==t[j] & d[j]==1 & touse[j]==1) {
				ds1[i]=ds1[i]+w[j]
			}
			if(ts[i]==t[j] & d[j]==2 & touse[j]==1) {
				ds2[i]=ds2[i]+w[j]
			}
		}
	}
	
	// The total cip.
	noobs=colsum(cs+ds1+ds2)
	Ytot=J(length(ts),1,0)
	Stot=J(length(ts),1,0)
	Stotlag=J(length(ts),1,0)
	Ytot[1]=noobs
	Stot[1]=1-(ds1[1]+ds2[1])/Ytot[1]
	Stotlag[1]=1
	for (i=2; i<=length(ts); i++) {
		Ytot[i]=Ytot[i-1]-ds1[i-1]-ds2[i-1]-cs[i-1]
		Stot[i]=Stot[i-1]*(1-(ds1[i]+ds2[i])/Ytot[i])
		Stotlag[i]=Stot[i-1]
	}
	citot=runningsum(Stotlag:*(ds1:/Ytot))	
	citottmax=J(1,ntmax,0)
	for(i=1; i<=ntmax; i++) {
		citottmax[i]=max(select(citot,(ts:<=tmax[i])))
		if(citottmax[i]==.) citottmax[i]=0
	}
	
	// The pseudo values.
	pseudocs=J(length(ts),ntmax,0)	
	pseudods1=J(length(ts),ntmax,0)
	pseudods2=J(length(ts),ntmax,0)
	S=J(length(ts),1,0)
	Slag=J(length(ts),1,0)
	obsno=0
	lastobsno=0
	for(i=1; i<=length(ts); i++) {
		for(j=0; j<=2; j++) {
			if( (j==0 & cs[i]>0) | (j==1 & ds1[i]>0) | (j==2 & ds2[i]>0) ) {
				if(j==0)      cs[i]=cs[i]-1
				else if(j==1) ds1[i]=ds1[i]-1
				else          ds2[i]=ds2[i]-1		
				Y=Ytot
				for(k=1; k<=i; k++) {
					Y[k]=Y[k]-1
				}
				S[1]=1-(ds1[1]+ds2[1])/Y[1]
				Slag[1]=1
				for (k=2; k<=length(ts); k++) {
					S[k]=S[k-1]*(1-(ds1[k]+ds2[k])/Y[k])
					Slag[k]=S[k-1]
				}
				ci=runningsum(Slag:*(ds1:/Y))
				citmax=J(1,ntmax,0)
				for(k=1; k<=ntmax; k++) {
					citmax[k]=max(select(ci,(ts:<=tmax[k])))
					if(citmax[k]==.) citmax[k]=0
					if(j==0)      pseudocs[i,k]=noobs*citottmax[k]-(noobs-1)*citmax[k]
					else if(j==1) pseudods1[i,k]=noobs*citottmax[k]-(noobs-1)*citmax[k]
					else          pseudods2[i,k]=noobs*citottmax[k]-(noobs-1)*citmax[k]
				}
				if(j==0)      cs[i]=cs[i]+1
				else if(j==1) ds1[i]=ds1[i]+1
				else   		  ds2[i]=ds2[i]+1 
				
				// Progress indication.
				lastobsno=obsno
				if(j==0)      obsno=obsno+cs[i]
				else if(j==1) obsno=obsno+ds1[i]
				else          obsno=obsno+ds2[i]

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
					if(d[i]==0)      pseudo[i,k]=pseudocs[j,k]
					else if(d[i]==1) pseudo[i,k]=pseudods1[j,k]
					else 			 pseudo[i,k]=pseudods2[j,k]
				}
			}
		}
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
