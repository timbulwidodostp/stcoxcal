*! version 1.0.1 PR 27nov2014
program define stcoxcal, rclass
version 13.0
// Time-dependent calibration plot and tests for PI from Cox model
syntax varlist(min=1 max=1 numeric) [if] [in] , TImes(numlist >0) ///
 [ noGRaph RESiduals SAVing(string) test TRend val(varname) * ]
quietly {
	tempname ests
	capture _estimates hold `ests' // in case a user model has been fit recently
	local tounhold = (c(rc) == 0)
	// Get H0 from PI
	marksample touse
	markout `touse' `val'
	preserve
	keep if `touse'==1
	drop `touse'
	tempvar xb
	if "`val'"!="" {
		local ifval1 if `val'!=0 & `val'<.
		count `ifval1'
		if r(N) < 10 {
			di as err "insufficients observations (" r(N) ") of `varlist' in the validation dataset"
			exit 2001
		}
		sum `varlist' if `val'==0, meanonly
		if r(N) < 10 {
			di as err "insufficients observations (" r(N) ") of `varlist' in the derivation dataset"
			exit 2001
		}
		gen double `xb' = `varlist' - r(mean)
		stcox if `val'==0, offset(`xb') estimate
	}
	else {
		local ifval1 if 1
		sum `varlist', meanonly
		if r(N) < 10 {
			di as err "insufficients observations (" r(N) ") of `varlist'"
			exit 2001
		}
		gen double `xb' = `varlist' - r(mean)
		stcox, offset(`xb') estimate
	}
	Drop _lnH0
	predict _lnH0, basechaz
	replace _lnH0 = ln(_lnH0)

	// Find best fit FP2 powers to log cumulative hazard function
	fp <_t>, replace all: regress _lnH0 <_t>
	local b0 = _b[_cons]
	local b1 = _b[_t_1]
	local b2 = _b[_t_2]
	tempname fpfp
	matrix `fpfp' = e(fp_fp)
	local p1 = `fpfp'[1, 1]
	local p2 = `fpfp'[1, 2]
	drop _lnH0

	// Predict event probs on validation data at selected times.
	local j 0
	foreach time of local times {
		local ++j
		Drop _clogF`j'
		Drop _F`j'
		fraceval num `time' "`p1' `p2'" "`b1' `b2'" `b0'
		gen double _clogF`j' = `xb' + r(fp) `ifval1' // predicted log cumulative hazard at `time'
		gen double _F`j' = 1 - exp(-exp(_clogF`j'))
	}

	// Predict pseudovalues on validation
	stpsurv `ifval1', at(`times') gen(_f) failure
	local nt : word count `times'
	if `nt' < 2 rename _f _f1

	// Reshape data in validation subset for analysis and plotting
	keep `ifval1'
	Drop _id
	Drop _times
	gen long _id = _n
	reshape long _f _F _clogF, i(_id) j(_times)
	if "`test'`trend'" != "" {
		// Estimate overall calibration slope and intercept on pseudovalues
		local glmopt link(cloglog) vce(cluster _id) irls noheader nolog

		// Test constants with slope constrained to 1
		capture noisily glm _f ibn._times, noconstant offset(_clogF) `glmopt'
		if c(rc)==0 {
			noi di as txt _n "[Test 1: intercepts (gamma0) = 0 with slope (gamma1) constrained to 1]"
			noi testparm i._times
			local P0 = chi2tail(r(df), r(chi2))
		}
		else {
			di as err "[could not fit constrained GLM to estimate constant on pseudovalues]"
		}
		capture noisily glm _f _clogF ibn._times, noconstant `glmopt'
		if c(rc)==0 {
			noi di as txt _n "[Test 2: slope (gamma1) = 1 with constants (gamma0) estimated]"
			noi test _clogF = 1
			local P1 = chi2tail(1, r(chi2))
			local gamma1 = _b[_clogF]
			local gamma1_se = _se[_clogF]
			noi di as txt _n "[Test 3: joint test of slope (gamma1) = 1 and all constants (gamma0) = 0]"
			testparm i._times
			noi test _clogF = 1, accum
			local P01 = chi2tail(r(df), r(chi2))
		}
		else {
			di as err "[could not fit GLM to estimate constant and slope on pseudovalues]"
		}
		// Checking for change in calibration slope with time point (interaction)
		if `nt' > 1 {
			if "`trend'" == "" {
				capture noisily glm _f i._times##c._clogF, `glmopt'
				if c(rc)==0 {
					noi di as txt _n "[Test 4: interaction between slopes (gamma1) and times]"
					noi testparm _times#c._clogF
					local Pint = chi2tail(r(df), r(chi2))
				}
				else {
					di as err "[could not fit GLM to estimate slope x time interaction]"
				}
			}
			else {
				capture noisily glm _f c._times##c._clogF, `glmopt'
				if c(rc)==0 {
					noi di as txt _n "[Test 4: interaction between slope (gamma1) and scores for times]"
					noi test c._times#c._clogF
					local Pint = chi2tail(r(df), r(chi2))
				}
				else {
					di as err "[could not fit GLM to estimate slope x time trend]"
				}
			}
			tempvar fitted
			predict `fitted'
		}
	}
	if "`graph'" != "nograph" {
		if "`residuals'" != "" {
			tempvar diff
			gen `diff' = _f - _F
			if "`fitted'" != "" replace `fitted' = `fitted' - _F
		}
		// Running line plots of pseudovalues on predicted event probs
		local j 0
		local gs
		foreach time of local times {
			local ++j
			if "`residuals'" != "" {
				if "`fitted'" == "" {
					running `diff' _F if _times==`j', ///
					 title("t = `time'") ci leg(off) name(g`j',replace) xtitle("") ytitle("") ///
					 xlabel(0(.25)1) yline(0) nopts nodraw `options'
				}
				else {
					running `diff' _F if _times==`j', ///
					 addplot(line `fitted' _F if _times==`j', sort lp(-) lwidth(medthick ..)) ///
					 title("t = `time'") ci leg(off) name(g`j',replace) xtitle("") ytitle("") ///
					 xlabel(0(.25)1) yline(0) nopts nodraw `options'
				}
			}
			else {
				running _f _F if _times==`j', ///
				 addplot(line _F `fitted' _F if _times==`j', sort lp(l -) lwidth(medthick ..)) ///
				 title("t = `time'") ci leg(off) name(g`j',replace) xtitle("") ytitle("") ///
				 xlabel(0(.25)1) nopts nodraw `options'
			}				
			local gs `gs' g`j'
		}
		if "`residuals'" != "" local ltitle "Observed minus predicted event probability"
		else local ltitle "Observed event probability"
		graph combine `gs', imargin(small) b2title("Predicted event probability") ///
		 l1title("`ltitle'") xcommon ycommon name(_g, replace)
	}
	if `"`saving'"' != "" {
		if "`fitted'" != "" {
			cap drop _fitted
			rename `fitted' _fitted
			local fitted _fitted
		}
		keep _f _F `fitted' _clogF _times _id
		save `"`saving'"', replace
	}
}
return local fp_pwrs `p1' `p2'
foreach thing in Pint P01 P1 P0 gamma1_se gamma1 {
	if "``thing''" != "" return scalar `thing' = ``thing''
}
restore
if `tounhold' _estimates unhold `ests'
end

program define Drop
args var
capture confirm var `var', exact
if c(rc)==0 drop `var'
end

* v 2.2.0 PR 09aug2011.
program define fraceval, rclass
	version 10.0
	args t X pwrs betas beta0 adjust sh sc	/* t=var | num */
	if "`9'"!="" {
		mac shift 6
		di as err "invalid `*'"
		exit 198
	}
	// Form of command: 1=FP stuff stored by fracpoly in e(), 0=get FP info from user input
	local stored=("`pwrs'"=="")
	if `stored' {
		if "`e(cmd)'" == "" | "`e(fp_cmd2)'" != "fracpoly" error 301
		if "`betas'`beta0'`adjust'`sh'`sc'"!="" {
			di as err "invalid first form of -fraceval- command"
			exit 198
		}
	}
	else {
		if "`adjust'"!="" 	confirm num `adjust'
		if "`sh'"!="" 		confirm num `sh'
		if "`sc'"!=""		confirm num `sc'
	}
	if "`t'"!="var" & "`t'"!="num" {
		di as err "invalid `t'"
		exit 198
	}
	if "`t'"=="var" {
		confirm var `X'
		tempvar x
		qui gen double `x'=`X'
		local gen qui gen double
		local replace qui replace
		local temp tempvar
	}
	else {
		tempname x
		cap scalar `x'=`X'
		if _rc {
			di as err "invalid " `X'
			exit 198
		}
		local gen scalar
		local replace scalar
		local temp tempname
	}
	tempname shift scale expx
	scalar `expx'=.
	if `stored' {
/*
	Get stored FP quantities. Calculate mean of X.
*/
		local np: word count `e(fp_xp)'
		tempname b0
		cap scalar `b0'=_b[_cons]
		if _rc scalar `b0'=0
		scalar `shift'=e(fp_shft)
		scalar `scale'=e(fp_sfac)
		if "`e(fp_xpx)'"!="" scalar `expx'=`e(fp_xpx)'
		// Check for adjust() option in cmdline
		local 0 `"`e(cmdline)'"'
		syntax [anything] [if] [in] [fweight  aweight  pweight  iweight] [, ADjust(string) *]
		if "`adjust'"=="no" local adjust ""
		else if "`adjust'"=="mean" | "`adjust'"=="" {
			// compute mean of original x in estimation sample
			qui sum `e(fp_x1)' if e(sample)
			local adjust = r(mean)
		}
	}
	else {
/*
	Quantities are in `pwrs', `beta0' and `betas'
*/
		local np: word count `pwrs'
		local nb: word count `betas'
		if `np'!=`nb' {
			di as err "numbers of powers and betas differ"
			exit 198
		}
		if "`beta0'"!="" {
			cap confirm num `beta0'
			if _rc==0 | (_rc>0 & "`beta0'"==".") {
				tempname b0
				scalar `b0'=`beta0'
			}
			else {
				confirm var `beta0'
				local b0 `beta0'
			}
		}
		else {
			tempname b0
			scalar `b0'=0
		}
		if "`sh'"=="" scalar `shift'=0
		else scalar `shift'=`sh'
		if "`sc'"=="" scalar `scale'=1
		else scalar `scale'=`sc'
	}
	forvalues i=1/`np' {
		if `stored' {
			tempname p`i' b`i'
			local w: word `i' of `e(fp_k1)'
			scalar `p`i''=`w'
			local w: word `i' of `e(fp_xp)'
			scalar `b`i''=_b[`w']
		}
		else {
			local w: word `i' of `pwrs'
			cap confirm num `w'
			if _rc==0 | (_rc>0 & "`w'"==".") {
				tempname p`i'
				scalar `p`i''=`w'
			}
			else {
				confirm var `w'
				if "`t'"!="var" {
					di as err "`w' invalid---" /*
					 */ "`w' is a var, `X' is a number"
					exit 198
				}
				local p`i' `w'
			}
			local w: word `i' of `betas'
			cap confirm num `w'
			if _rc==0 | (_rc>0 & "`w'"==".") {
				tempname b`i'
				scalar `b`i''=`w'
			}
			else {
				confirm var `w'
				if "`t'"!="var" {
					di as err "`w' invalid---" /*
					 */ "`w' is a var, `X' is a number"
					exit 198
				}
				local b`i' `w'
			}
		}
	}
/*
	Compute FP function.
*/
	tempname small
	`temp' fp h hlast lnx plast
	scalar `small' = 1e-6
	`replace' `x' = (`x' + `shift') / `scale'
	if (`expx' != .) `replace' `x' = exp(`expx' * `x')
	`gen' `lnx' = log(`x')
	`gen' `h' = .
	`gen' `hlast' = 1
	`gen' `plast' = 0
	if "`adjust'" != "" {
		tempname adj
		scalar `adj' = (`adjust' + `shift') / `scale'
		if (`expx' != .) scalar `adj' = exp(`expx' * `adj')
		`temp' a alast lna
		`gen' `lna' = log(`adj')
		`gen' `a' = .
		`gen' `alast' = 1
	}
	else {
		tempname a
		scalar `a' = 0
	}
	forvalues j = 1 / `np' {
		if ("`replace'" == "qui replace") local ifpj if `p`j''!=.
		`replace' `h' = cond(abs(`p`j'' - `plast') < `small', `lnx' * `hlast', ///
		  cond(abs(`p`j'') < `small', `lnx', ///
		  cond(abs(`p`j'' - 1) < `small', `x', `x'^`p`j''))) `ifpj'
		if "`adjust'" != "" {
			`replace' `a' = cond(abs(`p`j'' - `plast') < `small', `lna' * `alast', ///
			  cond(abs(`p`j'') < `small', `lna', ///
			  cond(abs(`p`j'' - 1) < `small', `adj', `adj'^`p`j''))) `ifpj'
			`replace' `alast' = `a'
		}
		// If appropriate, extract adjustment
		if `stored' {
			if !missing(e(fp_a`j')) scalar `a' = e(fp_a`j')
		}
		if (`j' == 1) `gen' `fp' = `b0' + `b`j'' * (`h' - `a')
		else `replace' `fp' = `fp' + `b`j'' * (`h' - `a')
		`replace' `hlast' = `h'
		`replace' `plast' = `p`j''
	}
	return local fpvar _fp
	if "`t'"=="num" {
		di as txt "FP(" as res `X' as txt ") = " /*
		 */ as res %9.0g `fp'
		return scalar fp = `fp'
	}
	else {
		cap drop _fp
		rename `fp' _fp
		di as txt "_fp created from `X'"
		return scalar fp = _fp[1]
	}
end
