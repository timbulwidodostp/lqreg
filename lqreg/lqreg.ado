*! version 2.0.1 5aug16 N.Orsini, M.Bottai
*! version 2.0.0 23jun16 N.Orsini, M.Bottai
*! version 1.0.0 21jun12 N.Orsini, M.Bottai

capture program drop lqreg
program lqreg, eclass byable(recall) prop(sw mi or)
	if _by() {
		local BY `"by `_byvars'`_byrc0':"'
	}
	if replay() {
		if ("`BY'"!="") error 190
		if ("`e(cmd)'"!="lqreg") error 301
		Replay `0'
		exit
	}
	local vv : di "version " string(_caller()) ":"
	`vv' `BY' Estimate `0'
end

capture program drop Estimate
program Estimate, eclass byable(recall) sort prop(sw mi or)
	local vc = _caller()
	local cmdline : copy local 0
	version 13, missing
	local options "Level(cilevel)"
	
	syntax varlist(numeric fv) [if] [in] [fw iw pw] [, `options' Quantiles(numlist)  ///
			WLSiter(integer 1) Reps(string) SEED(string) vce(string) ///
			 noLOg noDOts Level(cilevel) ///
			GENerate(string) ymin(string) ymax(string) CLuster(varlist) COEF  * ]

	_get_diopts diopts options, `options'	
 	local fvops = ("`s(fvops)'"=="true" | `vc'>12)

	local vceqreg "`vce'"
	
		// check significance level  
 
		local opts "level(`level')"
		
		if "`cluster'" != "" & "`reps'" == "" {
			di _n as err "specify the option reps(#) for bootstrap"
			exit 198
		}
		
		if "`weight'" != "" & "`reps'" != "" {
				 di _n as err "weights not allowed with bootstrap"
				 exit 198
		}
		
		marksample touse
		qui count if `touse'
		
 		local N = r(N)
		if (!`N') error 2000
		if (`N'<=2) error 2001
	
		if "`log'"!="" | "`dots'"!="" {
			local log "*"
		}

 		SetQ `quantiles'
 		tokenize "`r(quants)'"
		local nq 1
		while "``nq''" != "" {
			local q`nq' ``nq''	
			local nq = `nq' + 1
		}
		local nq = `nq' - 1
		
		if "`reps'"!= "" { 
						  tempfile BOOTRES
						  capture confirm integer number `reps' 
							if _rc != 0 {
								di as err "Reps must be an integer number > 1"
								exit 198
							}
						    if (`reps')<2 { 
							di as err "Reps must be an integer number > 1"
								exit 198
							}
		}

		tempname coefs vce VCE coefs0 coefs1 handle
		tempvar e bcons bwt
 
		quietly count if `touse'
		if r(N)<4 { 
			di in red "insufficient observations"
			exit 2001
		}

	gettoken depv indep : varlist
	
	_rmcoll `indep' `weights' if `touse', expand
	local indep `r(varlist)' 
	local varlist `depv' `indep'
	qui _regress `depv' `indep' if `touse'
	local colna : colna e(b)
	local ncolna : list sizeof colna
 
	if ("`log'"!="") local qui quietly
	
	// determine unit (minimal increment of depvar) from codebook   
	
	tempname p 
		
	scalar `p' = 1
	
	capture assert float(`depv') == float(round(`depv',1)) if `touse'
	if _rc == 0 {
		while _rc == 0 {
			scalar `p' = `p'*10
			capture assert float(`depv') == float(round(`depv',`p')) if `touse'
		}
		scalar `p' = `p'/10
	}
	else {
		while _rc {
			scalar `p' = `p'/10
			capture assert float(`depv') == float(round(`depv',`p')) if `touse'
		}
	}
	
		// Transform the bounded outcome

		tempvar logity flag   
		
		tempname depmin depmax epsilon 
		
		qui su `depv'		
		
		// Epsilon = (minimal increment of the dependent variable)/2
		
		scalar `epsilon' = `p'/2
		 
         if "`ymin'" == "" scalar `depmin' = r(min)-`epsilon'
         else scalar `depmin' = `ymin'
							
		 if "`ymax'" == "" scalar `depmax' = r(max)+`epsilon'
		 else scalar `depmax' = `ymax'


	qui gen `logity' = log((`depv'-`depmin')/(`depmax'-`depv'))                         

	if "`generate'" != "" {		
		local nw : word count `generate'
		
		if `nw' != 1 {
				di as err "specify one variable name"
				exit 198
		}
		else {
				gen double `generate' = `logity'
				label var `generate' "Logit `depv'"
		}
	}
 
		forvalues k=1/`nq' {
 			local myeq = "q" + string(`q`k''*100,"%2.0f")
			local eqnames `eqnames' `myeq'
		}
		local eqnames : list uniq eqnames
		local k : word count `eqnames'

	if `k' < `nq' {
		di as err "only `k' out of `nq' quantiles are " /*
		 */ "unique within a relative precision of c(epsfloat)"
		exit 498
	}
	
	local EQNAMES : copy local eqnames
	forval k = 1/`nq' {
		tempname coef`k' vce`k'
		gettoken myeq EQNAMES : EQNAMES
		local COLNA : copy local colna
		forval i = 1/`ncolna' {
			gettoken name COLNA : COLNA
			local result "`result' (`coef`k''[1,`i'])"
			local eqnams "`eqnams' `myeq'"
			local conams "`conams' `name'"
			tempvar v
			local vl "`vl' `v'"
		}
	}
 

	if "`vceqreg'"	== "" qui qreg `logity' `indep' [`weight'`exp'] if `touse', `opts' q(`q1')  
   else qui qreg `logity' `indep' [`weight'`exp'] if `touse', `opts' q(`q1') vce(`vceqreg')

  
		local nobs `e(N)'
		local tdf `e(df_r)'
		local rsd1 `e(sum_rdev)'
		local msd1 `e(sum_adev)'
		local vle "_cons"
		local vli "`bcons'"

		local this_vce = "`e(vce)'"

		mat `coefs' = e(b)
        mat `vce1' = e(V)
		
		local k 2
		while `k' <= `nq' {
		
	 if ("`vceqreg'"	== "") qui qreg `logity'  `indep' [`weight'`exp'] if `touse', `opts' q(`q`k'')
	 else qui qreg `logity' `indep' [`weight'`exp'] if `touse', `opts' q(`q`k'') vce(`vceqreg')	

			if e(N) != `nobs' {
				di in red /*
	*/ "`q0' quantile:  `nobs' obs. used" _n /*
	*/ "`q`k'' quantile:  `e(N)' obs. used" _n /*
	*/ "Same sample cannot be used to estimate both quantiles." /*
	*/ "Sample size probably too small."
				exit 498
			}
			if e(df_r) != `tdf' {
				di in red /*
	*/ "`q0' quantile:  " `nobs'-`tdf' " coefs estimated" _n /*
	*/ "`q`k'' quantile:  " `e(N)'-`e(df_r)' coefs estimated" _n /*
	*/ "Same model cannot be used to estimate both quantiles." /*
	*/ "Sample size probably too small."
				exit 498
			}
			local msd`k' `e(sum_adev)'
			local rsd`k' `e(sum_rdev)'
			mat `coefs' = `coefs', e(b)
	
			mat  `vce`k'' = e(V)
		   
			local k = `k' + 1
		}
		
		mat colnames `coefs' = `conams'
		mat coleq `coefs' = `eqnams'

		// Asymptotic standard errors using qreg
	
				tempname coefs2
				mat `coefs2' = `coefs'
				
				//  Create the diagonal var/cov matrix from qreg output on each quantile
				
				mat `VCE' = `vce1' 
				
				local k 2
				while `k' <= `nq' {
					mata: st_matrix("`VCE'", blockdiag(st_matrix("`VCE'"),st_matrix("`vce`k''")) )		
					local k = `k' + 1	
				}
				
				mat rownames `VCE' = `conams'
				mat roweq `VCE' = `eqnams'
				mat colnames `VCE' = `conams'
				mat coleq `VCE' = `eqnams'
					
				ereturn  post `coefs2' `VCE', obs(`nobs') dof(`tdf') depn(`depv') 
						
		// Bootstrap confidence intervals
						
		if ("`reps'" != "") {
				
				 preserve
			     qui keep if `touse'				 
			     qui fvrevar `indep', list
 			     qui keep `depv' `logity' `r(varlist)' `cluster'
			
				qui gen double `bwt' = .
				`log' di in gr "(bootstrapping " _c
				qui postfile `handle' `vl' using "`BOOTRES'", double
				quietly noisily {
				
					// set the seed
					if "`seed'" != "" {
						set seed `seed'
					}
					local seed `c(seed)'
			
					local j 1
					while `j'<=`reps' {

					
						if "`cluster'" != "" {
										bsample , weight(`bwt') cluster(`cluster') 
						}
						else {
										bsample , weight(`bwt') 
						}
						
						capture noisily {
							local k 1
							while `k'<=`nq' {
								qreg_c `logity' `indep' , `opts' q(`q`k'') wvar(`bwt')
								
								mat `coef`k'' = e(b)
								local k =`k' + 1
							}
						}
						local rc = _rc
						if (`rc'==0) {
							post `handle' `result'
							`log' di in gr "." _c
							local j=`j'+1
						}
						else {
							if _rc == 1 { exit 1 }
							`log' di in gr "*" _c
						}
					}
				}
				local rc = _rc 
				postclose `handle'
				if `rc' { 
					exit `rc'
				}

				qui use "`BOOTRES'", clear

				quietly mat accum `VCE' = `vl', dev nocons
				mat rownames `VCE' = `conams'
				mat roweq `VCE' = `eqnams'
				mat colnames `VCE' = `conams'
				mat coleq `VCE' = `eqnams'
				mat `VCE'=`VCE'*(1/(`reps'-1))

				restore
				
				
				ereturn post `coefs' `VCE', obs(`nobs') dof(`tdf') depn(`depv')  
			
				`log' noi di in gr ")"
				
				capture erase "`BOOTRES'"
				ereturn scalar reps = `reps'
				ereturn local vcetype "Bootstrap"
				
				 
		} // End Boostrap 

	
	 ereturn repost, esample(`touse') buildfvinfo 
		
		local k 1
		while `k' <= `nq' {
			local rounded : di %3.2f `q`k''
			ereturn scalar q`k' = `q`k''
			local k = `k' + 1
		}
		
	ereturn local depvar "`depv'"
	ereturn scalar N = `nobs'
	ereturn  scalar df_r = `tdf'
	ereturn  scalar ymin = `depmin'
	ereturn  scalar ymax = `depmax'
	if "`reps'" != "" {
			ereturn  local reps `reps'
			ereturn  local seed `seed'
	}	

	ereturn local cmd "lqreg"
	ereturn local depvar "`depv'"
	ereturn scalar n_q = `nq'
	ereturn local eqnames "`eqnames'"
	ereturn scalar convcode = 0
	ereturn local marginsnotok stdp stddp Residuals
	ereturn local predict "sqreg_p"
	ereturn local cmdline `"lqreg `cmdline'"'
	
	 if "`reps'" == ""  ereturn local vce = "`this_vce'"
     if "`this_vce'" == "robust" ereturn local vcetype = "Robust"
     if "`reps'" != "" ereturn local vcetype = "Bootstrap"

    Replay, `diopts' level(`level')  bs(`reps') `coef'
    error `e(convcode)'

end

capture program drop SetQ
program define SetQ /* <nothing> | # [,] # ... */ , rclass
	if "`*'"=="" {
		ret local quants ".5"
		exit
	}
	local orig "`*'"

	tokenize "`*'", parse(" ,")

	while "`1'" != "" {
		FixNumb "`orig'" `1'
		ret local quants "`return(quants)' `r(q)'"
		mac shift 
		if "`1'"=="," {
			mac shift
		}
	}
end

capture program drop FixNumb
program define FixNumb /* # */ , rclass
	local orig "`1'"
	mac shift
	capture confirm number `1'
	if _rc {
		Invalid "`orig'" "`1' not a number"
	}
	if `1' >= 1 {
		ret local q = `1'/100
	}
	else 	ret local q `1'
	if `return(q)'<=0 | `return(q)'>=1 {
		Invalid "`orig'" "`return(q)' out of range"
	}
end
		
capture program drop Invalid
program define Invalid /* "<orig>" "<extra>" */
	di in red "quantiles(`1') invalid"
	if "`2'" != "" {
		di in red "`2'"
	}
	exit 198
end

capture program drop PrForm
program define PrForm /* # */ , rclass
	local aa : di %8.2f `1'
	ret local pr `aa'
	if substr("`return(pr)'",1,1)=="0" {
		ret local pr = substr("`return(pr)'",2,.)
	}
end

capture program drop Replay
program define Replay
	syntax [, bs(string) COEF  * ]

		di _n in gr "Logistic Quantile Regression" _col(54) _c
	di in gr "Number of obs =" in ye %10.0g e(N)
	di in gr "Bounded Outcome: `e(depvar)'" in gr "(" in y `e(ymin)' in gr ", " in y  `e(ymax)' in gr ")" _col(54) _c
	if ("`bs'" != "") {
								di in gr "Bootstrap(" in y `e(reps)' in gr ") SEs"
	}
	else {
			di in gr "Asymptotic SEs"  
	}
	
	_get_diopts diopts options, `options'	

 	if "`coef'" == "" ereturn display,  `diopts' eform("Odds Ratio")
     else ereturn display, `diopts'
	 
	error e(convcode)
end

exit

 sysuse auto, clear
 set trace off
  lqreg mpg i.foreign,       cformat(%9.2f)
  exit
  eret list
  lqreg mpg i.foreign,   q(25 50)   vce(iid)
  eret list
    
  lqreg mpg i.foreign,   q(25 50)   vce(robust)
  eret list
exit
 
* Some examples from the SJ paper

use http://nicolaorsini.altervista.org/data/pa_luts, clear
gen cage = age-59
lqreg ipss i.tpac  
lqreg ipss i.tpac c.cage, cformat(%9.2f) q(25 50 75)
