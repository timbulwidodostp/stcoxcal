{smcl}
{* *! version 1.0.1 PR 27nov2014}{...}
{cmd:help stcoxcal}{right: ({browse "http://www.stata-journal.com/article.html?article=st0357":SJ14-4: st0357})}
{hline}


{title:Title}

{p2colset 5 17 20 2}{...}
{p2col :{hi:stcoxcal} {hline 2}}Calibration plots and tests for Cox proportional hazards model{p_end}
{p2colreset}{...}


{title:Syntax}

{p 8 12 2}
{cmd:stcoxcal}
{it:xbetavar}
{ifin}{cmd:,}
{cmdab:ti:mes(}{it:numlist}{cmd:)}
[{it:options}]


{marker stcoxcal_options}{...}
{synoptset 24 tabbed}{...}
{synopthdr:options}
{synoptline}
{p2coldent:* {cmdab:ti:mes(}{it:numlist}{cmd:)}}list time points at which to assess calibration{p_end}
{synopt :{cmdab:nogr:aph}}suppress graphs{p_end}
{synopt :{cmdab:res:iduals}}plot smoothed residuals (observed minus predicted event probabilities){p_end}
{synopt :{cmdab:sav:ing(}{it:filename}{cmd:)}}save created variables to file {it:filename}{p_end}
{synopt :{opt test}}test the overall calibration slope and factorial interaction of calibration slope with time{p_end}
{synopt :{cmdab:tr:end}}test the overall calibration slope and linear interaction of calibration slope with integer scores for time{p_end}
{synopt :{opt val(varname)}}evaluate calibration in independent data indicated by {it:varname} ~= 0{p_end}
{synopt :{it:graph_twoway_options}}options for {cmd:graph twoway}{p_end}
{synopt :{it:running_options}}options for {cmd:running}{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
* {opt times(numlist)} is required.{p_end}

{pstd}
You must have {cmd:running} and {cmd:stpsurv} installed before using
{cmd:stcoxcal}. {cmd:running} is available on the Statistical Software
Components archive (see {helpb ssc}).  {cmd:stpsurv} can be installed from the
{it:Stata Journal} archive via the command
{cmd:net install st0202, from(http://www.stata-journal.com/software/sj10-3)}.


{title:Description}

{pstd}
{cmd:stcoxcal} is a tool for examining the (possibly time-dependent) calibration
of a Cox model whose linear predictor ("prognostic index") is supplied in
{it:xbetavar}.

{pstd}
A "well-calibrated" model is one that accurately predicts survival or
event probabilities at all relevant follow-up times.
A model that includes covariates whose effects change (for example,
dwindle) over time is unlikely to be well calibrated.
Such a model will give a more or less biased prediction of survival
probabilities.  {cmd:stcoxcal} is designed to detect and display the lack of
calibration graphically.  It also includes tests of good calibration and
of time-dependent trends of miscalibration -- see Royston (Forthcoming) for
further information.

{pstd}
By default, {cmd:stcoxcal} examines calibration of a model on its "own" dataset.
With the {opt val()} option, {cmd:stcoxcal} can be used to examine model
calibration in an independent dataset (that is, for external validation).

{pstd}
Details of the methodology underlying {cmd:stcoxcal} and examples in real data
are given by Royston (Forthcoming).


{title:Options}

{phang}
{opt times(numlist)} lists times at
which model calibration is to be assessed.
{cmd:times()} is required.

{phang}
{opt nograph} suppresses the production of calibration plots.

{phang}
{opt residuals} plots the smoothed residuals (difference between observed
and predicted event probabilities) against the predicted event probabilities.
The default is to plot smoothed observed against predicted event
probabilities.

{phang}
{opt saving(filename)} saves five variables in the validation dataset
to file {it:filename}:

        {cmd:_id}    observation number in the original data
        {cmd:_times} integer scores (levels) 1, 2, ... of times specified in {opt times()}
        {cmd:_f}     pseudovalues for event probabilities
        {cmd:_F}     event probabilities predicted from a Cox model
        {cmd:_clogF} complementary log-log transformation of {cmd:_F}
    
{pmore}
These variables can be used by an expert to create plots
and to further analyze model calibration.  The data are in long
format, with a complete set of values for each level of {cmd:_times}.

{phang}
{opt test} tests whether the slope (on the log cumulative-hazard scale)
of the regression of pseudovalues for event probabilities on predicted
event probabilities over all time points in {opt times()} equals one.  A
nonsignificant p-value suggests good overall calibration, sometimes
called "calibration in the large".  {opt test} also tests the
interaction between the slopes and the times specified in {opt times()}.
A significant p-value suggests that calibration changes over time.
Typically, calibration declines as follow-up time increases.

{phang}
{opt trend} tests whether the slope (on the log cumulative-hazard scale)
of the regression of pseudovalues for event probabilities on predicted
event probabilities over all time points in {opt times()} equals one
(same as for {opt test}).  {opt trend} also tests the linear
interaction between the slopes and the integer scores for the times
specified in {opt times()}.  This may be more powerful than the interaction
test provided by {opt test}.

{phang}
{opt val(varname)} is for use in external validation.  {it:varname} is a
binary variable coded zero to define the "model derivation" dataset and
any other nonmissing value to define the "model validation" dataset.
Predictions of event probabilities at different times from the
derivation dataset are made in the validation dataset via the linear
predictor and a smoothed version of the baseline cumulative-hazard
function in the derivation dataset.  Royston and Altman (2013) call this
"strict" calibration.

{phang}
{it:graph_twoway_options} are options of {cmd:graph twoway}.  These may
be used to customize the appearance of the calibration plots.

{phang}
{it:running_options} are options of {cmd:running}.  These may be
used to customize the smoothing of pseudovalues.  The most relevant
option is likely to be {opt span(#)}.  See help on {helpb running} for further
information.


{title:Remarks}

{pstd}
Note that {cmd:stcoxcal} computes the baseline survival and cumulative
hazard functions internally.  As a preliminary, {cmd:stcoxcal} centers the
prognostic index supplied in {it:xbetavar} on zero.  If {opt val(varname)} is
provided, the mean of {it:xbetavar} in the subset defined by {it:varname} = 0
is subtracted from all values of {it:xbetavar}.  Otherwise, centering takes
place over the estimation sample.  Next, a Cox model is fit with no
covariates and with {it:xbetavar} offset from the linear predictor.  Again,
this is done either in the {it:varname} = 0 subset or in the estimation
sample.  Finally, the baseline cumulative hazard function is predicted and
smoothed for use with the calibration method described by Royston
(Forthcoming).

{pstd}
Because {it:xbetavar} or indeed the original covariates are not refitted to
the validation data, {cmd:stcoxcal} can be used in "partial validation" mode.
The prognostic index is created from a derivation model fit elsewhere
and imported for application in the available validation dataset.  Validation
is partial because the baseline cumulative hazard and survival functions are
estimated by {cmd:stcoxcal} on the validation data, whereas {it:xbetavar} is
calculated by the user on the validation data from regression coefficients
estimated externally.  Although imperfect, partial validation
nevertheless allows a useful evaluation of the predictive accuracy of a
predefined model when the baseline distribution function is (perforce)
tailored to the validation data.


{title:Examples}

{phang}{cmd:. webuse brcancer}{p_end}
{phang}{cmd:. stset rectime, failure(censrec) scale(365.24)}{p_end}
{phang}{cmd:. fp generate x1^(-2 -0.5)}{p_end}
{phang}{cmd:. fp generate x6^(0.5), scale}{p_end}
{phang}{cmd:. stcox x1_1 x1_2 x4a x4b x5e x6_1 hormon}{p_end}
{phang}{cmd:. predict xb, xb}{p_end}
{phang}{cmd:. stcoxcal xb, times(1(1)6) test}{p_end}
{phang}{cmd:. stcoxcal xb, times(1(1)6) trend}{p_end}

{phang}{cmd:. set seed 3143}{p_end}
{phang}{cmd:. gen byte random_half = (runiform() < 0.5)}{p_end}
{phang}{cmd:. stcox x1_1 x1_2 x4a x4b x5e x6_1 hormon if random_half==0}{p_end}
{phang}{cmd:. predict xb2, xb}{p_end}
{phang}{cmd:. stcoxcal xb2, val(random_half) times(1(1)6) test}{p_end}

{phang}{cmd:. stcox x1 x4a x4b x5 x6 hormon}{p_end}
{phang}{cmd:. predict xb3, xb}{p_end}
{phang}{cmd:. stcoxcal xb3, times(1(1)6) test}{p_end}


{title:Stored results}

{cmd:stcoxcal} stores the following in {cmd:r()}:

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Scalars}{p_end}
{synopt:{cmd:r(gamma1)}}estimate of gamma1 with gamma0 estimate{p_end}
{synopt:{cmd:r(gamma1_se)}}standard error of gamma1{p_end}
{synopt:{cmd:r(P0)}}p-value for test 1, of gamma0 = 0 given gamma1 = 1{p_end}
{synopt:{cmd:r(P1)}}p-value for test 2, of gamma1 = 1 with gamma0 estimated{p_end}
{synopt:{cmd:r(P01)}}p-value for test 3, joint test of (gamma0, gamma1) = (0, 1){p_end}
{synopt:{cmd:r(Pint)}}p-value for test 4, of interaction of gamma1 with time{p_end}
{p2colreset}{...}

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Macros}{p_end}
{synopt:{cmd:r(fp_pwrs)}}powers of {cmd:_t} in a second-degree fractional polynomial model for baseline log cumulative hazard{p_end}
{p2colreset}{...}


{title:References}

{phang}
Royston, P.  Forthcoming.  Tools for checking calibration of a Cox model in
external validation: Approach based on individual event probabilities.
{it:Stata Journal}.

{phang}
Royston, P., and D. G. Altman.  2013.  External validation of a Cox prognostic
model: Principles and methods. {it:BMC Medical Research Methodology} 13: 33.


{title:Author}

{phang}Patrick Royston{p_end}
{phang}MRC Clinical Trials Unit at UCL{p_end}
{phang}London, UK{p_end}
{phang}j.royston@ucl.ac.uk{p_end}


{marker also_see}{...}
{title:Also see}

{p 4 14 2}Article:  {it:Stata Journal}, volume 14, number 4: {browse "http://www.stata-journal.com/article.html?article=st0357":st0357}

{p 7 14 2}
Manual:  {manhelp fracpoly R}{p_end}

{p 7 14 2}
Online:  {helpb running}, {helpb stpsurv} (if installed){p_end}
