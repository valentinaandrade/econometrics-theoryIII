*Rats program used for:
* Hilde C. Bj?rnland (2009): ?Monetary policy and exchange rate overshooting: Dornbusch was right after all?.
* Journal of International Economics, 79(1), 2009, 64-77.
*
*** 0. LOADING PROCEDURES ***
clear all
source "C:\Program Files\Estima\WinRATS 6.2\shortandlongrun.src"
source "C:\Program Files\Estima\WinRATS 6.2\varirf.src"
source "C:\Program files\Estima\WinRATS 6.2\hpfilter.src"
*
CALENDAR 1978 1 4
ALLOCATE 22 2004:4
*open data C:\Documents and Settings\hilde\My Documents\mcvar\data\VARAUS.XLS
data(format=xls,org=obs) / r GDP CPISA FED RERinv Aldp rtwi oilp du00Q3 du84Q1
set lgdp / = log(GDP(t))
set lop / = log(oilp(t))
set lcpi / = LOG(CPISA(t))
set rdom / = r(t)
set rfed / = FED(t)
set lrexc / = LOG(RERinv(t))
set rfor / = rtwi(t)
*
set rdiff / = rdom(t)-rfor(t)
diff lop / dlop
diff lgdp / dlgdp
diff lcpi / dlcpi
diff lrexc / dlrexc
diff rdom / drdom
diff rfor / drfor
set anninfl / = Aldp
set pi = (cpisa/cpisa{4} - 1)*1
diff pi / dpi
*
*
smpl 1978:1 2004:4
set trend = T
set trendsq = T^2
****@HPFILTER  series  start  end  growth component
@HPFILTER(lambda=1600)  lgdp  1978:1  2004:4  HP_lgdp
set hp_gap / = lgdp(t) - HP_lgdp(t)
*
*
*
graph 1
# dlrexc
*
linreg lgdp / gap
# constant trend trendsq
*
open plot C:\Documents and Settings\hilde\My Documents\NB\graph.gsp
graph 2
# gap
# hp_gap
*
graph 2
# rdom
# rfor
 *
** Dummies
*set du93q1 = t == 1993:1
*
*
system(model=excpuzzle) 1 to 5
VAR rfor lgdp pi dlrexc rdom
lags 1 to 3
det constant du00Q3 du84Q1 trend
*det constant trend
end(system)
scratch 5 / resids
estimate(sigma,outsigma=V,noftests) 1983:1 2004:4 resids+1
*print 1983:1 1987:4 resids+1 resids+2  resids+3
*
*
*
dec rect lr(5,5) sr(5,5)
input lr
. . . . .
. . . . .
. . . . .
. . . . 0
. . . . .
input sr
. 0 0 0 0
. . 0 0 0
. . . 0 0
. . . . .
. . . . .
*
@ShortAndLong(lr=lr,sr=sr,masum=inv(%varlagsums)) %sigma f
*
impulse(model=excpuzzle,decomp=f,steps=24)
*
errors(model=excpuzzle,decomp=f,steps=24,impulses)
*
*
@VARIRF(model=excpuzzle,steps=20,DECOMP=f,byshock,errors)
*
* Options:
*    MODEL=VAR model whose IRF's are to be calculated
*    DECOMP=Factor of the covariance matrix [Choleski of %SIGMA]
*    STEPS=number of forecast steps[48]
*    BYSHOCK/[NOBYSHOCK] Produce a separate graph for each shock with responses to all
*    BYVARIABLE/[NOBYVARIABLE] Produce a separate graph for each variable showing its responses
*     to all shocks
*    ERRORS/[NOERRORS] Produce a forecast error variance decomposition
*    VLABELS=VECTOR[STRINGS] of labels for variables [labels of dependent variables]
*    LABELS=VECTOR[STRINGS] of labels for shocks [same as VLABELS]




















