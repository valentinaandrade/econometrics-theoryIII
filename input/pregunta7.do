********************************************************************************
* Tarea 1. Teoría Econometrica III
* Autor: Camilo Pérez N.
********************************************************************************

clear 
clear all
set more off
global path "C:\Users\Camilo\Documents\PUC\2 Semestre\Econometria III\Tareas\1"
global pathdat "$path\datos"
global pathgr "$path\graficos"
global patht "$path\tablas"

********************************************************************************
* Pregunta 7: empírica
********************************************************************************
import excel "$pathdat\datos.xlsx", sheet("data") cellrange(A1:I400) firstrow

rename Periodo date
label var inf_sa "Inflación total"
label var infsv_sa "Inflación SV"
label var tcn_sa "TCN"
label var ipsa_sa "IPSA"

g fecha=mofd(date)
drop if fecha==.
tsset fecha
format fecha %tm

*GRAFICOS
tsline inf_sa, graphr(c(white)) ytitle(IPC) xtitle(Fecha) tline(2001m1) tline(2020m3) lp(line dash)
graph export "$pathgr\1_serie_inf_sa.png", replace

tsline tcn_sa, graphr(c(white)) ytitle(TCN) xtitle(Fecha) tline(2001m1) tline(2020m3) lp(line dash)
graph export "$pathgr\2_serie_tcn_sa.png", replace

tsline ipsa_sa, graphr(c(white)) ytitle(IPSA) xtitle(Fecha) tline(2001m1) tline(2020m3) lp(line dash)
graph export "$pathgr\3_serie_ipsa_sa.png", replace

********************************************************************************
** INFLACION TOTAL (inf_sa): Serie SA en var m/m 
********************************************************************************

* Gráfico de la serie
tsline inf_sa if fecha>=tm(2001m1) & fecha<=tm(2020m3), graphr(c(white)) ytitle(IPC) xtitle(Fecha) lp(line dash)
graph export "$pathgr\IPC\1_serie_inf_sa.png", replace

* Autocorrelaciones para determinar MA
ac inf_sa if fecha>=tm(2001m1) & fecha<=tm(2020m3), graphr(c(white)) ytitle(Autocorrelaciones) xtitle(Tiempo)
graph save "$pathgr\IPC\2_serie_inf_sa_ac.gph", replace

* Autocorrelaciones parciales para determinar AR
pac inf_sa if fecha>=tm(2001m1) & fecha<=tm(2020m3), graphr(c(white)) ytitle(Autocorrelaciones Parciales) xtitle(Tiempo)
graph save "$pathgr\IPC\3_serie_inf_sa_pac.gph", replace

graph combine "$pathgr\IPC\3_serie_inf_sa_pac.gph" "$pathgr\IPC\2_serie_inf_sa_ac.gph",  graphr(c(white)) c(1)
graph export "$pathgr\IPC\4_pac_ac_serie_inf_sa.png", replace

* Test de raiz unitaria 
dfuller inf_sa if fecha>=tm(2001m1) & fecha<=tm(2020m3), dr reg lags(5)
dfuller inf_sa if fecha>=tm(2001m1) & fecha<=tm(2020m3), reg lags(5)

* Modelos
local m1_ar=1
local m1_ma=1

local m2_ar=2
local m2_ma=1

local m3_ar=2
local m3_ma=2

local m4_ar=3
local m4_ma=2

local m5_ar=3
local m5_ma=3

local m6_ar=4
local m6_ma=3

local m7_ar=4
local m7_ma=4

local m8_ar=5
local m8_ma=4

local m9_ar=5
local m9_ma=5

local m10_ar=5
local m10_ma=6

putexcel set "$patht\IPC.xlsx", replace

forv k=1/10 {
	di "ARMA(`m`k'_ar',`m`k'_ma')"
	putexcel set "$patht\IPC.xlsx", modify sheet("serie_inf_sa")
	local list_p "B C D E F G H I J K"
	local c=word("`list_p'",`k')
	putexcel `c'2="ARMA(`m`k'_ar',`m`k'_ma')"
	putexcel A3="Test F"
	putexcel A4="AIC"
	putexcel A5="BIC"
	putexcel A6="White Noise test"
	putexcel A7="Breusch-Godfrey"
	** Estimación
	arima inf_sa if fecha>=tm(2001m1) & fecha<=tm(2020m3), arima(`m`k'_ar',0,`m`k'_ma')
	** Test F
	test
	putexcel `c'3=`r(p)', nformat(#.00000)
	** Criterios de información
	estat ic
	mat a=r(S)
	putexcel `c'4=matrix(a[1,5])
	putexcel `c'5=matrix(a[1,6])
	** Grafica de raices
	estat aroots, graphr(c(white))
	graph export "$pathgr\IPC\serie_inf_sa_roots_m`k'.png", replace
	** Generación de residuos
	predict r_m_`k', r
	** Test de ruido blanco
	wntestq r_m_`k' if fecha>=tm(2001m1) & fecha<=tm(2020m3)
	putexcel `c'6=`r(p)', nformat(#.00000)
	** Estimación auxiliar para test Breusch-Godfrey
	reg r_m_`k' l(1/`m`k'_ar').inf_sa if fecha>=tm(2001m1) & fecha<=tm(2020m3)
	** Test BG con 10 rezagos
	estat bgodfrey , lags(1/10)
	mat z=r(p)
	putexcel `c'8=mat(z')	, nformat(#.00000)
	** Autocorrelaciones y autocorrelaciones parciales de los residuos de la estimación arma
	ac r_m_`k' if fecha>=tm(2001m1) & fecha<=tm(2020m3), graphr(c(white)) ytitle(Autocorrelaciones) xtitle(Tiempo)
	graph save "$pathgr\IPC\serie_inf_sa_m`k'_ac_res.gph", replace
	pac r_m_`k' if fecha>=tm(2001m1) & fecha<=tm(2020m3), graphr(c(white)) ytitle(Autocorrelaciones Parciales) xtitle(Tiempo)
	graph save "$pathgr\IPC\serie_inf_sa_m`k'_pac_res.gph", replace
	** Combinación de gráficos
	graph combine "$pathgr\IPC\serie_inf_sa_m`k'_pac_res.gph" "$pathgr\IPC\serie_inf_sa_m`k'_ac_res.gph", graphr(c(white)) c(1)
	graph export "$pathgr\IPC\serie_inf_sa_ac_pac_resid_m_`k'.png", replace
}

drop r_m_*

********************************************************************************
** TIPO DE CAMBIO NOMINAL (tcn_sa): Serie SA en var m/m 
********************************************************************************

* Gráfico de la serie
tsline tcn_sa, graphr(c(white)) ytitle(TCN) xtitle(Fecha) lp(line dash)
graph export "$pathgr\TCN\1_serie_tcn_sa.png", replace

* Autocorrelaciones para determinar MA
ac tcn_sa, graphr(c(white)) ytitle(Autocorrelaciones) xtitle(Tiempo)
graph save "$pathgr\TCN\2_serie_tcn_sa_ac.gph", replace

* Autocorrelaciones parciales para determinar AR
pac tcn_sa, graphr(c(white)) ytitle(Autocorrelaciones Parciales) xtitle(Tiempo)
graph save "$pathgr\TCN\3_serie_tcn_sa_pac.gph", replace

graph combine "$pathgr\TCN\3_serie_tcn_sa_pac.gph" "$pathgr\TCN\2_serie_tcn_sa_ac.gph",  graphr(c(white)) c(1)
graph export "$pathgr\TCN\4_pac_ac_serie_tcn_sa.png", replace

* Test de raiz unitaria 
dfuller tcn_sa, dr reg lags(1)
dfuller tcn_sa, reg lags(1)

* Modelos
local m1_ar=1
local m1_ma=1

local m2_ar=2
local m2_ma=1

local m3_ar=2
local m3_ma=2

local m4_ar=3
local m4_ma=2

local m5_ar=3
local m5_ma=3

local m6_ar=4
local m6_ma=3

local m7_ar=4
local m7_ma=4

local m8_ar=5
local m8_ma=4

local m9_ar=5
local m9_ma=5

local m10_ar=5
local m10_ma=6

local m11_ar=6
local m11_ma=6

local m12_ar=6
local m12_ma=7

putexcel set "$patht\TCN.xlsx", replace

forv k=1/12 {
	di "ARMA(`m`k'_ar',`m`k'_ma')"
	putexcel set "$patht\TCN.xlsx", modify sheet("serie_tcn_sa")
	local list_p "B C D E F G H I J K L M"
	local c=word("`list_p'",`k')
	putexcel `c'2="ARMA(`m`k'_ar',`m`k'_ma')"
	putexcel A3="Test F"
	putexcel A4="AIC"
	putexcel A5="BIC"
	putexcel A6="White Noise test"
	putexcel A7="Breusch-Godfrey"
	** Estimación
	arima tcn_sa, arima(`m`k'_ar',0,`m`k'_ma')
	** Test F
	test
	putexcel `c'3=`r(p)', nformat(#.00000)
	** Criterios de información
	estat ic
	mat a=r(S)
	putexcel `c'4=matrix(a[1,5])
	putexcel `c'5=matrix(a[1,6])
	** Grafica de raices
	estat aroots, graphr(c(white))
	graph export "$pathgr\TCN\serie_tcn_sa_roots_m`k'.png", replace
	** Generación de residuos
	predict r_m_`k', r
	** Test de ruido blanco
	wntestq r_m_`k'
	putexcel `c'6=`r(p)', nformat(#.00000)
	** Estimación auxiliar para test Breusch-Godfrey
	reg r_m_`k' l(1/`m`k'_ar').tcn_sa
	** Test BG con 10 rezagos
	estat bgodfrey , lags(1/10)
	mat z=r(p)
	putexcel `c'8=mat(z')	, nformat(#.00000)
	** Autocorrelaciones y autocorrelaciones parciales de los residuos de la estimación arma
	ac r_m_`k', graphr(c(white)) ytitle(Autocorrelaciones) xtitle(Tiempo)
	graph save "$pathgr\TCN\serie_tcn_sa_m`k'_ac_res.gph", replace
	pac r_m_`k', graphr(c(white)) ytitle(Autocorrelaciones Parciales) xtitle(Tiempo)
	graph save "$pathgr\TCN\serie_tcn_sa_m`k'_pac_res.gph", replace
	** Combinación de gráficos
	graph combine "$pathgr\TCN\serie_tcn_sa_m`k'_pac_res.gph" "$pathgr\TCN\serie_tcn_sa_m`k'_ac_res.gph", graphr(c(white)) c(1)
	graph export "$pathgr\TCN\serie_tcn_sa_ac_pac_resid_m_`k'.png", replace
}

drop r_m_*


********************************************************************************
** ÍNDICE DE PRECIOS SELECTIVO DE ACCIONES (ipsa_sa): Serie SA en var m/m 
********************************************************************************

* Gráfico de la serie
tsline ipsa_sa, graphr(c(white)) ytitle(TCN) xtitle(Fecha) lp(line dash)
graph export "$pathgr\IPSA\1_serie_ipsa_sa.png", replace

* Autocorrelaciones para determinar MA
ac ipsa_sa, graphr(c(white)) ytitle(Autocorrelaciones) xtitle(Tiempo)
graph save "$pathgr\IPSA\2_serie_ipsa_sa_ac.gph", replace

* Autocorrelaciones parciales para determinar AR
pac ipsa_sa, graphr(c(white)) ytitle(Autocorrelaciones Parciales) xtitle(Tiempo)
graph save "$pathgr\IPSA\3_serie_ipsa_sa_pac.gph", replace

graph combine "$pathgr\IPSA\3_serie_ipsa_sa_pac.gph" "$pathgr\IPSA\2_serie_ipsa_sa_ac.gph",  graphr(c(white)) c(1)
graph export "$pathgr\IPSA\4_pac_ac_serie_ipsa_sa.png", replace

* Test de raiz unitaria 
dfuller ipsa_sa, dr reg lags(7)
dfuller ipsa_sa, reg lags(7)

* Modelos
local m1_ar=1
local m1_ma=1

local m2_ar=2
local m2_ma=1

local m3_ar=2
local m3_ma=2

local m4_ar=3
local m4_ma=2

local m5_ar=3
local m5_ma=3

local m6_ar=4
local m6_ma=3

local m7_ar=4
local m7_ma=4

local m8_ar=5
local m8_ma=4

local m9_ar=5
local m9_ma=5

local m10_ar=5
local m10_ma=6

local m11_ar=6
local m11_ma=6

local m12_ar=6
local m12_ma=7

putexcel set "$patht\IPSA.xlsx", replace

forv k=1/12 {
	di "ARMA(`m`k'_ar',`m`k'_ma')"
	putexcel set "$patht\IPSA.xlsx", modify sheet("serie_ipsa_sa")
	local list_p "B C D E F G H I J K L M"
	local c=word("`list_p'",`k')
	putexcel `c'2="ARMA(`m`k'_ar',`m`k'_ma')"
	putexcel A3="Test F"
	putexcel A4="AIC"
	putexcel A5="BIC"
	putexcel A6="White Noise test"
	putexcel A7="Breusch-Godfrey"
	** Estimación
	arima ipsa_sa, arima(`m`k'_ar',0,`m`k'_ma')
	** Test F
	test
	putexcel `c'3=`r(p)', nformat(#.00000)
	** Criterios de información
	estat ic
	mat a=r(S)
	putexcel `c'4=matrix(a[1,5])
	putexcel `c'5=matrix(a[1,6])
	** Grafica de raices
	estat aroots, graphr(c(white))
	graph export "$pathgr\IPSA\serie_ipsa_sa_roots_m`k'.png", replace
	** Generación de residuos
	predict r_m_`k', r
	** Test de ruido blanco
	wntestq r_m_`k'
	putexcel `c'6=`r(p)', nformat(#.00000)
	** Estimación auxiliar para test Breusch-Godfrey
	reg r_m_`k' l(1/`m`k'_ar').ipsa_sa
	** Test BG con 10 rezagos
	estat bgodfrey , lags(1/10)
	mat z=r(p)
	putexcel `c'8=mat(z'), nformat(#.00000)
	** Autocorrelaciones y autocorrelaciones parciales de los residuos de la estimación arma
	ac r_m_`k', graphr(c(white)) ytitle(Autocorrelaciones) xtitle(Tiempo)
	graph save "$pathgr\IPSA\serie_ipsa_sa_m`k'_ac_res.gph", replace
	pac r_m_`k', graphr(c(white)) ytitle(Autocorrelaciones Parciales) xtitle(Tiempo)
	graph save "$pathgr\IPSA\serie_ipsa_sa_m`k'_pac_res.gph", replace
	** Combinación de gráficos
	graph combine "$pathgr\IPSA\serie_ipsa_sa_m`k'_pac_res.gph" "$pathgr\IPSA\serie_ipsa_sa_m`k'_ac_res.gph", graphr(c(white)) c(1)
	graph export "$pathgr\IPSA\serie_ipsa_sa_ac_pac_resid_m_`k'.png", replace
}

drop r_m_*

********************************************************************************
** PREDICCION DE MOELOS SELECCIONADOS
********************************************************************************

*IPC
arima inf_sa if fecha>=tm(2001m1) & fecha<=tm(2020m3), arima(5,0,4)
predict pred_inf if fecha>=tm(2001m1) & fecha<=tm(2020m3)
tsline inf_sa pred_inf if fecha>=tm(2001m1) & fecha<=tm(2020m3), graphr(c(white)) ytitle(IPC) xtitle(Fecha) lp(line dash)
graph export "$pathgr\4_fit_ipc_sa.png", replace

*IPSA
arima ipsa_sa, arima(6,0,7)
predict pred_ipsa
tsline ipsa_sa pred_ipsa, graphr(c(white)) ytitle(IPSA) xtitle(Fecha) lp(line dash)
graph export "$pathgr\5_fit_ipsa_sa.png", replace

*TCN
arima tcn_sa, arima(1,0,1)
predict pred_tcn
tsline tcn_sa pred_tcn, graphr(c(white)) ytitle(TCN) xtitle(Fecha) lp(line dash)
graph export "$pathgr\6_fit_tcn_sa.png", replace
