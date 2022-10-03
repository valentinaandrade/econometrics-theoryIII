Author: Markus Kirchner (MK)
Date:	16 June 2014
e-mail: mkirchner@bcentral.cl

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CONTENTS OF THIS FOLDER

------------------------------------------------------------------------------------------------------------------------------------------------

*** Main file ***

ModeloBASE				:	Demonstrates the use of the functions provided below, including the required settings.
					Estimation and calculation of impulse responses for a structural VAR with exclusion restrictions.

------------------------------------------------------------------------------------------------------------------------------------------------

*** Functions ***

CheckStability			:	Checks stability of the VAR coefficients.
ComputeConfidenceBands		:	Computes bootstrapped confidence bands for impulse responses of the VAR.
ComputeHistoricalDecomposition:	Provides historical decomposition based on the identified (Cholesky) VAR.
ComputeImpulseResponses		: 	Computes point estimates of impulse responses for the identified (Cholesky) VAR.
ComputeShockDecomposition	: 	Provides structural shocks for the identified (Cholesky) VAR.
EstimateVAR				:	Estimates the reduced-form VAR estimated by iterated FGLS or OLS.
GenerateCompanionMatrix		:	Generates VAR(1) companion matrix of the VAR(p) from autoregressive lag matrices.
GenerateLagMatrices		:	Generates matrices of the VAR(p) lag polynominal from autoregressive coefficients in simult. eqs. form.
GenerateSimultEqsMatrices	:	Generates matrices of the simulataneous equations form of the VAR(p) with exclusion restrictions.
LoadData				:	Loads data for the vector autoregression.
loc (borrowed)			: 	Locates a given string in a vector, written by Guenter Coenen, modified by MK.
mprint (borrowed)			:	Prints an (nobs x nvar) matrix in formatted form.
PlotImpulseResponses		:	Plots impulse responses for a VAR.
PlotRegressionResults		:	Plots data against fitted values as well as residuals from a VAR.
PrepareRegressionData		:	Prepares data for the reduced-form VAR that allows for general exclusion restrictions.
PrintTestResults			:	Prints results of various tests conducted on the VAR.
round2 (borrowed)			:	Round to a specified number of decimals.
shadedplot (borrowed)		: 	Produces shaded areas between two lines for impulse response error bands, modified by MK.
TestCoefficientRestrictions	:	Tests coefficient restrictions of the VAR using a likelihood ratio test.
TestLagLength			:	Tests for appropriate lag length of the VAR.
TestTrendSpecification		:	Tests for appropriate specification of the deterministic terms in the VAR.
vec (borrowed)			:	Create a matrix stacking the columns of Y.

------------------------------------------------------------------------------------------------------------------------------------------------

*** Contents of sub-folders ***
	
backup				:	Backup of Matlab code and data, as provided by MK (read only).
data					:	Data base including VAR data, Excel versions 1997-2003 for use in Matlab R2010a.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

INSTRUCTIONS

Note: Detailed instructions are provided in the main file and the documentation of the individual functions.

------------------------------------------------------------------------------------------------------------------------------------------------

1. Different VAR specifications are provided in "\data\base.xls".
2. Run the main file, change specifications and settings there.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CHANGE LOG

21-11-2011: Corrected a bug in LoadData, variable names now load correctly independent of their ordering in DemoEstimation relative to data.xls
07-02-2012: Corrected a mistake in EstimateBlockVAR, total number of parameters for information criteria is now correctly imputed
20-02-2012: Corrected ComputeImpulseResponses, ComputeConfidenceBands, PlotImpulseResponses, normalization now only for the shock of interest

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%