*************************************************************************************************************************************************************************** ;
*	BRIDGING EFFICACY ANALYSIS WITH PRESCREENED SAMPLES                                                                                                                     ;
*   Version 1.1, Uploaded 02/28/2023 												                                                                                        ;
*																					               										                                    ;
*	Code written by Wei Wang (2023)													               										                                    ;
* 																					                										                                ;
*	Reference:																		               										                                    ;
* 	Wang W, Nair RJ. Statistical Consideration in Bridging Study of Personalized Medicine with Biomarker-based Prescreened Samples                                          ;
* 	Submitted. 				                    										                                                                                    ;
* 																					                                                                                        ;
*************************************************************************************************************************************************************************** ;
* 																					                										                                ;
*	PROGRAM DESCRIPTION:															                                                                                        ;
*   ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- ;
* 	This program provided statistical methods to address prescreening bias and estimate the drug efficacy in CDx intended use population in the bridging efficacy analysis  ;
* 	using CTA-drug pivotal clinical trial outcome data and three-way concordance results between CDx, CTA and LDT. This SAS macro can be applied for the binary drug        ;
* 	efficacy of the treatment in cancer clinical trial, e.g. the objective response rate (ORR) defined as the proportion of Complete response (CR) and partial response     ;
* 	(PR), or proportion of death within pre-specified time window after cancer treatment                                                                                    ;
* 																					                                                                                        ;
*	MACRO VARIABLES:																                                                                                        ;
*   ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- ;
*	Three Concordance Table Input:                                                                                                               	                        ;
*	N111: 		    number of subjects with LDT = 1 (positive), CTA = 1 (positive) and CDx = 1 (positive)                                          	                        ;
*																					                                                                                        ;
*	N101: 		    number of subjects with LDT = 1 (positive), CTA = 0 (negative) and CDx = 1 (positive)                                          	                        ;
*																					                                                                                        ;
*	N110: 		    number of subjects with LDT = 1 (positive), CTA = 1 (positive) and CDx = 0 (negative)                                          	                        ;
*																					                                                                                        ;
*	N100: 		    number of subjects with LDT = 1 (positive), CTA = 0 (negative) and CDx = 0 (negative)                                          	                        ;
*																					                                                                                        ;
*	N011: 		    number of subjects with LDT = 0 (positive), CTA = 1 (positive) and CDx = 1 (positive)                                          	                        ;
*																					                                                                                        ;
*	N001: 		    number of subjects with LDT = 0 (positive), CTA = 0 (negative) and CDx = 1 (positive)                                          	                        ;
*																					                                                                                        ;
*	N010: 		    number of subjects with LDT = 0 (positive), CTA = 1 (positive) and CDx = 0 (negative)                                          	                        ;
*																					                                                                                        ;
*	N000: 		    number of subjects with LDT = 0 (positive), CTA = 0 (negative) and CDx = 0 (negative)                                          	                        ;
*																					                                                                                        ;
*	P_D1: 		    prevalence of LDT = 1 (positive) in the intended use population                                                             	                        ;
*																					                                                                                        ;
*	N_TREAT: 		number of triple positve subjects (LDT = 1, CTA = 1 and CDx = 1) in the drug clinical trial                                    	                        ;
*																					                                                                                        ;
*	N_RES: 		    number of subjects with positive outcome (e.g. ORR, or proportion of death within pre-specified time window after cancer treatment) in the triple       ;
*					positives of the drug clinical trial                                                                                                                    ;
*																					                                                                                        ;
*	SEEDBOOT:		the random-number seed to generate bootstrap samples 	                                                                                                ;
*																					                                                                                        ;
*	PRINTALL:		Outpout data set print out or not	                                                                                                                    ;
*					F: output data set is not printed out								                                                                                    ;
*					T: output data set is printed out, default												                                                                ;
*																					                                                                                        ;
*	OUT:	        Output dataset                                                                                                                                          ;
*																					                                                                                        ;
*	OUTPUT DATA VARIABLES:  														                                                                                        ;
*   ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- ;
*																					                                                                                        ;
*	PROP: c value, %, indicating that response rate for the positive discordant patients, i.e., (CTA- or LDT-)nCDx+, d((LDT- or CTA-)nCDx+),                                ;
*	      is assumed as c% of the observed response rate in the positive concordant patients, i.e., LDT+nCTA+nCDx+, d(LDT+nCTA+nCDx+)		                                ;
*																					                                                                                        ;
*	RES: the response rate in the CDx+ intended use population, dCDx+ under different c values, %                                                                           ;
*																					                                                                                        ;
*	RES_LB: the lower bound of the 95% confidence intervals for the response rate in the CDx+ intended use population, dCDx+ under different c values                       ;
*																					                                                                                        ;
*	RES_UB: the upper bound of the 95% confidence intervals for the response rate in the CDx+ intended use population, dCDx+ under different c values                       ;
*																					                                                                                        ;
*	EXAMPLE CODE:																	                                                                                        ;
*																					                                                                                        ;
*   %include 'prescreenbridge.sas'												                                                                                            ;
*																					                                                                                        ;
*   %PRESCREENBRIDGE   (N111             = 126,                                                                                                                             ;
*                       N101             = 2,                                                                                                                               ;
*                       N110             = 1,                                                                                                                               ;
*                       N100             = 3,                                                                                                                               ;
*                       N011             = 4,                                                                                                                               ;
*                       N001             = 2,                                                                                                                               ;
*                       N010             = 0,                                                                                                                               ;
*                       N000             = 84,                                                                                                                               ;
*                       P_D1             = 0.3824,                                                                                                                          ;
*                       N_TREAT          = 120,                                                                                                                             ;
*                       N_RES            = 42,                                                                                                                              ;
*                       SEEDBOOT         = 1234667,                                                                                                                         ;
*                       PRINTALL         = T,                                                                                                                               ;
*                       OUT              = out1                                                                                                                             ;
*                       );                                                                                                                                                  ;
* 																					                                                                                        ;
*************************************************************************************************************************************************************************** ;
* 																					                                                                                        ;
*************************************************************************************************************************************************************************** ;
***********************************************************************************Final Macro******************************************************************************;
*************************************************************************************************************************************************************************** ;

dm  'log;clear;out;clear;';

OPTIONS ls=130  ps=57 NOCENTER DATE PAGENO=1; 

*******************************************************************************************************
                                             Main Macro;
*******************************************************************************************************;

%macro PRESCREENBRIDGE   (N111             = , 
                          N101             = ,
                          N110             = ,
                          N100             = ,
                          N011             = ,
                          N001             = ,
                          N010             = ,
                          N000             = ,
                          P_D1             = ,
                          N_TREAT          = , 
                          N_RES            = , 
                          seedboot         = 1237893,
                          Printall         = T,
                          Out              = 
                          );

data multinomial1;
x1 = &N111;
x2 = &N101;
x3 = &N110;
x4 = &N100;
x5 = &N011;
x6 = &N001;
x7 = &N010;
x8 = &N000;
repl = 1;
run;

data aa1;
set multinomial1;
id = 0;
if repl le 1000;
do i = 1 to x1;
LDT = 1;
CTA = 1;
CDX = 1;
id = id + 1;
output;
end;
do j = 1 to x2;
LDT = 1;
CTA = 0;
CDX = 1;
id = id + 1;
output;
end;
do k = 1 to x3;
LDT = 1;
CTA = 1;
CDX = 0;
id = id + 1;
output;
end;
do l = 1 to x4;
LDT = 1;
CTA = 0;
CDX = 0;
id = id + 1;
output;
end;
do m = 1 to x5;
LDT = 0;
CTA = 1;
CDX = 1;
id = id + 1;
output;
end;
do n = 1 to x6;
LDT = 0;
CTA = 0;
CDX = 1;
id = id + 1;
output;
end;
do o = 1 to x7;
LDT = 0;
CTA = 1;
CDX = 0;
id = id + 1;
output;
end;
do p = 1 to x8;
LDT = 0;
CTA = 0;
CDX = 0;
id = id + 1;
output;
end;
drop i j k l m n o p;
run;

data aa2;
set aa1;
tem = 1;
run;

%let seed1 = %eval(&seedboot + 1234566);

proc surveyselect data= aa2 out=boots
     method = urs
	 samprate = 1 outhits rep = 1000 seed = &seed1;
	 strata repl;
run;

data aa3;
retain repl replicate id;
set aa2 (in = a) boots ;
if a then replicate = 0;
*if replicate = 0;
drop NumberHits ExpectedHits SamplingWeight;
run;

proc sort data = aa3;
by repl replicate tem;
run;

data aa4;
set aa3;
by repl replicate tem;
retain n111 n110 n101 n100 n011 n010 n001 n000;
if first.tem then do;
n111 = 0; 
n101 = 0;
n110 = 0;
n100 = 0;
n011 = 0;
n001 = 0;
n010 = 0;
n000 = 0;
end;
if ldt = 1 and cta = 1 and cdx = 1 then n111 = n111+1;
if ldt = 1 and cta = 0 and cdx = 1 then n101 = n101+1;
if ldt = 1 and cta = 1 and cdx = 0 then n110 = n110+1;
if ldt = 1 and cta = 0 and cdx = 0 then n100 = n100+1;
if ldt = 0 and cta = 1 and cdx = 1 then n011 = n011+1;
if ldt = 0 and cta = 0 and cdx = 1 then n001 = n001+1;
if ldt = 0 and cta = 1 and cdx = 0 then n010 = n010+1;
if ldt = 0 and cta = 0 and cdx = 0 then n000 = n000+1;
if last.replicate;
drop id tem ldt cta cdx x1 -- x8;
run;

data aa5;
set aa4;
p_d1 = &p_d1;
p = n111/(n111+n110+n101+n100) * p_d1/((n111+n101)/(n111+n110+n101+n100) * p_d1 + (n011+n001)/(n011+n010+n001+n000)*(1 - p_d1));
p_no = n111/(n111+n110) * p_d1/(n111/(n111+n110) * p_d1 + (n001+n101)/(n101+n001+n100+n000)*(1 - p_d1));
run;

data simula1;
x1 = &n_res;
x2 = &n_treat - &n_res;
repl = 1;
run;

data bb1;
set simula1;
id = 0;
if repl le 1000;
do i = 1 to x1;
res = 1;
id = id + 1;
output;
end;
do j = 1 to x2;
res = 0;
id = id + 1;
output;
end;
drop i j;
run;

data bb2;
set bb1;
tem = 1;
run;

%let seed2 = %eval(&seedboot + 2234566);

proc surveyselect data= bb2 out=boots1
     method = urs
	 samprate = 1 outhits rep = 1000 seed = &seed2;
	 strata repl;
run;

data bb3;
retain repl replicate id;
set bb2 (in = a) boots1 ;
if a then replicate = 0;
*if replicate = 0;
drop NumberHits ExpectedHits SamplingWeight;
run;

proc sort data = bb3;
by repl replicate tem;
run;

data bb4;
set bb3;
by repl replicate tem;
retain n1 n0;
if first.tem then do;
n1 = 0; 
n0 = 0;
end;
if RES = 1 then n1 = n1+1;
if RES = 0 then n0 = n0+1;
if last.replicate;
drop id tem RES x1 x2;
run;

data bb5;
set bb4;
p_res = n1/(n1 + n0);
run;

data cc1;
merge aa5 bb5;
by repl replicate;
run;

data cc2;
set cc1;
do i = 1 to 101;
res = p_res * p + (i - 1) * 0.01 * p_res * (1 - p); 
res_no = p_res * p_no + (i - 1) * 0.01 * p_res * (1 - p_no); 
output;
end;
run;

proc sort data = cc2;
by i replicate;
run;

proc univariate data = cc2;
  var res res_no;
  output out=cc3 pctlpts=2.5 97.5 pctlpre = res resno pctlname = _lb _ub ;
  by i;
  where replicate ne 0;
run;

data cc4;
set cc2;
if replicate = 0;
run;

data cc5;
retain i res res_lb res_ub res_no resno_lb resno_ub;
merge cc4 cc3;
drop replicate -- p_res repl;
run;

data cc6;
retain prop res res_lb res_ub;
set cc5;
prop = (i - 1);
drop i res_no resno_lb resno_ub;
run;

data &out;
set cc6;
run;

%if &printall=T %then %do;

  Title 'Response Rate in the CDx+ Intended Use Population';

  proc print data = &out noobs label;
    label prop = 'c Value, %';
    label res = 'Response Rate in CDx+ Subjects';
    label res_lb = '95% Lower Limit for Response Rate in CDx+ Subjects';
    label res_ub = '95% Upper Limit for Response Rate in CDx+ Subjects';
	var prop res res_lb res_ub;
  run;

%end;  *-end printall=T option;

proc datasets lib = work;
delete aa1 aa2 aa3 aa4 aa5 bb1 bb2 bb3 bb4 bb5 boots boots1 cc1 cc2 cc3 cc4 cc5 cc6 multinomial1 simula1;
run;

%mend;

%PRESCREENBRIDGE   (N111             = 126, 
                          N101             = 2,
                          N110             = 1,
                          N100             = 3,
                          N011             = 4,
                          N001             = 2,
                          N010             = 0,
                          N000             = 84,
                          P_D1             = 0.3824,
                          N_TREAT          = 120, 
                          N_RES            = 42, 
                          seedboot         = 1237893,
                          Printall         = T,
                          Out              = temp1
                          );
