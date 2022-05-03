* Encoding: UTF-8.

* Make sure the Active Dataset shows your data file (.sav).
* Steps:
* 1. Open EndoS.sps.
*     Select all, and Run Selection (click the triangle green button in the SPSS menu - next to the binocular button..
*     If you do this, you have just invoked the macro in the SPSS memory, which you can call during your SPSS session.
* 2. Highlight the code below, and Run Selection (click the triangle green button in the SPSS menud).

TSLS dv = lwage
/iv = educ 
/ivm = motheduc fatheduc
/olsoutp= 0 
/reshaus= 1 
/overz = 1
/bptest= 0
/wktest=1
/robse= 99.

* Documentation.
* dv = dependent variable.
* iv = endogeneous regressor.
* ivm = list of instrument(s).
* olsoutp,  1= OLS outputs printed, 0 = not print.
* reshaus, 1 = run Hausman's specifcation test; 0 = not run. 
*overz, 1 = run overidentification test, 0 = not run
*bptest, 1 = run Breusch-Pagan test for heteroskedasticity, 0=not run.
*wktest, 1 = run weak instrument's test, 0=not run.
*robse, 99 = no robust std error, 0 = HC0, 1=HC1, 2=HC2, 3=HC3, 4=HC4.
