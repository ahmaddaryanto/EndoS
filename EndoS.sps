* Encoding: UTF-8.

DATASET ACTIVATE DataSet3.

* Encoding: UTF-8.
/*Two-Stage Least Square. 
/*Macro created by Ahmad Daryanto*/.
/*July 2018. 

set printback=off mprint=off errors=off .
preserve.
SET WORKSPACE=134576.
DEFINE OLSMACRO (Xols = !charend('/')
/Yols = !charend('/')
/numcol = !charend('/')
/numrow = !charend('/')
/xnames = !charend('/')
/ynames = !charend('/')).
/*OLS.
compute nmvars={!xnames}.
print !ynames/Title="Dependent variable"/format=A9.
compute bols=(inv(sscp(!Xols)))*t(!Xols)*!Yols.
compute e=!Yols-!Xols*bols.
compute e2=e(:,1)&*e(:,1).
compute sser=csum(sscp(e)).
compute mse=(1/(!numrow-!numcol))*sser.
compute vbols=mse*inv(sscp(!Xols)).
compute sbols=sqrt(diag(vbols)).
compute tbols=bols/sbols.
compute dffols=!numrow-!numcol.
compute Fols=tbols&*tbols.
compute pFols=1-fcdf(Fols,1,dffols).
     *--95% CI--.
compute LBols=bols-1.96*sbols.
compute UBols=bols+1.96*sbols.
compute olsout={bols,sbols, tbols,pFols, LBols, UBols}.
*end of calculation.
*===Preparing input ANOVA table of OLS.
*computing mean square regression.
compute meanY=ones*t(csum(!Yols)/!numrow).
compute e_reg=meanY-!Xols*bols.
compute ssreg=csum(sscp(e_reg)).
compute sumsq=T({ssreg,sser}).
compute dfa=T({!numcol-1,!numrow-!numcol}).
compute mse_a=sumsq/dfa.
compute Fval=(ssreg/(!numcol-1))/(sser/(!numrow-!numcol)).
Compute pF_a=1-fcdf(Fval,!numcol-1,!numrow-!numcol).
compute F_a=T({Fval,-999}).
compute pFa=T({pF_a,-999}).
*--computing model summary.
compute sereg=sqrt(mse_a(2)).
Compute Rsq=ssreg/(sser+ssreg).
Compute Rsqadj=1-((1-Rsq)*(!numrow-1)/(!numrow-!numcol-1)).
        compute Rsqs=({Rsq,Rsqadj,sereg}).
        print Rsqs/title="Model summary"/cnames=nmvars1/format=f9.4.
                    compute rlabel={"df(m)","df(res)","F","Sig."}.
                    compute out=({t(dfa),F_a(1),pFa(1)}).
                     print out/title=" "/cnames=rlabel/format=f9.4. 
        *--OLS output.
      
        print /title="Note: Adj = Adjusted R^2, SER = Std. error of regression.".
        print olsout/title ="======="/rnames=nmvars/cnames=cnms/format=F9.3.
!ENDDEFINE.
/*&&&&&&&&&&&&&&&&&&&&.


DEFINE CraggD (Y = !charend('/')
/X = !charend('/')
/Z = !charend('/')).
compute  #Y=!Y.
compute #X=!X.
compute #Z=!Z.
compute Indt=IDENT(nrow(#Z),nrow(#Z)).
compute Z_={#X,#Z}.
compute MZ_=(Indt-Z_*(inv(sscp(Z_)))*t(Z_)).
compute MX=(Indt-#X*(inv(sscp(#X)))*t(#X)).
compute Y_T=MX*#Y.
Compute Z_T=MX*#Z.
Compute PZ_T=Z_T*(inv(sscp(Z_T)))*t(Z_T).
compute numT=nrow(#Z).
compute K1=ncol(#X).
compute K2=ncol(#Z).
compute Sigv=(t(#Y)*MZ_*#Y)/(numT-K2-K1).
compute sqrSigv=CHOL(Sigv).
compute invSigv=inv(sqrSigv).
/*compute fac=(numT-K2-K1)/K2.
compute matG=t(invSigv) *t(Y_T) * PZ_T *Y_T* invSigv/K2.
compute gmin=cmin(Eval(matG)).
print gmin/title="Cragg-Donald F-statistic"/format=f9.4.
!ENDDEFINE.
/*&&&&&&&&&&&&&&&&&&&&.
DEFINE TSLS (iv = !charend('/')
/dv = !charend('/')
/ivm=!charend('/')
/olsoutp=!charend('/')!default(0)
/reshaus=!charend('/')!default(0)
/overz=!charend('/')!default(0)
/bptest=!charend('/')!default(0)
/wktest=!charend('/')!default(0)
/robse=!charend('/')!default(99)).
SET MXLOOPS = 10000001.
SET PRINTBACK = OFF.
MATRIX.
get mat/variables=!dv !iv /names=nms /MISSING=OMIT.
get matz/variables=!ivm/names=nms2/MISSING=OMIT.
compute n=nrow(mat).
compute ones=make(n,1,1).
compute Z={ones,matZ}.
compute Y=mat(:,1).
compute X={ones,mat(:,2:ncol(mat))}.
compute k=ncol(X).
print/ title="Two-Stage Least Squares".
print/title="written by Ahmad Daryanto".
print/title="https://sites.google.com/site/ahmaddaryanto/".
print/title="-----------------------------------------".
compute dvname=t(nms(:,1)).
compute ivnames=t(nms(:,2:k)).
    DO if (!olsoutp=1).
    print/title="OLS regression".
    print n/title="Sample size".
    print dvname/title="Dependent variable (DV)"/ format=A10.
    print t(ivnames)/title="Independent variable (IVs)"
            / format=A10.
    End if.
*===============================================.
*            OLS Regression of Without Instruments .
*===============================================.
*dv in original metrix.
compute b=(inv(sscp(X)))*t(X)*Y.
/*---.
/*notes: sser=sum squares of error.ssreg=sum of squares of regression.
*===computing residuals, standard error of b, t value and p-value of OLS ==.
compute e=Y-X*b.
compute e2=e(:,1)&*e(:,1).
compute sser=csum(sscp(e)).
compute mse=(1/(n-k))*sser.
compute vb=mse*inv(sscp(X)).
compute sb=sqrt(diag(vb)).
compute tb=b/sb.
compute dff=n-k.
compute F=tb&*tb.
compute pF=1-fcdf(F,1,dff).
compute pF=1-fcdf(F,1,dff).
     *--95% CI--.
compute LB=b-1.96*sb.
compute UB=b+1.96*sb.
*for output with robust std errors. 
* No correction (default).
do if (!robse=99).
compute  vbh=vb.
   end if.
*HC0.  
do if (!robse=0).
    compute  vbh=inv(sscp(X))*t(X)*mdiag(e2)*X*inv(sscp(X)).
end if.
* HC1.
do if (!robse=1).
compute  vbh=inv(sscp(X))*t(X)*mdiag(e2)*X*inv(sscp(X)).
    compute vbh=vbh*N/(N-k).
end if.
*HC2.
do if (!robse=2).
compute hat=X*inv(sscp(X))* t(X).
compute dhat=e2&/(ones-diag(hat)).
compute vbh=inv(sscp(X))*t(X)*mdiag(dhat)*X*inv(sscp(X)).
end if.
*HC3.
do if (!robse=3).
compute hat=X*inv(sscp(X))* t(X).
compute hat2=(ones-diag(hat))&*(ones-diag(hat)).
compute dhat=e2&/hat2.
compute vbh=inv(sscp(X))*t(X)*mdiag(dhat)*X*inv(sscp(X)).
end if.
*HC4.
do if (!robse=4).
compute hat=X*inv(sscp(X))* t(X).
compute fours=make(n,1,4).
compute mh={fours,n*diag(hat)/k}.
compute dummy=rmin(mh).
compute hat2=(ones-diag(hat))&**dummy.
compute dhat=e2&/hat2.
compute vbh=inv(sscp(X))*t(X)*mdiag(dhat)*X*inv(sscp(X)).
end if.
    compute sbh=sqrt(diag(vbh)).
    compute tbh=b/sbh.
    compute dff=n-k.
    compute Fh=tbh&*tbh.
    compute pFh=1-fcdf(Fh,1,dff).
    compute pF=1-fcdf(Fh,1,dff).
     *--95% CI--.
   compute LBh=b-1.96*sb.
   compute UBh=b+1.96*sb.
   compute olsouth={b,sbh, tbh,pFh, LBh, UBh}.
*end of calculation.
*===Preparing input ANOVA table of OLS.
*computing mean square regression.
compute meanY=ones*t(csum(Y)/n).
compute e_reg=meanY-X*b.
compute ssreg=csum(sscp(e_reg)).
compute sumsq=T({ssreg,sser}).
compute dfa=T({k-1,n-k}).
compute mse_a=sumsq/dfa.
compute Fval=(ssreg/(k-1))/(sser/(n-k)).
Compute pF_a=1-fcdf(Fval,k-1,n-k).
compute F_a=T({Fval,-999}).
compute pFa=T({pF_a,-999}).
*--computing model summary.
compute sereg=sqrt(mse_a(2)).
Compute Rsqo=ssreg/(sser+ssreg).
Compute Rsqadj=1-((1-Rsqo)*(n-1)/(n-k-1)).
        compute nmv={"SER"}.
        compute nmvars1 = {"R^2"; "Adj";nmv}.
        compute Rsqs=({Rsqo,Rsqadj,sereg}).
        compute nmvars = t(nms(1,2:ncol(mat))).
        compute nmvars = {"constant"; nmvars}.
        compute cnms={"b","se", "t", "sig", "95%LB", "95%UB"}.
DO IF (!olsoutp=1).
        print Rsqs/title="Model summary"/cnames=nmvars1/format=f9.4.
                    compute rlabel={"df(m)","df(res)","F","Sig."}.
                    compute out=({t(dfa),F_a(1),pFa(1)}).
                    
        print /title="Note: Adj = Adjusted R^2, SER = Std. error of regression.".
       /* print {sumsq,dfa, mse_a,F_a,pFa} /space=3
          /*title '------- ANOVA TABLE  --------'.
          /*clabel "SS" "df" "MS" "F" "Sig".
          /*rlabel "Model" "Residual".
          /*format f10.3 .
        *--OLS output.
print olsouth/title ="OLS outputs"/rnames=nmvars/cnames=cnms/format=F9.3.
*Correcting for standard errors--robust standard errors.
do if (!robse=99).
print/title="* Note: standard errors are assumed to be homoskedastic--no adjustments.".
end if.
do if (!robse=0).
print/title="* Note: standard errors are adjusted using HC0 variant (Eicker-Huber–White standard "+
    "errors).".
print/title="* The adjustments are not recommended for sample sizes < 250 (Long and Ervin, 2000).".
end if.
do if (!robse=1).
print/title="* Note: standard errors are adjusted using HC1 variant.".
end if.
do if (!robse=2).
print/title="* Note: standard errors are adjusted using HC2 variant.".
end if.
do if (!robse=3).
print/title="* Note: standard errors are adjusted using HC3 variant.".
end if.
do if (!robse=4).
print/title="* Note: standard errors are adjusted using HC4 variant.".
end if.
        
END IF.
/*===============================================.
/*            Two-Stage Least Squares .
/*===============================================.
/*variable names manipulations.
    compute fnames={"constant"}.
    compute nmx = nms(1,2:ncol(mat)).
    compute nmz = nms2(1,1:ncol(matz)).
    compute fnames={fnames,nmx}.
/*matrix containing external instruments.
/*external instruments (MExo) are instruments that do not act as their own instruments.
/*MEndo contains endogeneous variables.
   Compute Unexo={"constant"}.
    Compute MExo={ones}.
    LOOP i=2 to ncol(z).
    compute NEqual=0.
        LOOP j=2 TO k. 
            Compute D=Z(:,i)-X(:,j).
                    DO IF CSSQ(D) > 0.
                        compute dummy=1.                 
                    ELSE IF CSSQ(D) = 0.
                        compute dummy=0.
                        END IF.      
             Compute NEqual=NEqual+dummy.
        END LOOP.
        Do if NEqual=k-1.    
            Compute MExo={MExo,Z(:,i)}.
            Compute fnames={fnames,nmz(i-1)}.
            Compute Unexo={unexo,nmz(i-1)}.
        End if.  
    END LOOP.
    compute MExo=MExo.
    /*........
    /*combining matrix iv+exo.
    compute G ={X,Mexo(:,2:ncol(MExo))}.
/*------extracting endogeneous variables.
   Compute Ennames={"constant"}.
    Compute MEndo={ones}.
    LOOP i=2 to k.
    compute NEqual=0.
        LOOP j=2 TO ncol(z). 
            Compute D=X(:,i)-Z(:,j).
                    DO IF CSSQ(D) > 0.
                        compute dummy=1.                 
                    ELSE IF CSSQ(D) = 0.
                        compute dummy=0.
                        END IF.      
             Compute NEqual=NEqual+dummy.
        END LOOP.
        Do if NEqual=ncol(z)-1.    
            Compute MEndo={MEndo,x(:,i)}.
            Compute Ennames={Ennames,nmx(i-1)}.
        End if.  
    END LOOP.
    compute MEndo=MEndo.
    compute Ennames=Ennames(2:ncol(Ennames)).
/*print/title="ali baba88".
/*print z/format=f9.3.
/*print MEndo/format=f9.3.
/*------Extracting Internal instruments variables--MINST.
Compute Minst={ones}. 
    LOOP i=1 to k-1.
    compute #Equal=0.
        LOOP j=1 TO ncol(z)-1. 
              Compute D=Z(:,j+1)-X(:,i+1).
                    DO IF CSSQ(D) > 0.
                        compute dummy=0.                 
                    ELSE IF CSSQ(D) = 0.
                        compute dummy=1.
                    END IF.      
             Compute #Equal=#Equal+dummy.
        END LOOP.
        Do if #Equal>0.   
            Compute MINST={MINST,X(:,i+1)}.      
        End if. 
    END LOOP.
compute MINST=MINST.
/*print MINST/title="MINST"/format=F9.3.
/*test of order of condition.
Do IF (ncol(Z) >= ncol(X)).
        *===computing residuals, standard error of b, t value and p-value of 2SLS ==.
         compute PZ=Z*inv(sscp(Z))*t(Z).
        compute biv=inv(t(X)*PZ*X)*t(X)*PZ*Y.            
        compute eiv=Y-X*biv.
        compute res_iv=meanY-PZ*X*biv.  
        compute ssregiv=csum(sscp(res_iv)).      
        compute sseriv=csum(sscp(eiv)).
        compute mseiv=(1/(n-k))*sseriv.
        /*vbiv uses a simplied formula, see CT, p.43, assuming errors are homoskedastic.
        compute vbiv=mseiv*inv(t(X)*PZ*X).
        compute sbiv=sqrt(diag(vbiv)).
        compute tbiv=biv/sbiv.
        compute dff=n-k.
        compute F=tbiv&*tbiv.
        compute pF=1-fcdf(F,1,dff).
        compute pF=1-fcdf(F,1,dff).
             *--95% CI--.
        compute LB=biv-1.96*sbiv.
        compute UB=biv+1.96*sbiv.
/*=========bp test============.
/*************.
/*Breusch-Pagan test for heteroskedasticity.
DO IF !bptest=1.
compute e=eiv.
compute var_e=sscp(e)/n.
*residuals are scaled.
compute g=e2/var_e.
compute bp=inv(t(X)*PZ*X)*t(X)*PZ*g.
compute ep=g-X*bp.
compute e2p=ep(:,1)&*ep(:,1).
compute sserp=csum(e2p).
print/title="============================================".
print/title="Breusch-Pagan and Koenker tests using 2SLS residuals".
print/title="============================================".
print/title="The tests use scaled residuals from the 2SLS regression.".   
*--Computing LM statistics .
compute meanY=ones*t(csum(g)/n).
compute e_regp=X*bp-meanY.
compute ssregp=csum(sscp(e_regp)).
Compute total=sserp+ssregp.
Compute Rsqp=ssregp/total.
compute F=(ssregp/(k-1))/(sserp/(n-k)).
Compute pF=1-fcdf(Fval,k-1,n-k).
/* test statistics by Breusch-Pagan.
compute np=ncol(mat)-1.
Compute LMb=0.5*ssregp.
compute sigb=1-chicdf(LMb,np).
/* test statistics by Koenker.
Compute LMk=n*Rsqp.
compute sigk=1-chicdf(LMk,np).
compute LM=T({LMb,LMk}).
compute sig=T({sigb,sigk}).
print{LM,sig}
  /title '------- Breusch-Pagan and Koenker test statistics and sig-values --------'
  /clabel "LM"  "Sig"
  /rlabel "BP" "Koenker"
  /format f10.3 .
print/title="Null hypothesis: heteroskedasticity not present (homoskedasticity).".
print/title="If sig-value less than 0.05, reject the null hypothesis.". 
print/title="Note: Breusch-Pagan test is a large sample test and assumes the residuals to be "+
    "normally distributed.".
END IF.
/*======end of bp test===========.
*for output with robust std errors. 
* No correction (default).
do if (!robse=99).
compute  vbh=vbiv.
   end if.
*HC0.  
do if (!robse=0).
    compute  vbh=inv(t(X)*PZ*X)*t(X)*PZ*mdiag(e2)*PZ*X*inv(t(X)*PZ*X).
end if.
* HC1.
do if (!robse=1).
compute  vbh=inv(t(X)*PZ*X)*t(X)*PZ*mdiag(e2)*PZ*X*inv(t(X)*PZ*X).
    compute vbh=vbh*N/(N-k).
end if.
*HC2.
do if (!robse=2).
compute hat=PZ*X*inv(t(X)*PZ*X)* t(X)*PZ.
compute dhat=e2&/(ones-diag(hat)).
compute vbh=inv(t(X)*PZ*X)*t(X)*PZ*mdiag(dhat)*PZ*X*inv(t(X)*PZ*X).
end if.
*HC3.
do if (!robse=3).
compute hat=PZ*X*inv(t(X)*PZ*X)* t(X)*PZ.
compute hat2=(ones-diag(hat))&*(ones-diag(hat)).
compute dhat=e2&/hat2.
compute vbh=inv(t(X)*PZ*X)*t(X)*PZ*mdiag(dhat)*PZ*X*inv(t(X)*PZ*X).
end if.
*HC4.
do if (!robse=4).
compute hat=PZ*X*inv(t(X)*PZ*X)* t(X)*PZ.
compute fours=make(n,1,4).
compute mh={fours,n*diag(hat)/k}.
compute dummy=rmin(mh).
compute hat2=(ones-diag(hat))&**dummy.
compute dhat=e2&/hat2.
compute vbh=inv(t(X)*PZ*X)*t(X)*PZ*mdiag(dhat)*PZ*X*inv(t(X)*PZ*X).
end if.
    compute sbh=sqrt(diag(vbh)).
    compute tbh=biv/sbh.
    compute dff=n-k.
    compute Fh=tbh&*tbh.
    compute pFh=1-fcdf(Fh,1,dff).
    compute pF=1-fcdf(Fh,1,dff).
     *--95% CI--.
   compute LBh=biv-1.96*sbh.
   compute UBh=biv+1.96*sbh.
   compute tslsout={biv,sbh, tbh,pFh, LBh, UBh}.
*end of calculation.
        print/title=" ".
        print/title="=========================================".
        print/title="Two-Stage Least Squares".  
        print n/title="Sample size".      
        compute temp=t(nms(:,1)).
        print temp/title="Dependent Variable"
        / format=A10.
        print Ennames/title="Endogeneous (instrumented) variable "
            / format=A10.
        compute temp=nms2(1:ncol(matZ)).
        print temp/title="Instrument(s)"
        / format=A10.
        print Unexo(2:ncol(Unexo))/title="External instrument(s)"
        / format=A10.
       print /title="Note: External instrument(s) is truly exogeneous, not part of the original "+
    "model.".
        *===Preparing input ANOVA table of TSLS.
        *computing mean square regression.
        computing sst=sscp(y-meanY).
        /*compute ssregiv=sst-sseriv.
         compute sumsq=T({ssregiv,sseriv}).
        compute dfa=T({k-1,n-k}).
        compute mse_a=sumsq/dfa.
        compute Fval=(ssregiv/(k-1))/(sseriv/(n-k)).
        Compute pF_a=1-fcdf(Fval,k-1,n-k).
        compute F_a=T({Fval,-999}).
        compute pFa=T({pF_a,-999}).
        *--computing model summary.
       compute sereg=sqrt(mse_a(2)).
       compute Rsq=1-(sseriv/sst).
       compute Rsqadj=1-((1-Rsq)*(n-1)/(n-k-1)).
       compute nmv={"SER"}.
        compute nmvars1 = {"R^2"; "Adj";nmv}.
        compute Rsqs=({Rsq,Rsqadj,sereg}).
        print Rsqs/title="Model summary"/cnames=nmvars1/format=f9.4.
                    compute rlabel={"df(m)","df(res)","F","Sig."}.
                    compute out=({t(dfa),F_a(1),pFa(1)}).
                    print out/title=" "
                        /cnames=rlabel
                         /format=f9.4. 
        print /title="Note: Adj = Adjusted R^2, SER = Std. error of regression.".
        /*----calculate standard error of residuals....
        /*compute sd = sqrt(1/n*sser).
        /*print sd/title="Residual standard error"/format=f9.3.
        /*print {sumsq,dfa, mse_a,F_a,pFa} /space=3
          /*title '------- ANOVA TABLE --------'
          /*clabel "SS" "df" "MS" "F" "Sig"
          /*rlabel "Model" "Residual"
          /*format f10.3 .
        /*compute tslsout={biv,sbiv, tbiv,pF,LB,UB}.
         /*   compute tslsout=olsouth.
        print tslsout/title ="TSLS outputs"/rnames=nmvars/cnames=cnms/format=F9.3.
        print/title="-----------------------------------------".
do if (!robse=99).
print/title="* Note: standard errors are assumed to be homoskedastic--no adjustments.".
end if.
do if (!robse=0).
print/title="* Note: standard errors are adjusted using HC0 variant (Eicker-Huber–White standard "+
    "errors).".
print/title="* The adjustments are not recommended for sample sizes < 250 (Long and Ervin, 2000).".
end if.
do if (!robse=1).
print/title="* Note: standard errors are adjusted using HC1 variant.".
end if.
do if (!robse=2).
print/title="* Note: standard errors are adjusted using HC2 variant.".
end if.
do if (!robse=3).
print/title="* Note: standard errors are adjusted using HC3 variant.".
end if.
do if (!robse=4).
print/title="* Note: standard errors are adjusted using HC4 variant.".
end if.
       print/title="-----------------------------------------".
/*====================================.
DO IF (!reshaus=1).                
/*........Hausman specification test: Two-stage OLS regression using residuals.
                    print/title="Hausman's specification test: Two-stage OLS regression using "+
    "residuals as predictors.".
                    /*variable names manipulations.
                        compute fnames={"constant"}.
                        compute nmx = nms(1,2:ncol(mat)).                        
                        compute fnames={fnames,nmx}.
                        compute pnames={"resid1", "resid2","resid3","resid4"+
                            "resid5","resid6","resid7","resid8","resid9","resid10"}.
                    Compute G=X.
/*print/title="ali baba".
/*print MEndo/format=f9.3.
                    compute resid=make(n,ncol(MEndo)-1,99).
                    LOOP i=1 to ncol(MEndo)-1.
                           compute YEndo=MEndo(:,i+1).
                           compute bh=(inv(sscp(z)))*t(z)*YEndo.
                           compute resid(:,i)=YEndo-z*bh.
                           compute G={G,resid(:,i)}.
                           compute fnames={fnames,pnames(i)}.
                    END LOOP.
                   Release YEndo.
                    Release bh.
                    /* OLS.
 
                    print/title="Note: Results obtained via two stages.". 
                    print/title="Stage 1: Endogeneous variables are regressed on the IVs"+
                     " and external instruments.". 
                    print/title="Stage 2: DV is regressed"+
                    " on IVs and residuals obtained from stage 1.".
                    print/title="The test assumes that the errors in stage 2 are homoskedastic.".
 compute temp=t(nms(:,1)).
 OLSMACRO Yols=Y/Xols=G/numrow=n/numcol=ncol(G)/xnames=fnames/ynames=temp.
/*=====Join test for resid variables======.
                     
                    print/title="----".
                    compute df1=ncol(resid).
                    compute df2=n-ncol(G).
                    compute Fstat=((Rsq-Rsqo)/(1-Rsqo))*((n-k)/ncol(resid)).
                    compute PF=1-fcdf(Fstat,df1,df2).
                    compute rlabel={"df1","df2","F","Sig."}.
                    compute out=({df1,df2,Fstat,PF}).
                    print out/title="Joint F test for the significance of the 'resid' variables"
                        /cnames=rlabel
                         /format=f9.3. 
                    print/title="H0: Instrumented variables are exogeneous (OLS is efficient)"+
                   " (i.e., estimates of all resids are equal to zero).".
                   print/title="H1: Instrumented variables are endogeneous.".
                    print/title="-----------------------------------------".
END IF.
DO IF (!overz)=1.
/*--- J-Statistics.           
            compute bj=(inv(sscp(Z)))*t(Z)*eiv.           
            compute ej=eiv-Z*bj.
            compute ebar=ones*t(csum(ej))/n.        
            compute sser=csum(sscp(ej)).
            compute ssreg=csum(sscp(ebar-Z*bj)).
            compute Rsq=ssreg/(sser+ssreg).
            compute Fj=(Rsq/(1-Rsq))*((n-k)/(n-1)).           
            compute Jstat=n*Fj.
            do if (ncol(Z)=ncol(X)).
                  compute dfj=0.01. 
           else if (ncol(Z) > ncol(X)).
                 compute dfj=ncol(Z)-ncol(X).
            end if.
            compute PFj=1-chicdf(Jstat,dfj).
            compute rlabel={"Value", "Sig."}.
            compute out=t({Jstat,PFj}).
            print out/title="Overidentifying restrictions test (The J-statistic)"
                 /rnames=rlabel/
                 /format=f9.3.
             print /title="H0: All instruments are exogenous.".
             print/title="The test assumes a large sample size, strong instruments, and "+
    "homoskedastic errors.".
     print/title="-----------------------------------------".
*=== end of J-statistics===.
END IF.
Else if (ncol(Z) < ncol(X)).
print/ title="==================================================".
print/title="Warning: Number of instruments should be at least equal to the number of endogenous "+
    "variables".
prin/title="Identification condition is not satisfied.".
print/title="Thus, Two-stage least square procedure is terminated.".
End if.
*Weak instrument test.
/*============WEAK INSTRUMENT TEST=============.
compute Endov=MENDO(:,2:ncol(MENDO)).
compute ExZ=MEXO(:,2:ncol(MEXO)).
DO IF (!wktest=1).
print /title="Weak instrument tests".
/*Call Cragg .
CraggD Y=Endov/X=MINST/Z=ExZ.
compute sizec={0.10,0.15,0.20,0.25}.
compute textc={"Size","Values"}.
compute sizec={0.10,0.15,0.20,0.25}.
compute textc={"Size","Values"}.
compute whichrow=ncol(ExZ).
/*Weak Instrument Test Based on TSLS Size.
 Do IF ncol(ENDOV)=1.
/*for n=2.
Compute CVS={16.38,8.96,6.66,5.53;19.93,11.59,8.75,7.25;22.30,12.83,9.54,7.80;+
                       24.58,13.96,10.26,8.31;26.87,15.09,10.98,8.84;29.18,16.23,11.72,9.38;+
                       31.50,17.38,12.48,9.93;33.84,18.54,13.24,10.50;36.19,19.71,14.01,11.07;+
                       38.54,20.88,14.78,11.65;40.90,22.06,15.56,12.23;43.27,23.24,16.35,12.82;+
                       45.64,24.42,17.14,13.41;48.01,25.61,17.93,14.00;50.39,26.80,18.72,14.60;+
                        52.77,27.99,19.51,15.19;55.15,29.19,20.31,15.79;57.53,30.38,21.10,16.39}.
compute whrow={sizec;CVS(whichrow,:)}.
print t(whrow)/title="Stock-Yogo critical values for the weak instrument test based on TSLS size" +
    " at 5% significance level:"/cnames=textc/format=f9.2.
ELSE IF ncol(ENDOV)=2 AND ncol(ExZ)>1.
Compute CVS={0,0,0,0;7.03,4.58,3.95,3.63;13.43,8.18,6.40,5.45;16.87,9.93,7.54,6.28; +
                        19.45,11.22,8.38,6.89;21.68,12.33,9.10,7.42;23.72,13.34,9.77,7.91; +
                        23.72,13.34,9.77,7.91;25.64,14.31,10.41,8.3;27.51,15.24,11.03,8.85; +
                        29.32,16.16,11.65,9.31;31.11, 17.06,12.25,9.77;32.88,17.95,12.86,10.22;+
                        34.62,18.84,13.45,10.68;36.36,19.72,14.05,11.13;38.08,20.60,14.65,11.58;+
                        39.80,21.48,15.24,12.03;41.51,22.35,15.83,12.49;43.22,23.22,16.42,12.94;+
                        44.92,24.09,17.02,13.39;46.62,24.96,17.61,13.84;48.31,25.82,18.20,14.29;+
                        50.01,26.69,18.79,14.74;51.70,27.56,19.38,15.19;53.39,28.42,19.97,15.64;+
                        55.07,29.29,20.56,16.10;56.76,30.15,21.15,16.55;58.45,31.02,21.74,17.00 ;+
                        60.13,31.88,22.33,17.45;61.82,32.74,22.92,17.90;63.51,33.61,23.51,18.35;+
                        59.92,31.58,21.90,16.99;62.30,32.77,22.70,17.60;64.69,33.97,23.50,18.20;+
                        67.07,35.17,24.30,18.80;69.46,36.37,25.10,19.41;71.85,37.57,25.90,20.01;+
                        74.24,38.77,26.71,20.61;76.62,39.97,27.51,21.22;79.01,41.17,28.31,21.83;+
                        81.40,42.37,29.12,22.43;83.79,43.57,29.92,23.04;86.17,44.78,30.72,23.65}.
compute whrow={sizec;CVS(whichrow,:)}.
print t(whrow)/title="Stock-Yogo critical values for the weak instrument test based on TSLS size" +
    " at 5% significance level:"/cnames=textc/format=f9.2.
END IF.
/*------.
/*Weak Instrument Test Based on TSLS Bias.
/*for n=1, k2 at least =3.
compute whichrow=whichrow-2.
Do IF ncol(ENDOV)=1 AND ncol(ExZ)>2.
Compute CVBS={13.91,9.08,6.46,5.39;16.85,10.27,6.71,5.34;18.37,10.83,6.77,5.25;19.28,11.12,6.76,
    5.15;+
                          19.86,11.29,6.73,5.07;20.25,11.39,6.69,4.99;20.53,11.46,6.65,4.92;20.74,
    11.49,6.61,4.86;+
                          20.90,11.51,6.56,4.80;21.01,11.52,6.53,4.75;21.10,11.52,6.49,4.71;21.28,
    11.52,6.45,4.67;+
                          21.23,11.51,6.42,4.63;21.28,11.50,6.39,4.59;21.31,11.49,6.36,4.56;21.34,
    11.48,6.33,4.53;+
                          21.36,11.46,6.31,4.51;21.38,11.45,6.28,4.48;21.39,11.44,6.26,4.46;21.40,
    11.42,6.24,4.43;+
                          21.41,11.41,6.22,4.41;21.41,11.40,6.20,4.39;21.42,11.38,6.18,4.37;21.42,
    11.37,6.16,4.35;+
                          21.42,11.36,6.14,4.34;21.42,11.34,6.13,4.32;21.42,11.33,6.11,4.31;21.42,
    11.32,6.09,4.29}.
compute whrow={sizec;CVBS(whichrow,:)}.
print t(whrow)/title="Stock-Yogo critical values for the weak instrument test based on TSLS "+
    "relative bias" +
    " at 5% significance level:"/cnames=textc/format=f9.2.
ELSE IF ncol(ENDOV)=1 AND ncol(ExZ)>3.
Compute CVBS={11.04,7.56,5.57,4.73;13.97,8.78,5.91,4.79;15.72,9.48,6.08,4.78;16.88,9.92,6.16,4.76;+
                          17.70,10.22,6.20,4.73;18.30,10.43,6.22,4.69;18.76,10.58,6.23,4.66;19.12,
    10.69,6.23,4.62;+
                          19.40,10.78,6.22,4.59;19.64,10.84,6.21,4.56;19.83,10.89,6.20,4.53;19.98,
    10.93,6.19,4.50;+
                          20.12,10.96,6.17,4.48;20.23,10.99,6.16,4.45;20.33,11.00,6.14,4.43;20.41,
    11.02,6.13,4.41;+
                          20.48,11.03,6.11,4.39;20.54,11.04,6.10,4.37;20.60,11.05,6.08,4.35;20.65,
    11.05,6.07,4.33;+
                          20.69,11.05,6.06,4.32;20.73,11.06,6.05,4.30;20.76,11.06,6.03,4.29;20.79,
    11.06,6.02,4.27;+
                          20.82,11.05,6.01,4.26;20.84,11.05,6.00,4.24;20.86,11.05,5.99,4.23}.
compute whrow={sizec;CVBS(whichrow,:)}.
print t(whrow)/title="Stock-Yogo critical values for the weak instrument test based on TSLS "+
    "relative bias" +
    " at 5% significance level:"/cnames=textc/format=f9.2.
ELSE IF ncol(ENDOV)=1 AND ncol(ExZ)>4.
Compute CVBS={9.53,6.61,4.99,4.30;12.20,7.77,5.35,4.40;13.95,8.50,5.56,4.44;15.18,9.01,5.69,4.46;+
                         16.10,9.37,5.78,4.46;16.80,9.64,5.83,4.45;17.35,9.85,5.87,4.44;17.80,10.01,
    5.90,4.42;+
                         18.17,10.14,5.92,4.41;18.47,10.25,5.93,4.39;18.73,10.33,5.94,4.37;18.94,
    10.41,5.94,4.36;+
                         19.13,10.47,5.94,4.34;19.29,10.52,5.94,4.32;19.44,10.56,5.94,4.31;19.56,
    10.60,5.93,4.29;+
                         19.67,10.63,5.93,4.28;19.77,10.65,5.92,4.27;19.86,10.68,5.92,4.25;19.94,
    10.70,5.91,4.24;+
                         20.01,10.71,5.90,4.23;20.07,10.73,5.90,4.21;20.13,10.74,5.89,4.20;20.18,
    10.75,5.88,4.19;
                         20.23,10.76,5.88,4.18;20.27,10.77,5.87,4.17}.
compute whrow={sizec;CVBS(whichrow,:)}.
print t(whrow)/title="Stock-Yogo critical values for the weak instrument test based on TSLS "+
    "relative bias" +
    " at 5% significance level:"/cnames=textc/format=f9.2.
ELSE IF ncol(ExZ)<3.
print /title="Stock-Yogo critical values for relative bias are not available for models with "+
    "less than 3 instruments.".
END IF.
DO IF (ncol(ENDOV)=1).
print/ title="Note: For one endogeneous variable, Cragg-Donald F-statistic is the"+
                   " F-value of first stage regression. ".
print/ title="For one endogeneous regressor, instruments are weak" +
                    " if the F-statistics is less than 10 (Staiger and Stock, 1997).".
END IF.

print/title="Note: Tables for the critical values (Stock & Yogo, 2002) are reproduced by "+
    "permission--"+
    "email communication (Stock, 16-11-2019).".
print /title="==============================".
END IF.

/*=========END OF WEAK INSTRUMENT TEST.
Release MEXO.
Release MENDO.
Release X.
Release Z.
END MATRIX.
!ENDDEFINE.
restore
