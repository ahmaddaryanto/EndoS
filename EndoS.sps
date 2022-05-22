* Encoding: UTF-8.
/* ENDoS v2 for ENDogeneity analysis.
/*Macro created by Ahmad Daryanto*/.
/* March 2022.
/* Update: 
/* omit missing cases.
/* CIs for TSLS estimates corrected.

SET WORKSPACE=134576.
DEFINE OLSMACRO (Xols = !charEND('/')
/Yols = !charEND('/')
/numcol = !charEND('/')
/numrow = !charEND('/')
/xnames = !charEND('/')
/ynames = !charEND('/')).
/*OLS.
COMPUTE nmvars={!xnames}.
PRINT !ynames/Title="Dependent variable"/format=A9.
COMPUTE bols=(inv(sscp(!Xols)))*t(!Xols)*!Yols.
COMPUTE e=!Yols-!Xols*bols.
COMPUTE e2=e(:,1)&*e(:,1).
COMPUTE sser=csum(sscp(e)).
COMPUTE mse=(1/(!numrow-!numcol))*sser.
COMPUTE vbols=mse*inv(sscp(!Xols)).
COMPUTE sbols=sqrt(diag(vbols)).
COMPUTE tbols=bols/sbols.
COMPUTE dffols=!numrow-!numcol.
COMPUTE Fols=tbols&*tbols.
COMPUTE pFols=1-fcdf(Fols,1,dffols).
     *--95% CI--.
COMPUTE LBols=bols-1.96*sbols.
COMPUTE UBols=bols+1.96*sbols.
COMPUTE olsout={bols,sbols, tbols,pFols, LBols, UBols}.
*END of calculation.
*===Preparing input ANOVA table of OLS.
*computing mean square regression.
COMPUTE meanY=ones*t(csum(!Yols)/!numrow).
COMPUTE e_reg=meanY-!Xols*bols.
COMPUTE ssreg=csum(sscp(e_reg)).
COMPUTE sumsq=T({ssreg,sser}).
COMPUTE dfa=T({!numcol-1,!numrow-!numcol}).
COMPUTE mse_a=sumsq/dfa.
COMPUTE Fval=(ssreg/(!numcol-1))/(sser/(!numrow-!numcol)).
COMPUTE pF_a=1-fcdf(Fval,!numcol-1,!numrow-!numcol).
COMPUTE F_a=T({Fval,-999}).
COMPUTE pFa=T({pF_a,-999}).
*--computing model summary.
COMPUTE sereg=sqrt(mse_a(2)).
COMPUTE Rsq=ssreg/(sser+ssreg).
COMPUTE Rsqadj=1-((1-Rsq)*(!numrow-1)/(!numrow-!numcol-1)).
        COMPUTE Rsqs=({Rsq,Rsqadj,sereg}).
        PRINT Rsqs/title="Model summary"/cnames=nmvars1/format=f9.4.
                    COMPUTE rlabel={"df(m)","df(res)","F","Sig."}.
                    COMPUTE out=({t(dfa),F_a(1),pFa(1)}).
                     PRINT out/title=" "/cnames=rlabel/format=f9.4. 
        *--OLS output.
      
        PRINT /title="Note: Adj = Adjusted R^2, SER = Std. error of regression.".
        PRINT olsout/title ="======="/rnames=nmvars/cnames=cnms/format=F9.3.
!ENDDEFINE.
/*&&&&&&&&&&&&&&&&&&&&.


DEFINE CraggD (Y = !charEND('/')
/X = !charEND('/')
/Z = !charEND('/')).
COMPUTE  #Y=!Y.
COMPUTE #X=!X.
COMPUTE #Z=!Z.
COMPUTE Indt=IDENT(nrow(#Z),nrow(#Z)).
COMPUTE Z_={#X,#Z}.
COMPUTE MZ_=(Indt-Z_*(inv(sscp(Z_)))*t(Z_)).
COMPUTE MX=(Indt-#X*(inv(sscp(#X)))*t(#X)).
COMPUTE Y_T=MX*#Y.
COMPUTE Z_T=MX*#Z.
COMPUTE PZ_T=Z_T*(inv(sscp(Z_T)))*t(Z_T).
COMPUTE numT=nrow(#Z).
COMPUTE K1=ncol(#X).
COMPUTE K2=ncol(#Z).
COMPUTE Sigv=(t(#Y)*MZ_*#Y)/(numT-K2-K1).
COMPUTE sqrSigv=CHOL(Sigv).
COMPUTE invSigv=inv(sqrSigv).
/*COMPUTE fac=(numT-K2-K1)/K2.
COMPUTE matG=t(invSigv) *t(Y_T) * PZ_T *Y_T* invSigv/K2.
COMPUTE gmin=cmin(Eval(matG)).
PRINT gmin/title="Cragg-Donald F-statistic"/format=f9.4.
!ENDDEFINE.

/*&&&&&&&&&&&&&&&&&&&&.
/*&&&&&&&&&&&&&&&&&&&&.
DEFINE TSLS (iv = !charEND('/')
/dv = !charEND('/')
/ivm=!charEND('/')
/olsoutp=!charEND('/')!default(0)
/reshaus=!charEND('/')!default(0)
/overz=!charEND('/')!default(0)
/bptest=!charEND('/')!default(0)
/wktest=!charEND('/')!default(0)
/robse=!charEND('/')!default(99)).
SET MXLOOPS = 10000001.
SET PRINTBACK = OFF.
MATRIX.
GET mat/variables=!dv !iv /names=nms /MISSING=-999.999.
GET matz/variables=!ivm/names=nms2/MISSING=-999.999.

/*================.
/* Retaining complete cases.

COMPUTE rawmat={mat,matz}.

COMPUTE k=0.
LOOP i=1 to nrow(rawmat).
    COMPUTE ENDofrow=0.
     COMPUTE flag=0.
      LOOP j=1 to ncol(rawmat).
          DO IF rawmat(i,j) = -999.999.
             COMPUTE flag=1.              
          END IF.
      END LOOP.
      COMPUTE ENDofrow=1.
      
        /* counting
        DO IF ((flag=0) AND (ENDofrow=1)). 
            COMPUTE k=k+1.
         END IF.
END LOOP.

/* prepare the matrix to store complete cases.
COMPUTE mcomp=make(k,ncol(rawmat),0).

COMPUTE k=1.
/* Store non-missing.
LOOP i=1 to nrow(rawmat).
     COMPUTE ENDofrow=0.
     COMPUTE flag=0.
      LOOP j=1 to ncol(rawmat).
          DO IF rawmat(i,j) = -999.999.
               COMPUTE flag=1.  
          END IF.
      END LOOP.
      COMPUTE ENDofrow=1.
      
        /* 
        DO IF ((flag=0) AND (ENDofrow=1)). 
           COMPUTE mcomp(k,:)=rawmat(i,:).
            COMPUTE k=k+1.
         END IF.
                
        /*

END LOOP.

/*====Preparing matrix==.

PRINT/title="=========================================".
COMPUTE nomit=nrow(rawmat)-nrow(mcomp).
PRINT nomit/title="Number of missing cases (omitted in the analysis)".
PRINT/title="=========================================".

COMPUTE n=nrow(mcomp).     
COMPUTE ones=make(n,1,1).
COMPUTE Y=mcomp(:,1).
COMPUTE X={ones,mcomp(:,2:ncol(mat))}.
COMPUTE Z={ones,mcomp(:,(ncol(mat)+1):ncol(mcomp))}.
COMPUTE k=ncol(X).
/*======.

PRINT/ title="Two-Stage Least Squares".
PRINT/title="written by Ahmad Daryanto".
PRINT/title="https://sites.google.com/site/ahmaddaryanto/".
PRINT/title="cite: Daryanto, A. (2020). ENDoS: An SPSS macro to assess ENDogeneity. The Quantitative Methods for Psychology, 16(1), 56-70.".
PRINT/title="-----------------------------------------".
COMPUTE dvname=t(nms(:,1)).
COMPUTE ivnames=t(nms(:,2:k)).
    DO IF (!olsoutp=1).
    PRINT/title="OLS regression".
    PRINT n/title="Sample size (complete cases)f".
    PRINT dvname/title="Dependent variable (DV)"/ format=A10.
    PRINT t(ivnames)/title="Independent variable (IVs)"
            / format=A10.
    END IF.

*===============================================.
*            OLS Regression of Without Instruments .
*===============================================.
*dv in original metrix.
COMPUTE b=(inv(sscp(X)))*t(X)*Y.
/*---.
/*notes: sser=sum squares of error.ssreg=sum of squares of regression.
*===computing residuals, standard error of b, t value and p-value of OLS ==.
COMPUTE e=Y-X*b.
COMPUTE e2=e(:,1)&*e(:,1).
COMPUTE sser=csum(sscp(e)).
COMPUTE mse=(1/(n-k))*sser.
COMPUTE vb=mse*inv(sscp(X)).
COMPUTE sb=sqrt(diag(vb)).
COMPUTE tb=b/sb.
COMPUTE dff=n-k.
COMPUTE F=tb&*tb.
COMPUTE pF=1-fcdf(F,1,dff).
COMPUTE pF=1-fcdf(F,1,dff).
     *--95% CI--.
COMPUTE LB=b-1.96*sb.
COMPUTE UB=b+1.96*sb.
*for output with robust std errors. 
* No correction (default).
DO IF (!robse=99).
COMPUTE  vbh=vb.
   END IF.
*HC0.  
DO IF (!robse=0).
    COMPUTE  vbh=inv(sscp(X))*t(X)*mdiag(e2)*X*inv(sscp(X)).
END IF.
* HC1.
DO IF (!robse=1).
COMPUTE  vbh=inv(sscp(X))*t(X)*mdiag(e2)*X*inv(sscp(X)).
    COMPUTE vbh=vbh*N/(N-k).
END IF.
*HC2.
DO IF (!robse=2).
COMPUTE hat=X*inv(sscp(X))* t(X).
COMPUTE dhat=e2&/(ones-diag(hat)).
COMPUTE vbh=inv(sscp(X))*t(X)*mdiag(dhat)*X*inv(sscp(X)).
END IF.
*HC3.
DO IF (!robse=3).
COMPUTE hat=X*inv(sscp(X))* t(X).
COMPUTE hat2=(ones-diag(hat))&*(ones-diag(hat)).
COMPUTE dhat=e2&/hat2.
COMPUTE vbh=inv(sscp(X))*t(X)*mdiag(dhat)*X*inv(sscp(X)).
END IF.
*HC4.
DO IF (!robse=4).
COMPUTE hat=X*inv(sscp(X))* t(X).
COMPUTE fours=make(n,1,4).
COMPUTE mh={fours,n*diag(hat)/k}.
COMPUTE dummy=rmin(mh).
COMPUTE hat2=(ones-diag(hat))&**dummy.
COMPUTE dhat=e2&/hat2.
COMPUTE vbh=inv(sscp(X))*t(X)*mdiag(dhat)*X*inv(sscp(X)).
END IF.
    COMPUTE sbh=sqrt(diag(vbh)).
    COMPUTE tbh=b/sbh.
    COMPUTE dff=n-k.
    COMPUTE Fh=tbh&*tbh.
    COMPUTE pFh=1-fcdf(Fh,1,dff).
    COMPUTE pF=1-fcdf(Fh,1,dff).
     *--95% CI--.
   COMPUTE LBh=b-1.96*sb.
   COMPUTE UBh=b+1.96*sb.
   COMPUTE olsouth={b,sbh, tbh,pFh, LBh, UBh}.
*END of calculation.
*===Preparing input ANOVA table of OLS.
*computing mean square regression.
COMPUTE meanY=ones*t(csum(Y)/n).
COMPUTE e_reg=meanY-X*b.
COMPUTE ssreg=csum(sscp(e_reg)).
COMPUTE sumsq=T({ssreg,sser}).
COMPUTE dfa=T({k-1,n-k}).
COMPUTE mse_a=sumsq/dfa.
COMPUTE Fval=(ssreg/(k-1))/(sser/(n-k)).
COMPUTE pF_a=1-fcdf(Fval,k-1,n-k).
COMPUTE F_a=T({Fval,-999}).
COMPUTE pFa=T({pF_a,-999}).
*--computing model summary.
COMPUTE sereg=sqrt(mse_a(2)).
COMPUTE Rsqo=ssreg/(sser+ssreg).
COMPUTE Rsqadj=1-((1-Rsqo)*(n-1)/(n-k-1)).
        COMPUTE nmv={"SER"}.
        COMPUTE nmvars1 = {"R^2"; "Adj";nmv}.
        COMPUTE Rsqs=({Rsqo,Rsqadj,sereg}).
        COMPUTE nmvars = t(nms(1,2:ncol(mat))).
        COMPUTE nmvars = {"constant"; nmvars}.
        COMPUTE cnms={"b","se", "t", "sig", "95%LB", "95%UB"}.
DO IF (!olsoutp=1).
        PRINT Rsqs/title="Model summary"/cnames=nmvars1/format=f9.4.
                    COMPUTE rlabel={"df(m)","df(res)","F","Sig."}.
                    COMPUTE out=({t(dfa),F_a(1),pFa(1)}).
                    
        PRINT /title="Note: Adj = Adjusted R^2, SER = Std. error of regression.".
       /* PRINT {sumsq,dfa, mse_a,F_a,pFa} /space=3
          /*title '------- ANOVA TABLE  --------'.
          /*clabel "SS" "df" "MS" "F" "Sig".
          /*rlabel "Model" "Residual".
          /*format f10.3 .
        *--OLS output.
PRINT olsouth/title ="OLS outputs"/rnames=nmvars/cnames=cnms/format=F9.3.
*Correcting for standard errors--robust standard errors.
DO IF (!robse=99).
PRINT/title="* Note: standard errors are assumed to be homoskedastic--no adjustments.".
END IF.
DO IF (!robse=0).
PRINT/title="* Note: standard errors are adjusted using HC0 variant (Eicker-Huber–White standard "+
    "errors).".
PRINT/title="* The adjustments are not recommENDed for sample sizes < 250 (Long and Ervin, 2000).".
END IF.
DO IF (!robse=1).
PRINT/title="* Note: standard errors are adjusted using HC1 variant.".
END IF.
DO IF (!robse=2).
PRINT/title="* Note: standard errors are adjusted using HC2 variant.".
END IF.
DO IF (!robse=3).
PRINT/title="* Note: standard errors are adjusted using HC3 variant.".
END IF.
DO IF (!robse=4).
PRINT/title="* Note: standard errors are adjusted using HC4 variant.".
END IF.
        
END IF.
  
/*===============================================.
/*            Two-Stage Least Squares .
/*===============================================.
/*variable names manipulations.
    COMPUTE fnames={"constant"}.
    COMPUTE nmx = nms(1,2:ncol(mat)).
    COMPUTE nmz = nms2(1,1:ncol(matz)).
    COMPUTE fnames={fnames,nmx}.
/*matrix containing external instruments.
/*external instruments (MExo) are instruments that DO not act as their own instruments.
/*MENDo contains ENDogeneous variables.
   COMPUTE Unexo={"constant"}.
    COMPUTE MExo={ones}.
    LOOP i=2 to ncol(z).
    COMPUTE NEqual=0.
        LOOP j=2 TO k. 
            COMPUTE D=Z(:,i)-X(:,j).
                    DO IF CSSQ(D) > 0.
                        COMPUTE dummy=1.                 
                    ELSE IF CSSQ(D) = 0.
                        COMPUTE dummy=0.
                        END IF.      
             COMPUTE NEqual=NEqual+dummy.
        END LOOP.
        Do IF NEqual=k-1.    
            COMPUTE MExo={MExo,Z(:,i)}.
            COMPUTE fnames={fnames,nmz(i-1)}.
            COMPUTE Unexo={unexo,nmz(i-1)}.
        END IF.  
    END LOOP.
    COMPUTE MExo=MExo.
    /*........
    /*combining matrix iv+exo.
    COMPUTE G ={X,Mexo(:,2:ncol(MExo))}.
/*------extracting ENDogeneous variables.
   COMPUTE Ennames={"constant"}.
    COMPUTE MENDo={ones}.
    LOOP i=2 to k.
    COMPUTE NEqual=0.
        LOOP j=2 TO ncol(z). 
            COMPUTE D=X(:,i)-Z(:,j).
                    DO IF CSSQ(D) > 0.
                        COMPUTE dummy=1.                 
                    ELSE IF CSSQ(D) = 0.
                        COMPUTE dummy=0.
                        END IF.      
             COMPUTE NEqual=NEqual+dummy.
        END LOOP.
        Do IF NEqual=ncol(z)-1.    
            COMPUTE MENDo={MENDo,x(:,i)}.
            COMPUTE Ennames={Ennames,nmx(i-1)}.
        END IF.  
    END LOOP.
    COMPUTE MENDo=MENDo.
    COMPUTE Ennames=Ennames(2:ncol(Ennames)).
/*------Extracting Internal instruments variables--MINST.
COMPUTE Minst={ones}. 
    LOOP i=1 to k-1.
    COMPUTE #Equal=0.
        LOOP j=1 TO ncol(z)-1. 
              COMPUTE D=Z(:,j+1)-X(:,i+1).
                    DO IF CSSQ(D) > 0.
                        COMPUTE dummy=0.                 
                    ELSE IF CSSQ(D) = 0.
                        COMPUTE dummy=1.
                    END IF.      
             COMPUTE #Equal=#Equal+dummy.
        END LOOP.
        Do IF #Equal>0.   
            COMPUTE MINST={MINST,X(:,i+1)}.      
        END IF. 
    END LOOP.
COMPUTE MINST=MINST.
/*PRINT MINST/title="MINST"/format=F9.3.
/*test of order of condition.
Do IF (ncol(Z) >= ncol(X)).
        *===computing residuals, standard error of b, t value and p-value of 2SLS ==.
         COMPUTE PZ=Z*inv(sscp(Z))*t(Z).
        COMPUTE biv=inv(t(X)*PZ*X)*t(X)*PZ*Y.            
        COMPUTE eiv=Y-X*biv.
        COMPUTE res_iv=meanY-PZ*X*biv.  
        COMPUTE ssregiv=csum(sscp(res_iv)).      
        COMPUTE sseriv=csum(sscp(eiv)).
        COMPUTE mseiv=(1/(n-k))*sseriv.
        /*vbiv uses a simplied formula, see CT, p.43, assuming errors are homoskedastic.
        COMPUTE vbiv=mseiv*inv(t(X)*PZ*X).
        COMPUTE sbiv=sqrt(diag(vbiv)).
        COMPUTE tbiv=biv/sbiv.
        COMPUTE dff=n-k.
        COMPUTE F=tbiv&*tbiv.
        COMPUTE pF=1-fcdf(F,1,dff).
        COMPUTE pF=1-fcdf(F,1,dff).
             *--95% CI--.
        COMPUTE LB=biv-1.96*sbiv.
        COMPUTE UB=biv+1.96*sbiv.
/*=========bp test============.
/*************.
 
/*Breusch-Pagan test for heteroskedasticity.
DO IF !bptest=1.
COMPUTE e=eiv.
COMPUTE var_e=sscp(e)/n.
*residuals are scaled.
COMPUTE g=e2/var_e.
COMPUTE bp=inv(t(X)*PZ*X)*t(X)*PZ*g.
COMPUTE ep=g-X*bp.
COMPUTE e2p=ep(:,1)&*ep(:,1).
COMPUTE sserp=csum(e2p).
PRINT/title="============================================".
PRINT/title="Breusch-Pagan and Koenker tests using 2SLS residuals".
PRINT/title="============================================".
PRINT/title="The tests use scaled residuals from the 2SLS regression.".   
*--Computing LM statistics .
COMPUTE meanY=ones*t(csum(g)/n).
COMPUTE e_regp=X*bp-meanY.
COMPUTE ssregp=csum(sscp(e_regp)).
COMPUTE total=sserp+ssregp.
COMPUTE Rsqp=ssregp/total.
COMPUTE F=(ssregp/(k-1))/(sserp/(n-k)).
COMPUTE pF=1-fcdf(Fval,k-1,n-k).
/* test statistics by Breusch-Pagan.
COMPUTE np=ncol(mat)-1.
COMPUTE LMb=0.5*ssregp.
COMPUTE sigb=1-chicdf(LMb,np).
/* test statistics by Koenker.
COMPUTE LMk=n*Rsqp.
COMPUTE sigk=1-chicdf(LMk,np).
COMPUTE LM=T({LMb,LMk}).
COMPUTE sig=T({sigb,sigk}).
PRINT{LM,sig}
  /title '------- Breusch-Pagan and Koenker test statistics and sig-values --------'
  /clabel "LM"  "Sig"
  /rlabel "BP" "Koenker"
  /format f10.3 .
PRINT/title="Null hypothesis: heteroskedasticity not present (homoskedasticity).".
PRINT/title="If sig-value less than 0.05, reject the null hypothesis.". 
PRINT/title="Note: Breusch-Pagan test is a large sample test and assumes the residuals to be "+
    "normally distributed.".
END IF.
/*======END of bp test===========.
*for output with robust std errors. 
* No correction (default).
DO IF (!robse=99).
COMPUTE  vbh=vbiv.
   END IF.
*HC0.  
DO IF (!robse=0).
    COMPUTE  vbh=inv(t(X)*PZ*X)*t(X)*PZ*mdiag(e2)*PZ*X*inv(t(X)*PZ*X).
END IF.
* HC1.
DO IF (!robse=1).
COMPUTE  vbh=inv(t(X)*PZ*X)*t(X)*PZ*mdiag(e2)*PZ*X*inv(t(X)*PZ*X).
    COMPUTE vbh=vbh*N/(N-k).
END IF.
*HC2.
DO IF (!robse=2).
COMPUTE hat=PZ*X*inv(t(X)*PZ*X)* t(X)*PZ.
COMPUTE dhat=e2&/(ones-diag(hat)).
COMPUTE vbh=inv(t(X)*PZ*X)*t(X)*PZ*mdiag(dhat)*PZ*X*inv(t(X)*PZ*X).
END IF.
*HC3.
DO IF (!robse=3).
COMPUTE hat=PZ*X*inv(t(X)*PZ*X)* t(X)*PZ.
COMPUTE hat2=(ones-diag(hat))&*(ones-diag(hat)).
COMPUTE dhat=e2&/hat2.
COMPUTE vbh=inv(t(X)*PZ*X)*t(X)*PZ*mdiag(dhat)*PZ*X*inv(t(X)*PZ*X).
END IF.
*HC4.
DO IF (!robse=4).
COMPUTE hat=PZ*X*inv(t(X)*PZ*X)* t(X)*PZ.
COMPUTE fours=make(n,1,4).
COMPUTE mh={fours,n*diag(hat)/k}.
COMPUTE dummy=rmin(mh).
COMPUTE hat2=(ones-diag(hat))&**dummy.
COMPUTE dhat=e2&/hat2.
COMPUTE vbh=inv(t(X)*PZ*X)*t(X)*PZ*mdiag(dhat)*PZ*X*inv(t(X)*PZ*X).
END IF.
    COMPUTE sbh=sqrt(diag(vbh)).
    COMPUTE tbh=biv/sbh.
    COMPUTE dff=n-k.
    COMPUTE Fh=tbh&*tbh.
    COMPUTE pFh=1-fcdf(Fh,1,dff).
    COMPUTE pF=1-fcdf(Fh,1,dff).
     *--95% CI--.
   COMPUTE LBh=biv-1.96*sbh.
   COMPUTE UBh=biv+1.96*sbh.
   COMPUTE tslsout={biv,sbh, tbh,pFh, LBh, UBh}.
*END of calculation.
        PRINT/title=" ".
        PRINT/title="=========================================".
        PRINT/title="Two-Stage Least Squares".  
        PRINT n/title="Sample size (complete cases)".      
        COMPUTE temp=t(nms(:,1)).
        PRINT temp/title="Dependent Variable"
        / format=A10.
        PRINT Ennames/title="Endogeneous (instrumented) variable "
            / format=A10.
        COMPUTE temp=nms2(1:ncol(matZ)).
        PRINT temp/title="Instrument(s)"
        / format=A10.
        PRINT Unexo(2:ncol(Unexo))/title="External instrument(s)"
        / format=A10.
       PRINT /title="Note: External instrument(s) is truly exogeneous, not part of the original "+
    "model.".
        *===Preparing input ANOVA table of TSLS.
        /*computing mean square regression.
        COMPUTE sst=sscp(y-meanY).
        /*COMPUTE ssregiv=sst-sseriv.
         COMPUTE sumsq=T({ssregiv,sseriv}).
        COMPUTE dfa=T({k-1,n-k}).
        COMPUTE mse_a=sumsq/dfa.
        COMPUTE Fval=(ssregiv/(k-1))/(sseriv/(n-k)).
        COMPUTE pF_a=1-fcdf(Fval,k-1,n-k).
        COMPUTE F_a=T({Fval,-999}).
        COMPUTE pFa=T({pF_a,-999}).
        *--computing model summary.
       COMPUTE sereg=sqrt(mse_a(2)).
       COMPUTE Rsq=1-(sseriv/sst).
       COMPUTE Rsqadj=1-((1-Rsq)*(n-1)/(n-k-1)).
       COMPUTE nmv={"SER"}.
        COMPUTE nmvars1 = {"R^2"; "Adj";nmv}.
        COMPUTE Rsqs=({Rsq,Rsqadj,sereg}).
        PRINT Rsqs/title="Model summary"/cnames=nmvars1/format=f9.4.
                    COMPUTE rlabel={"df(m)","df(res)","F","Sig."}.
                    COMPUTE out=({t(dfa),F_a(1),pFa(1)}).
                    PRINT out/title=" "
                        /cnames=rlabel
                         /format=f9.4. 
        PRINT /title="Note: Adj = Adjusted R^2, SER = Std. error of regression.".
        /*----calculate standard error of residuals....
        /*COMPUTE sd = sqrt(1/n*sser).
        /*PRINT sd/title="Residual standard error"/format=f9.3.
        /*PRINT {sumsq,dfa, mse_a,F_a,pFa} /space=3
          /*title '------- ANOVA TABLE --------'
          /*clabel "SS" "df" "MS" "F" "Sig"
          /*rlabel "Model" "Residual"
          /*format f10.3 .
        /*COMPUTE tslsout={biv,sbiv, tbiv,pF,LB,UB}.
         /*   COMPUTE tslsout=olsouth.
        PRINT tslsout/title ="TSLS outputs"/rnames=nmvars/cnames=cnms/format=F9.3.
        PRINT/title="-----------------------------------------".
DO IF (!robse=99).
PRINT/title="* Note: standard errors are assumed to be homoskedastic--no adjustments.".
END IF.
DO IF (!robse=0).
PRINT/title="* Note: standard errors are adjusted using HC0 variant (Eicker-Huber–White standard "+
    "errors).".
PRINT/title="* The adjustments are not recommENDed for sample sizes < 250 (Long and Ervin, 2000).".
END IF.
DO IF (!robse=1).
PRINT/title="* Note: standard errors are adjusted using HC1 variant.".
END IF.
DO IF (!robse=2).
PRINT/title="* Note: standard errors are adjusted using HC2 variant.".
END IF.
DO IF (!robse=3).
PRINT/title="* Note: standard errors are adjusted using HC3 variant.".
END IF.
DO IF (!robse=4).
PRINT/title="* Note: standard errors are adjusted using HC4 variant.".
END IF.
       PRINT/title="-----------------------------------------".
/*====================================.
DO IF (!reshaus=1).                
/*........Hausman specIFication test: Two-stage OLS regression using residuals.
                    PRINT/title="Hausman's specIFication test: Two-stage OLS regression using "+
    "residuals as predictors.".
                    /*variable names manipulations.
                        COMPUTE fnames={"constant"}.
                        COMPUTE nmx = nms(1,2:ncol(mat)).                        
                        COMPUTE fnames={fnames,nmx}.
                        COMPUTE pnames={"resid1", "resid2","resid3","resid4"+
                            "resid5","resid6","resid7","resid8","resid9","resid10"}.
                    COMPUTE G=X.
/*PRINT MENDo/format=f9.3.
                    COMPUTE resid=make(n,ncol(MENDo)-1,99).
                    LOOP i=1 to ncol(MENDo)-1.
                           COMPUTE YENDo=MENDo(:,i+1).
                           COMPUTE bh=(inv(sscp(z)))*t(z)*YENDo.
                           COMPUTE resid(:,i)=YENDo-z*bh.
                           COMPUTE G={G,resid(:,i)}.
                           COMPUTE fnames={fnames,pnames(i)}.
                    END LOOP.
                   RELEASE YENDo.
                    RELEASE bh.
                    /* OLS.
 
                    PRINT/title="Note: Results obtained via two stages.". 
                    PRINT/title="Stage 1: ENDogeneous variables are regressed on the IVs"+
                     " and external instruments.". 
                    PRINT/title="Stage 2: DV is regressed"+
                    " on IVs and residuals obtained from stage 1.".
                    PRINT/title="The test assumes that the errors in stage 2 are homoskedastic.".
 COMPUTE temp=t(nms(:,1)).
 OLSMACRO Yols=Y/Xols=G/numrow=n/numcol=ncol(G)/xnames=fnames/ynames=temp.
/*=====Join test for resid variables======.
                     
                    PRINT/title="----".
                    COMPUTE df1=ncol(resid).
                    COMPUTE df2=n-ncol(G).
                    COMPUTE Fstat=((Rsq-Rsqo)/(1-Rsqo))*((n-k)/ncol(resid)).
                    COMPUTE PF=1-fcdf(Fstat,df1,df2).
                    COMPUTE rlabel={"df1","df2","F","Sig."}.
                    COMPUTE out=({df1,df2,Fstat,PF}).
                    PRINT out/title="Joint F test for the signIFicance of the 'resid' variables"
                        /cnames=rlabel
                         /format=f9.3. 
                    PRINT/title="H0: Instrumented variables are exogeneous (OLS is efficient)"+
                   " (i.e., estimates of all resids are equal to zero).".
                   PRINT/title="H1: Instrumented variables are ENDogeneous.".
                    PRINT/title="-----------------------------------------".
END IF.
DO IF (!overz)=1.
/*--- J-Statistics.           
            COMPUTE bj=(inv(sscp(Z)))*t(Z)*eiv.           
            COMPUTE ej=eiv-Z*bj.
            COMPUTE ebar=ones*t(csum(ej))/n.        
            COMPUTE sser=csum(sscp(ej)).
            COMPUTE ssreg=csum(sscp(ebar-Z*bj)).
            COMPUTE Rsq=ssreg/(sser+ssreg).
            COMPUTE Fj=(Rsq/(1-Rsq))*((n-k)/(n-1)).           
            COMPUTE Jstat=n*Fj.
            DO IF (ncol(Z)=ncol(X)).
                  COMPUTE dfj=0.01. 
           else IF (ncol(Z) > ncol(X)).
                 COMPUTE dfj=ncol(Z)-ncol(X).
            END IF.
            COMPUTE PFj=1-chicdf(Jstat,dfj).
            COMPUTE rlabel={"Value", "Sig."}.
            COMPUTE out=t({Jstat,PFj}).
            PRINT out/title="OveridentIFying restrictions test (The J-statistic)"
                 /rnames=rlabel/
                 /format=f9.3.
             PRINT /title="H0: All instruments are exogenous.".
             PRINT/title="The test assumes a large sample size, strong instruments, and "+
    "homoskedastic errors.".
     PRINT/title="-----------------------------------------".
*=== END of J-statistics===.
END IF.
Else IF (ncol(Z) < ncol(X)).
PRINT/ title="==================================================".
PRINT/title="Warning: Number of instruments should be at least equal to the number of ENDogenous "+
    "variables".
prin/title="IdentIFication condition is not satisfied.".
PRINT/title="Thus, Two-stage least square procedure is terminated.".
END IF.
*Weak instrument test.
/*============WEAK INSTRUMENT TEST=============.
COMPUTE ENDov=MENDO(:,2:ncol(MENDO)).
COMPUTE ExZ=MEXO(:,2:ncol(MEXO)).
DO IF (!wktest=1).
PRINT /title="Weak instrument tests".
/*Call Cragg .
CraggD Y=ENDov/X=MINST/Z=ExZ.
COMPUTE sizec={0.10,0.15,0.20,0.25}.
COMPUTE textc={"Size","Values"}.
COMPUTE sizec={0.10,0.15,0.20,0.25}.
COMPUTE textc={"Size","Values"}.
COMPUTE whichrow=ncol(ExZ).
/*Weak Instrument Test Based on TSLS Size.
 Do IF ncol(ENDOV)=1.
/*for n=2.
COMPUTE CVS={16.38,8.96,6.66,5.53;19.93,11.59,8.75,7.25;22.30,12.83,9.54,7.80;+
                       24.58,13.96,10.26,8.31;26.87,15.09,10.98,8.84;29.18,16.23,11.72,9.38;+
                       31.50,17.38,12.48,9.93;33.84,18.54,13.24,10.50;36.19,19.71,14.01,11.07;+
                       38.54,20.88,14.78,11.65;40.90,22.06,15.56,12.23;43.27,23.24,16.35,12.82;+
                       45.64,24.42,17.14,13.41;48.01,25.61,17.93,14.00;50.39,26.80,18.72,14.60;+
                        52.77,27.99,19.51,15.19;55.15,29.19,20.31,15.79;57.53,30.38,21.10,16.39}.
COMPUTE whrow={sizec;CVS(whichrow,:)}.
PRINT t(whrow)/title="Stock-Yogo critical values for the weak instrument test based on TSLS size" +
    " at 5% signIFicance level:"/cnames=textc/format=f9.2.
ELSE IF ncol(ENDOV)=2 AND ncol(ExZ)>1.
COMPUTE CVS={0,0,0,0;7.03,4.58,3.95,3.63;13.43,8.18,6.40,5.45;16.87,9.93,7.54,6.28; +
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
COMPUTE whrow={sizec;CVS(whichrow,:)}.
PRINT t(whrow)/title="Stock-Yogo critical values for the weak instrument test based on TSLS size" +
    " at 5% signIFicance level:"/cnames=textc/format=f9.2.
END IF.
/*------.
/*Weak Instrument Test Based on TSLS Bias.
/*for n=1, k2 at least =3.
COMPUTE whichrow=whichrow-2.
Do IF ncol(ENDOV)=1 AND ncol(ExZ)>2.
COMPUTE CVBS={13.91,9.08,6.46,5.39;16.85,10.27,6.71,5.34;18.37,10.83,6.77,5.25;19.28,11.12,6.76,
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
COMPUTE whrow={sizec;CVBS(whichrow,:)}.
PRINT t(whrow)/title="Stock-Yogo critical values for the weak instrument test based on TSLS "+
    "relative bias" +
    " at 5% signIFicance level:"/cnames=textc/format=f9.2.
ELSE IF ncol(ENDOV)=1 AND ncol(ExZ)>3.
COMPUTE CVBS={11.04,7.56,5.57,4.73;13.97,8.78,5.91,4.79;15.72,9.48,6.08,4.78;16.88,9.92,6.16,4.76;+
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
COMPUTE whrow={sizec;CVBS(whichrow,:)}.
PRINT t(whrow)/title="Stock-Yogo critical values for the weak instrument test based on TSLS "+
    "relative bias" +
    " at 5% signIFicance level:"/cnames=textc/format=f9.2.
ELSE IF ncol(ENDOV)=1 AND ncol(ExZ)>4.
COMPUTE CVBS={9.53,6.61,4.99,4.30;12.20,7.77,5.35,4.40;13.95,8.50,5.56,4.44;15.18,9.01,5.69,4.46;+
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
COMPUTE whrow={sizec;CVBS(whichrow,:)}.
PRINT t(whrow)/title="Stock-Yogo critical values for the weak instrument test based on TSLS "+
    "relative bias" +
    " at 5% signIFicance level:"/cnames=textc/format=f9.2.
ELSE IF ncol(ExZ)<3.
PRINT /title="Stock-Yogo critical values for relative bias are not available for models with "+
    "less than 3 instruments.".
END IF.
DO IF (ncol(ENDOV)=1).
PRINT/ title="Note: For one ENDogeneous variable, Cragg-Donald F-statistic is the"+
                   " F-value of first stage regression. ".
PRINT/ title="For one ENDogeneous regressor, instruments are weak" +
                    " IF the F-statistics is less than 10 (Staiger and Stock, 1997).".
END IF.

PRINT/title="Note: Tables for the critical values (Stock & Yogo, 2002) are reproduced by "+
    "permission--"+
    "email communication (Stock, 16-11-2019).".
PRINT /title="==============================".
END IF.
/*=========END OF WEAK INSTRUMENT TEST.
RELEASE MEXO.
RELEASE MENDO.
RELEASE X.
RELEASE Z.
END MATRIX.
!ENDDEFINE.
TSLS dv = lwage
/iv = educ 
/ivm = motheduc fatheduc
/olsoutp= 0 
/reshaus= 0 
/overz = 0
/bptest= 0
/wktest=0
/robse= 99.
SET WORKSPACE=24576.