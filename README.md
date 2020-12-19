# EndoS
Daryanto, A. (2020). EndoS: An SPSS macro to assess endogeneity. The Quantitative Methods for Psychology, 16(1), 56-70. https://doi.org/10.20982/tqmp.16.1.p056
What's new in this version:
the confidence interval: 
compute tcrit = IDF.T(0.975, dffols).     
compute LBols=bols-tcrit*sbols.
compute UBols=bols+tcrit*sbols.
in the old version: it used 1.96 for tcrit.
