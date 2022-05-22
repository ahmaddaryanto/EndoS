# EndoS
# Daryanto, A. (2020). EndoS: An SPSS macro to assess endogeneity. The Quantitative Methods for Psychology, 16(1), 56-70. https://doi.org/10.20982/tqmp.16.1.p056

Never version 2022:
* In the original version, it used 1.96 as the z-critical value of the 95% CI.
It is replaced with 97.5% percentile of the t distribution i.e., IDF.T(0.975, df) in the current version. For a large sample these values are identical. 
* Missing cases omitted
