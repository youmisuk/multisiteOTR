# Designing Optimal, Data-Driven Policies from Multisite Randomized Trials

Youmi Suk<sup>1</sup> and Chan Park<sup>2</sup>

<sup>1</sup> Department of Human Development, Teachers College Columbia University  
<sup>2</sup> Department of Statistics and Data Science, The Wharton School, University of Pennsylvania


## Overview

Optimal treatment regimes (OTRs) have been widely employed in computer science and personalized medicine to provide data-driven, optimal recommendations to individuals. However, previous research on OTRs has primarily focused on settings that are independent and identically distributed, with little attention given to the unique characteristics of educational settings, where students are nested within schools and there are hierarchical dependencies. The main goal of this study is to design OTRs from multisite randomized trials, a commonly used experimental design in education and psychology to evaluate educational programs. We investigate modifications to popular OTR methods, specifically Q-learning and weighting methods, in order to improve their performance in multisite randomized trials. A total of 12 modifications, 6 for Q-learning and 6 for weighting, are proposed by utilizing different multilevel models, moderators, and augmentations. Simulation studies reveal that all Q-learning modifications improved performance in multisite randomized trials, with modifications that incorporate random treatment effects showing the most promise in handling cluster-level moderators. Among weighting methods, the modification that incorporates cluster dummies into moderator variables and augmentation terms performed best across simulation conditions. The proposed modifications are demonstrated through an application to estimate an OTR of conditional cash transfer programs using a multisite randomized trial in Colombia to maximize educational attainment.

For more details of our proposed methods, see [our paper](). 
Here, we provide `R` codes to reproduce our simulation study and replicate our data analysis using data about conditional cash transfer (CCT) programs. 

## Simulation Study

* `DataGeneratingModels.R`  

   This `R` file includes data generating codes for data from a multisite randomized trial with a cluster-level unmeasured covariate.

* `Qlearn.R`   

   This `R` file includes a function named `Qlearn` to implement our proposed modifcations for Q-learning as well as the baseline Q-learning method.

* `SimulationCodes.R`
 
   This `R` file includes simulation codes with our proposed modifcations for Q-learning and weighting methods where the parameter `beta1` represents the cofficient of a cross-level interaction effect between treatment status and a cluster-level unmeasured covariate. For more information on simulation condtions, please refer to [our paper](XXXXXXXXXXXXXXXXX).



## CCT Data Study

* Data on conditional cash transfer (CCT) programs
  
  For our empirical analysis, we used data collected by [Barrera-Osorio et al. (2011)](https://doi.org/10.1257/app.3.2.167). The data can be downloaded from ICPSR by clicking [here](https://doi.org/10.3886/E113783V1). For more information on the data, please refer to the detiled report by [Barrera-Osorio et al. (2011)](https://doi.org/10.1257/app.3.2.167) and the codebook provided [here](https://doi.org/10.3886/E113783V1). 

* `DataAnalysisCodes.R` 
 
   This `R` file can be used to replicate our data analysis.

Please note that these supplemental materials are provided for the purpose of reproducibility and should be used in accordance with academic ethical guidelines. Any reference to these materials in your work should properly cite the original sources.
