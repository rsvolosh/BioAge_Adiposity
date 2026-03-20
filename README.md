
# Greater than the sum of its parts: combining epigenetic clocks to characterize the association of biological age acceleration and adiposity in young Filipino adults

<!-- badges: start -->
<!-- badges: end -->

## Abstract

### Background

Diverse epigenetic clocks are known to capture health risks associated with increased adiposity, but their estimates have never been combined to represent a holistic estimate of biological age acceleration (BAA). There is also a gap in research using epigenetic clocks to study adiposity in lower-middle income Asian countries.

### Methods and Findings

Data from 1,745 participants (21.7±0.3 years old, 45% female) of the Cebu (Philippines) Longitudinal Health and Nutrition Survey were analyzed. BAA was calculated using PCHorvath 2, PCHannum, PCPhenoAge, PCGrimAge, PCDNAmTL, and DunedinPACE. After ascertaining suitability for factor analysis (Kaiser-Meyer-Olkin 0.81), factor analysis was used to create PCFactorAge. Analogously, FactorAge was created using Horvath, Hannum, PhenoAge, GrimAge, DNAmTL, and DunedinPACE. BMI, waist circumference (WC), and waist-to-height ratio (WHtR) were used to represent adiposity. Linear regression was used to test the association of each adiposity measure with each BAA measure. 
BMI, WC, and WHtR were positively associated with both BAA combinations: 5 kg/m2 higher BMI corresponded to 0.097 (p=0.015) standard deviation (SD) increase in FactorAge and 0.099 (p=0.004) SD increase in PCFactorAge; 10 cm increase in WC—with 0.091 (p=0.005) SD increase in FactorAge and 0.094 (p<0.001) SD increase in PCFactorAge; 0.1 increase in WHtR—with 0.164 (p=0.001) SD increase in FactorAge and 0.163 (p<0.001) SD increase in PCFactorAge. Additionally, WHtR was associated with meaningful increases in PhenoAge, PCPhenoAge, PCHorvath 2, PCHannum, PCGrimAge, and DunedinPACE. WC was positively associated with PCHorvath 2, PCHannum, PCPhenoAge, and DunedinPACE. BMI was positively associated with PCHannum, PCPhenoAge, and DunedinPACE.

### Conclusions

Our study presents a novel approach to creating a BAA estimate using multiple epigenetic clocks and shows that adiposity measures predict this factor in a young Filipino cohort.


## Implementation Guide

Code in the “imputation” folder must be run first.

Then, the main data set is to be merged using “prep_data.R”.

The main analyses are in “analysis-tables-figures.Rmd”. The knitted file
is “analysis-tables-figures.html”.

“methods_code.R” shows how the analytical methods (specifically,
applying factor analysis to multiple epigenetic clocks to find the
latent biological age acceleration variable) described in this paper can
be applied to other data sets.
