
# About 
This family of functions assists with the cleaning and analysis of spatial data to execute geographic regression discontinuity designs. Of the six functions, three are for cleaning and preparing data for analysis, and the other three are for executing the analysis. Below, I describe all six functions at a high-level. The actual scripts contain more detailed information about each function specifically, which should be used as a reference. Of course, please contact me with any problems that may arise. 

# Cleaning Data

## define_RD

This function takes as inputs the polygons (areas) that define places that are treated (action) and areas that are not treated (control) as well as the study area (typically the union of the control and action areas) as well as the spatial data (typically point data) which includes the units for analysis (individuals, towns, ect). From there, the function defines the distance of each unit from the study border and categegorizes the unit as being under either treatement or control. Further, the function may optionally define a unit's proximity to a border segement point. Segments are points along the border line that seperate treatment and control areas. One may wish to use segments if the border line is very long and if the researcher is interested in only comparing units that are near to eachother. 

## makes_study_border 

This function finds the border between the treatment and control polygons. This function is placed inside of the define_RD function, but is broken out of that function for individual use in case a researcher is intrested in extracting the border line for purpose of data visualization. 

## make_diagnostic_map

This function creates a diagnostic heatmap through kriging. In brief, the function uses raw spatial data to predict values across space, and then plots the predictions along with the border line to visualize values of the outcome variable across space. This is the two-dimensional spatial analog to the standard univariate RD plot which shows a break in the outcome variable (a discontinuity) at the value of the assignment variable which defines treatment. Should there be a change in the outcome of interest, one should be able to see a change in color at the borderline. 

# Analyzing Data

The following functions implement methodological advise from Kelly (2020) and Keele and Titiunik (2015). Please reference those papers for more detail on the background of the problems and proposed solutions that these functions are intended to address. 

## grdrobust 

This function performs a semi-parametric geographic RD while computing spatial heteroskedasticity and autocorrlation consistent (SHAC) standard errors. To construct SHAC errors, the function finds the spatial structure of the residuals to construct an empirical range of spatial autocorrelation. This is a useful innovation, because it allows researchers to define a range of autocorrelation that is not arbitrary guess, which reduces researcher degrees of freedom in selecting a range to perform their analysis. This enables reserchers to depart from the common practice of reporting SHAC errors at various ranges without theoretical guidance.  In doing so, the function also finds the bandwidth for the analysis, allows researchers to optionally include fixed effects for border segments, using distance or latitude-longitude for the running variable, to specify different polynomial degrees of the assignment variable, or to include covariates. Please note the function takes time to run because the process of finding the range of autocorrelation with MLE is computionally expensive.

## noise_test 

This function tests for whether spatial noise outperforms the observed outcome variable in GRD. 

## segment_RD 

This function computes the geographic RD within border segments, and then reports the pooled and individual estimates at border points. 
