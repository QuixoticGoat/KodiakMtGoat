The United States Fish and Wildlife Service (FWS) GitHub project code is provided on an
"as is" basis and the user assumes responsibility for its use. FWS has relinquished control
of the information and no longer has responsibility to protect the integrity, confidentiality, or
availability of the information. Any reference to specific commercial products, processes,
or services by service mark, trademark, manufacturer, or otherwise, does not constitute or
imply their endorsement, recommendation or favoring by FWS. The FWS seal and logo
shall not be used in any manner to imply endorsement of any commercial product or
activity by FWS or the United States Government.

# Kodiak Mountain Goat Resource Selection
An analysis of mountain goat resource selection on Kodiak Island, Alaska. We compared covariate values at GPS collar fixes (used) to those at available locations. We considered individual collared mountain goats as the sampling unit. For the population level analysis, we evaluated the relative probability of use using conditional logistic mixed effects models with a matched case-control design (Fortin et al. 2005, Duchesne et al. 2010). Models included a random intercept for each collared mountain goat to account for the unbalanced design and correlation among individuals and allow for population level inference. We quantify resource selection separately for summer and winter seasons. We defined the summer season as June-September and the winter season as December-April based on the typical presence of snow cover.

For every collared mountain goat, we draw 10 random steps in the study area for every used fix from an empirical distribution of step lengths and turning angles based on all other mountain goats with similarly programmed collars (i.e., separately for ATS and Telonics collars). We applied the sample of steps to each fix to generate 10 matched (case-control) available fixes for each used fix. We then extracted spatial covariates values to the used and available fixes. We withheld a random sample of 20% of the used and available dataset (“testing data”) for model validation (detailed below).

## Covariates
Similar to other large herbivores, mountain goat select habitats that minimize predation pressure but maximize opportunities of high quality forage consumption (Gross et al. 2002, White 2006). With these needs in mind and using information on Kodiak mountain goat diets (Hjeljord 1973), we selected a priori terrain and vegetation covariates that we hypothesized affect mountain goat resource selection (Table X). We derived terrain covariates [slope, aspect and vector terrain ruggedness  (VTR, Sappington et al. 2007)] from one arc-second (30-m) USGS National Elevation data (NED). We defined aspect as a continuous numeric index between 0 (north) and 1 (south). We obtained habitat covariates using a landscape cover classification of the Kodiak Archipelago derived from 30-m LandSat ETM+ imagery classification (Fleming and Spencer 2007). We selected six habitat classes from the coarsest hierarchical classification level that we hypothesized were most relevant to mountain goats: meadow (including forbs and graminoid dominant grassland and meadows), tundra (alpine tundra and lowland heath), shrub (dominated by alder, salmonberry and willow), forest (dominated by Sitka spruce, birch and cottonwoods), rock (solid and fragmented), and water (fresh). Based on field observations, we expected that females would select steeper slopes in more rugged terrain in the summer than males because of their unique need to provide proximate escape terrain for kids.   Therefore, we included sex as an individual-level covariate. 

We followed a multi-tiered approach to model selection (Lowrey et al. 2017). In the first tier, we fit univariate models of terrain covariates (slope and vector ruggedness measure) derived from USGS NEDs. We used AICc for model selection. In the second tier, we began with the top candidate terrain model as a base model and then evaluated all combinations of habitat type covariates, which included Forest, Shrubs, Tundra/Heath, Meadow, Water, Snow and Rock. 

To assess the predictive performance of the most supported models, we evaluated the correlations between the frequencies of occurrence of the testing data and their relative RSF scores using Spearman’s rank correlation coefficients (Boyce et al. 2002). High correlation indicates a well performing model.

## **Analysis** Folder Contents

*ImportFormat.R*  
Imports GPS collar data from ATS and Telonics collars. Reformats the data, removes erronious locations, and merges it into a single dataframe. Creates a spatial dataframe and plots it as a visual check for errors.

*FormatSpatialCov.R*  
Imports and formats spatial covariate shapefiles and rasters. Scales covariates when appropriate. Creates binomial rasters for each habitat class in the land cover classification (LCC). Creates and saves a raster stack of the formatted habitat covariate rasters. Imports and plots the raster stack. Saves plot as a pdf.

*ExtractCovVelox.R*  
Defines "available" fixes for standard RSF modeling by creating a 99% kernel density spatial polygon, using all fixes. Generates a random sample of fixes within the polygon and combines these data to the "used" GPS location dataframe, called dfRSF. Imports a raster stack of habitat covariates. Extract the covariate data from the raster stack to the used/available data. Export the result as a shapefile.

*ConditionalRSFmodel.R*  
Step 1: Examines resource selection based on mixed effects conditional logistic regression (CollarID as a random effect). To do this, it first creates a two trajectory objects using the "used" GPS locations (one for ATS collars, one for Telonics collars) that include the length of each move, the time interval between successive relocations and other values. Load a study area boundary (mask). Draws 10 random "available" steps within the study area bounday for each "used" step from a distribution of used fix step lengths and turning angles. Merges the ATS and Telonics datasets back together. Calculates the available fix coordinates from the random step lengths and turning angles. Saves the results.  

Step 2: Extracts spatial covariates raster values to the used and available fixes. Loads the used/available dataframe and the raster stack with the covariates. Converts the dataframe to a spatial dataframe. Extracts the covariate values to the used/available spatial data frame and saves the results as a dataframe.  

Step 3: Splits the used/available data frame into testing and training dataframes. 

Step 4: Run the conditional mixed effects model selection for summer and winter. Loads the used/available "training" dataframe from Step 3. Subsets summer data. Evaluates spatial grain of the VRM and slope covariates with univariate models via delta AICc. Saves the list of models and AICc table. Using the top model from the previous step as the base model, evaluates a candidate list of terrian models via delta AICc. Saves list of models and AICc table. Repeats process for the habitat covariates. Saves list of models and AICc tables and saves the "best" summer model. Subsets winter data and repeats the steps above to determine the "best" winter model.  

Step 5: Evaluate fit of the best summer and winter models using the testing data.

*RSFmodel.R*  
Examines resource selection based on mixed effects logistic regression (CollarID as a random effect). Checks for covariance among the habitat covariates. Saves the resulting table as a .csv. Splits the used/available data frame into testing and training dataframes. Creates a list of terrian models, based on topography. Checks models for overdispersion. Runs AIC model selection to determine the best model. Uses this as the base model for a second set of candidate models based on habitat covariates. Runs AIC model selection to determine the best model. Checks models for overdispersion. Creates a table of parameter estimates and CIs from the top model. Adds predicted values (based on the top model) to the used/available data frame. Evaluates the goodness of fit of the best model with Hosmer-Lemeshow Goodness of Fit (GOF) test. 

*RSFmodel_summer.R*  
Examines summer (June-Sept) resource selection based on mixed effects logistic regression (CollarID as a random effect). Checks for covariance among the habitat covariates. Saves the resulting table as a .csv. Splits the used/available data frame into testing and training dataframes. Creates a list of terrian models, based on topography. Checks models for overdispersion. Runs AIC model selection to determine the best model. Uses this as the base model for a second set of candidate models based on habitat covariates. Runs AIC model selection to determine the best model. Checks models for overdispersion. Creates a table of parameter estimates and CIs from the top model. Adds predicted values (based on the top model) to the used/available data frame. Evaluates the goodness of fit of the best model with Hosmer-Lemeshow Goodness of Fit (GOF) test. 

*RSFmodel_winter.R*  
Examines winter (Dec-Apr) resource selection based on mixed effects logistic regression (CollarID as a random effect). Checks for covariance among the habitat covariates. Saves the resulting table as a .csv. Splits the used/available data frame into testing and training dataframes. Creates a list of terrian models, based on topography. Checks models for overdispersion. Runs AIC model selection to determine the best model. Uses this as the base model for a second set of candidate models based on habitat covariates. Runs AIC model selection to determine the best model. Checks models for overdispersion. Creates a table of parameter estimates and CIs from the top model. Adds predicted values (based on the top model) to the used/available data frame. Evaluates the goodness of fit of the best model with Hosmer-Lemeshow Goodness of Fit (GOF) test. 

*RSFmaps.R*  
Creates an RSF surface raster using the exponential equation from the top model (as determined via model selection in RSFmodel.R) and saves it. Bins the RSF surface into quantiles and saves it. Plots the raster surface.

*PlotRSFSims.R*  
Simulates (bootstrapped) fixed effects posterior distributions of the top model (as determined via model selection in RSFmodel.R). Plots it. Predicts values based on the top model.
