##########################################
## R-code to conduct LOOCV with the model-data fusion
##########################################
## -- Michele Nguyen; Maricar Rabonza (2022)
##########################################

#-- We assess the performance of the unweighted kriging vs weighted kriging using leave-one-out cross-validation (LOOCV). Cross-validation is a popular method to compare the test performance between models. The total number of points used for testing is maximised via iterative sampling and model fitting. LOOCV is a type of cross-validation where just a single observation is held out for testing in every iteration; thus estimating the out-of-sample estimation performance (Shao, 1993). Detailed steps for the LOOCV procedure are provided in the Supplementary Information, Section S6. 

#-- We can also check if the kriging interpolation overfits the training points. To confirm this, we compare the LOOCV errors to pre-kriging training errors calculated using the same performance metrics (root mean square error, chi-square, and mean square log error) from the process-based model fit. Since in LOOCV we are conducting out-of-sample estimation, LOOCV errors which are higher than the training errors from the Tephra2 fit would suggest that kriging overfits the data.

###############
#-- INPUTS: For this code, 2 components are needed:
# (1) .CSV file containing the observed values and predicted values at each sampled site. Coordinates and the dataset group should also be provided (e.g. Dataset 1 or Dataset 2)
# (2) .OUT file which contains the modelled values at each grid coordinate over the area of interest.

#-- We provide sample inputs with this code, namely:
# "sample-results-MSL-weighted.csv"
# "tephra2-modelled-output.out"

###############
#-- OUTPUT: Table of LOOCV errors (.CSV file)


####--------------------------------------
# Set libraries
####--------------------------------------
rm(list = ls())
library(raster) # handling rasters
library(ggplot2) # the ultimate plotter
library(sp) # raster plotting
library(gstat) # helps with geostatistics and kriging
library(gridExtra) # grid.arrange
library(lemon) # grid_arrange_shared_legend
library(viridis) # for colour-blind palette
library(fields) # For image.plot legends
library(TMB)
library(stats)
library(Matrix)
library(geoR)
library(automap)
library(knitr)
library(RColorBrewer)
library(geoR)
library(rgdal)
library(rgeos)

####--------------------------------------
# Load results from inversion and the forward model (.CSV files)
####--------------------------------------

# (1) Predictions and Observed values at sampled sites
res_data <- read.csv("sample-results-MSL-weighted.csv")
# (2) Forward model gridded output in table form
forward_model <- "tephra2-modelled-output.out" 

####--------------------------------------
# User-defined inputs
####--------------------------------------

# Choose cost function
loss_metric <- "MSL" # Options: MSE, Chi, MSL
# Choose variogram model
vgm_model <- "matern" # exponential, matern, gaussian
# Define map projection to work with
sr_utm49 <- "+proj=utm +zone=49 ellps=WGS84 +init=epsg:23889"
# Factor to reduce resolution of maps
fac_res <- 2 # 2 for a 124 x 125 grid; 5 for 50x50 grid.


####--------------------------------------
# Calculate residuals according to the cost function 
####--------------------------------------

if(loss_metric == "MSE"){
    res_data$Residual <- res_data$Prediction - res_data$Observation}

if(loss_metric == "Chi"){
    res_data$Residual <- res_data$Residual/sqrt(res_data$Observation)}

if(loss_metric == "MSL"){
    res_data$Residual <- log10((res_data$Prediction+1)/(res_data$Observation+1))}


####--------------------------------------
####--------------------------------------
# LEAVE-ONE-OUT CROSS-VALIDATION STARTS HERE
## Un-weighted and weighted model data-fusion
####--------------------------------------

test_id <- which(res_data$Dataset == "Dataset_1") # Test on Dataset 1 only.

for (i in test_id){
    
    training_data <- res_data[-i, ] # Transformed scale for kriging.
    test_data <- res_data[i, ]
    
    #----------
    # Estimate spatial correlation from Dataset 1 training points.
    #----------
    
    ml_fit <- likfit(coords = training_data[training_data$Dataset == "Dataset_1", c("Easting", "Northing")]/1000, data = training_data$Residual[training_data$Dataset == "Dataset_1"], ini=c(0.5, 0.5), fix.nug = FALSE, cov.model = vgm_model, kappa = 1)
    # kappa parameter only works for "matern" cov.model.
    
    ml_spatial_correlation <- function(dists){
        cor <- cov.spatial(dists, cov.model = ml_fit$cov.model, cov.pars = ml_fit$cov.pars, kappa = ml_fit$kappa)
        cor[dists == 0] <- cor[dists == 0] + ml_fit$nugget 
        cor <- cor/(ml_fit$cov.par[1] + ml_fit$nugget)
        return(cor)
    }
    
    x <- seq(0, 10, by = 0.1)
    plot(x, ml_spatial_correlation(x), type = "l", ylab = "Spatial correlation",
         xlab = "Distance", main = "Estimated spatial correlation structure")
    
    #----------
    # Setting up for prediction
    #----------
    
    ### x1_dist = test location only
    x1_dist <- test_data[, c("Easting", "Northing")]
    
    #----------
    # Kriging
    #----------
    
    kriging_pred <- krige.conv(coords = training_data[, c("Easting", "Northing")]/1000, 
                               data = training_data$Residual, locations = x1_dist/1000,
                               krige = krige.control(type.krige = "sk", 
                                                     obj.model = ml_fit))
    
    kriging_df <- cbind(x1_dist, kriging_pred$predict, kriging_pred$krige.var)
    kriging_df <- as.data.frame(kriging_df); colnames(kriging_df) <- c("Easting", "Northing", "Kriging_Mean", "Kriging_Variance")
    
    if(loss_metric == "MSE"){
        
        # Subtract kriged residual from forward prediction.
        kriging_df$Kriging_Mean <- test_data$Prediction - kriging_df$Kriging_Mean
        
    }
    
    if(loss_metric == "Chi"){
        
        kriging_df$Kriging_Mean <- kriging_df$Kriging_Mean*sqrt(test_data$Prediction)
        kriging_df$Kriging_Variance <- kriging_df$Kriging_Variance*test_data$Prediction
        # Subtract kriged residual from forward prediction.
        kriging_df$Kriging_Mean <- test_data$Prediction - kriging_df$Kriging_Mean  
        
    }
    
    if(loss_metric == "MSL"){
        
        # Mean and variance of log_e(Actual + 1)
        lognormal_mu <- (log10(test_data$Prediction+1) - kriging_df$Kriging_Mean)/log10(exp(1)) 
        lognormal_sig <- kriging_df$Kriging_Variance/((log10(exp(1)))^2)
        
        # Mean and variance of Actual using characteristics of a lognormal rv
        kriging_df$Kriging_Mean <- exp(lognormal_mu + 0.5*lognormal_sig) - 1
        kriging_df$Kriging_Variance <- (exp(lognormal_sig) - 1)*exp(2*lognormal_mu + lognormal_sig)
        
    }
    
    #----------
    # Weighted fusion (Worden et al. method)
    #----------
    
    ngrid <- nrow(x1_dist)
    
    ### mu_im1 = mean MMI-equivalent from TEPHRA2/OpenQuake.
    
    mu_im1 <- rep(0, ngrid)
    
    phi_im1 <- rep(sqrt(ml_fit$cov.par[1] + ml_fit$nugget), ngrid)
    tau_im1 <- rep(0, ngrid)
    
    # Observations
    im2 <- training_data$Residual
    
    # Additional sigma for Dataset 2 observations
    obs_extra_sig <- rep(0, length(training_data$Residual))
    obs_extra_sig[training_data$Dataset == "Dataset_2"] <- sqrt(ml_fit$nugget + ml_fit$sigmasq)
    
    # Observation locations
    x2_dist <- training_data[, c("Easting", "Northing")]
    
    # Observation mu
    mu_im2 <- rep(0, length(training_data$Residual)) 
    # Observation uncertainty
    phi_im2 <- rep(sqrt(ml_fit$cov.par[1] + ml_fit$nugget), nrow(x2_dist))
    tau_im2 <-  rep(0, nrow(x2_dist))
    
    # The raw residuals
    zeta <- im2 - mu_im2
    
    # Do the bias correction
    J <- rep(1, length(zeta))
    
    # The Omega factors apply to the bias as well as the MVN,
    # but here we use phi because we dont have var_delta_B yet...
    
    omega_bias <- phi_im2 / sqrt(phi_im2^2 + obs_extra_sig^2)
    
    omega22_bias <- omega_bias %*% t(omega_bias)
    diag(omega22_bias) <- 1.0
    
    J <- J * omega_bias
    
    # Build the covariance matrix of the residuals and its inverse
    sigma22 <- phi_im2 %*% t(phi_im2) 
    # This gives us sigma_Y[1]^2, sigma_Y[2]^2 and sigma_Y[1]*sigma_Y[2] 
    # in the appropriate matrix positions.
    dist22 <- as.matrix(dist(x2_dist))/1000 # Distance matrix, divide by 1000 to work in km (as per variogram estimates).
    
    cov_x2x2 <- ml_spatial_correlation(dist22)
    # This is the correlation matrix.
    cov_x2x2 <- cov_x2x2 * sigma22 # This is Sigma, the original covariance matrix.
    cov_x2x2 <- cov_x2x2 * omega22_bias # This is Sigma', the modified covariance matrix.
    cov_x2x2_inv <- solve(cov_x2x2) # This is the inverse of Sigma'.
    
    # The event term and its variance (for the observation points)
    var_delta_B_im2 <- 1.0 / ((1.0 / tau_im2^2) + t(J) %*% cov_x2x2_inv %*% J) 
    # Equation (12) but with J not all ones.
    delta_B_im2 <- var_delta_B_im2 * (t(J) %*% (cov_x2x2_inv %*% (omega_bias * zeta))) # Equation (11).
    
    # The total within-event standard deviation (observations)
    sigma_delta_W_im2 <- sqrt(phi_im2^2 + var_delta_B_im2) # Within-event + between-event.
    
    # normalized residuals
    x2 <- (zeta - delta_B_im2)/sigma_delta_W_im2 
    # Working on a heteroskedasticity-invariant scale.
    
    # The event term, its variance and the total within-event standard
    # deviation (for the predictions)
    var_delta_B_im1 <- 1.0 / ((1.0 / tau_im1^2) + t(J) %*% cov_x2x2_inv %*% J) 
    delta_B_im1 <- var_delta_B_im1 * (t(J) %*% (cov_x2x2_inv %*% (omega_bias * zeta)))
    
    sigma_delta_W_im1 <- sqrt(phi_im1^2 + var_delta_B_im1)
    
    # Correlation adjustment for additional uncertainty
    omega_1 <- matrix(rep(1, length(phi_im1)), nrow = 1, ncol = length(phi_im1))
    omega_2 <- t((sigma_delta_W_im2 /
                      sqrt(sigma_delta_W_im2^2 + obs_extra_sig^2)))
    
    # Create the Omega' sub-matrices (Equation (46)-(47))
    omega_prime_12 <- t(omega_1) %*% omega_2
    omega_prime_22 <- t(omega_2) %*% omega_2
    diag(omega_prime_22) <- 1.0
    
    # Adjust the normalized residuals by the omega factors (Equation (48))
    x2_adjusted <- x2 * t(omega_2)
    
    # Distance matrices
    
    # Error: cannot allocate vector of size 28.8 Gb for full prediction grid. 
    
    dist_mega <- as.matrix(dist(rbind(x1_dist, x2_dist)))/1000
    dist11 <- dist_mega[1:nrow(x1_dist), 1:nrow(x1_dist)]
    dist22 <- dist_mega[(nrow(x1_dist)+1):nrow(dist_mega),
                        (nrow(x1_dist)+1):nrow(dist_mega)]
    dist12 <- dist_mega[1:nrow(x1_dist), (nrow(x1_dist)+1):nrow(dist_mega)]
    
    cov_x1x1 <-  ml_spatial_correlation(dist11)
    
    cov_x2x2 <-  ml_spatial_correlation(dist22)
    cov_x2x2 <- cov_x2x2* omega_prime_22
    
    cov_x1x2 <-  ml_spatial_correlation(dist12)
    cov_x1x2 <- cov_x1x2* omega_prime_12
    
    # Inverse of cov_x2x2
    cov_x2x2_inv <- solve(cov_x2x2)
    
    # Solve for the normalized conditional mean and covariance
    
    # Regression coefficient matrix
    rcmatrix <- cov_x1x2 %*% cov_x2x2_inv
    
    # Conditional normalized mean
    mu_x1_x2 <- rcmatrix %*% x2_adjusted # mu[Y[1]] = 0 because normalised. 
    
    # Conditional normalized covariance
    cov_x1x1_x2 = cov_x1x1 - rcmatrix %*% t(cov_x1x2)
    
    # Return to observational units
    
    # de-normalize the conditional mean
    mu_im1_im2 <- mu_im1 + (sigma_delta_W_im1 * mu_x1_x2 + delta_B_im1) 
    # Unnormalise and add mu_im1 because zeta substracted it before
    
    # de-normalize the conditional covariance
    phi_x1x1 <- sigma_delta_W_im1 %*% t(sigma_delta_W_im1)
    cov_im1im1_im2 <- phi_x1x1 * cov_x1x1_x2
    
    temp_diag <- diag(cov_im1im1_im2)
    
    temp_diag[temp_diag < 0] <- 0
    
    sig_im1im1_im2_diag <- sqrt(temp_diag)
    
    # Test performance
    Res_Pred <- mu_im1_im2
    
    worden_df <- kriging_df; colnames(worden_df) <- c("Easting", "Northing", "Worden_Mean", "Worden_Variance")
    
    if(loss_metric == "MSE"){
        
        # Subtract kriged residual from forward prediction.
        worden_df$Worden_Mean <- test_data$Prediction - Res_Pred
        worden_df$Worden_Variance <- sig_im1im1_im2_diag^2
        
    }
    
    if(loss_metric == "Chi"){
        
        worden_df$Worden_Mean <- Res_Pred*sqrt(test_data$Prediction)
        worden_df$Worden_Variance <- (sig_im1im1_im2_diag^2)*test_data$Prediction
        
        # Subtract kriged residual from forward prediction.
        worden_df$Worden_Mean <- test_data$Prediction - worden_df$Worden_Mean 
        
    }
    
    if(loss_metric == "MSL"){
        
        # Mean and variance of log_e(Actual + 1)
        lognormal_mu <- (log10(test_data$Prediction+1) - Res_Pred)/log10(exp(1)) 
        lognormal_sig <- (sig_im1im1_im2_diag^2)/((log10(exp(1)))^2)
        
        # Mean and variance of Actual using characteristics of a lognormal rv
        worden_df$Worden_Mean <- exp(lognormal_mu + 0.5*lognormal_sig) - 1
        worden_df$Worden_Variance <- (exp(lognormal_sig) - 1)*exp(2*lognormal_mu + lognormal_sig)
        
    }
    
    worden_df$Worden_Mean[worden_df$Worden_Mean < 0] <- 0
    
    ########### Pre-kriging, simple-kriging (unweighted fusion) and weighted fusion errors.
    
    test_data$Kriging_Mean <- kriging_df$Kriging_Mean
    test_data$Kriging_Variance <- kriging_df$Kriging_Variance
    test_data$Worden_Mean <- worden_df$Worden_Mean
    test_data$Worden_Variance <- worden_df$Worden_Variance
    
    if(i == 1){
        
        loo_results <- test_data
        
    }else{
        loo_results <- rbind(loo_results, test_data)
    }
    
    
}

loo_results

# =======
# END OF LEAVE-ONE-OUT CROSS-VALIDATION 
####--------------------------------------
####----------------------------------------------------------------------------


####--------------------------------------
# Save test results
####--------------------------------------
write.csv(loo_results, file = paste(loss_metric, "_", vgm_model, "_loo_results.csv", sep = ""), row.names = FALSE)


####--------------------------------------
# Display LOOCV errors
####--------------------------------------

# RMSE
sqrt(mean((loo_results$Observation - loo_results$Prediction)^2))
# Post-kriging
sqrt(mean((loo_results$Observation - loo_results$Kriging_Mean)^2))
# Post-Worden
sqrt(mean((loo_results$Observation - loo_results$Worden_Mean)^2))

# Chi-squared
mean(((loo_results$Observation - loo_results$Prediction)^2)/loo_results$Observation)
# Post-kriging
mean(((loo_results$Observation - loo_results$Kriging_Mean)^2)/loo_results$Observation)
# Post-Worden
mean(((loo_results$Observation - loo_results$Worden_Mean)^2)/loo_results$Observation)

# MSLE
mean(log10((loo_results$Prediction+1)/(loo_results$Observation+1))^2)
# Post-kriging
mean(log10((loo_results$Kriging_Mean+1)/(loo_results$Observation+1))^2)
# Post-Worden
mean(log10((loo_results$Worden_Mean+1)/(loo_results$Observation+1))^2)
