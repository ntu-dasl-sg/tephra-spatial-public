rm(list = ls())

# Set libraries
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

# Factor to reduce resolution of maps:
fac_res <- 2 # 2 for a 124 x 125 grid; 5 for 50x50 grid.

# Define map projection to work with
sr_utm49 <- "+proj=utm +zone=49 ellps=WGS84 +init=epsg:23889"

data.path <- "/Users/maricar/Github_repositories/tephra_v3.0/Runs/20220801_widerrange/Output CSV/"
#"D:/Documents/DASL Proj Tephra Spatial/Data/"
graphics.path <- "/Users/maricar/Github_repositories/tephra_v3.0/Runs/20220801_widerrange/Graphics/"
#"D:/Documents/DASL Proj Tephra Spatial/Graphics/"

# Choose loss metric and variogram model:
loss_metric <- "MSL" # MSE, Chi, MSL
vgm_model <- "matern" # exponential, matern, gaussian

##########-----
#--- SET FILENAMES:

# (1) Residual data with dataset information
res_pts <- "s03_msl_wt_ds2half_resid.csv" # ============ Change here.

# Choose which residual data to use (here: Residual = Prediction - Load)
res_data <- read.csv(paste(data.path, res_pts, sep = ""), head=TRUE, sep=",")
full_data <- res_data # Save a copy in the original (untransformed) scale.

if(loss_metric == "Chi"){
    
    res_data$Residual <- res_data$Residual/sqrt(res_data$Load)
    
}

if(loss_metric == "MSL"){
    
    res_data$Residual <- log10((res_data$Prediction+1)/(res_data$Load+1))
    
}

# (2) Forward model gridded output in table form

forward_model <- "msl_uw.out" # ============ Change here.

# ====== LEAVE-ONE-OUT CROSS-VALIDATION STARTS HERE =======

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
    # Worden et al. method
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
    
    ########### Pre-kriging, vanilla-kriging and Worden et al. errors.
    
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

# Save test results
write.csv(loo_results, file = paste(data.path, loss_metric, "_", vgm_model, "_loo_results.csv", sep = ""), row.names = FALSE)

# To read in previous results:
# loo_results <- read.csv(file = paste(data.path, loss_metric, "_", vgm_model, "_loo_results.csv", sep = ""), row.names = FALSE)

# RMSE
sqrt(mean((loo_results$Load - loo_results$Prediction)^2))
# Post-kriging
sqrt(mean((loo_results$Load - loo_results$Kriging_Mean)^2))
# Post-Worden
sqrt(mean((loo_results$Load - loo_results$Worden_Mean)^2))

# Chi-squared
mean(((loo_results$Load - loo_results$Prediction)^2)/loo_results$Load)
# Post-kriging
mean(((loo_results$Load - loo_results$Kriging_Mean)^2)/loo_results$Load)
# Post-Worden
mean(((loo_results$Load - loo_results$Worden_Mean)^2)/loo_results$Load)

# MSLE
mean(log10((loo_results$Prediction+1)/(loo_results$Load+1))^2)
# Post-kriging
mean(log10((loo_results$Kriging_Mean+1)/(loo_results$Load+1))^2)
# Post-Worden
mean(log10((loo_results$Worden_Mean+1)/(loo_results$Load+1))^2)


# ===================== Conduct spatial interpolation on all data.

#----------
# Load forward model prediction for prediction grid
#----------

# Read file
pred.chi2 <- read.table(paste(data.path, forward_model, sep = ""))  
# Rename columns - Third column has nothing (just zeroes) so remove that
names(pred.chi2) <- (c("easting", "northing", "na", "load"))
# Remove unnecessary rows - which is 3rd column
keep <- c("easting", "northing", "load") #keep these rows
pred.chi2 <- pred.chi2[keep]
# Convert prediction data frame to raster file
r.pred <- rasterFromXYZ(pred.chi2[,c(1,2,3)]) # select columns for x,y,z
# Define projection of raster
proj4string(r.pred) <- CRS(sr_utm49)  
# Crop raster to reduce file size
e <- extent(590000, 660000, 9100000, 9162000)
#e <- extent(535000, 670000, 9100000, 9180000) # For isopach.
# e <- extent(min(res_data$Easting), max(res_data$Easting), min(res_data$Northing), max(res_data$Northing)) 
r.pred <- crop(r.pred, e)

# Reduce resolution of forward solutions if necessary:
r.pred <- aggregate(r.pred, fact=fac_res)

#----------
# Estimate spatial correlation from Dataset 1 training points.
#----------

emp_vgm <- variog(coords = res_data[res_data$Dataset == "Dataset_1", c("Easting", "Northing")]/1000, data = res_data$Residual[res_data$Dataset == "Dataset_1"], max.dist = 50)

ml_fit <- likfit(coords = res_data[res_data$Dataset == "Dataset_1", c("Easting", "Northing")]/1000, data = res_data$Residual[res_data$Dataset == "Dataset_1"],  
                 #fix.nug = TRUE, nugget = 0.0104928*5, ini=c(0.5, 0.5), cov.model = vgm_model, kappa = 1)
                 fix.nug = FALSE, cov.model = vgm_model, kappa = 1, ini=c(0.5, 0.5)) # nugget = 0.0104928.
# kappa parameter only works for "matern" cov.model.

png(file = paste(graphics.path, loss_metric, "_", vgm_model, "_ml_fit.png", sep = ""),
    height = 1200,
    width = 1500,
    res = 300)
plot(emp_vgm, ylab = expression(gamma), xlab = "Distance (km)", main = "MSLE")
lines(ml_fit)
dev.off()

# check for directionality

vario.0 <- variog(coords = res_data[res_data$Dataset == "Dataset_1", c("Easting", "Northing")]/1000, data = res_data$Residual[res_data$Dataset == "Dataset_1"], max.dist = 50, dir=0, tol=pi/8)
vario.45 <- variog(coords = res_data[res_data$Dataset == "Dataset_1", c("Easting", "Northing")]/1000, data = res_data$Residual[res_data$Dataset == "Dataset_1"], max.dist = 50, dir=pi/4, tol=pi/8)
vario.90 <- variog(coords = res_data[res_data$Dataset == "Dataset_1", c("Easting", "Northing")]/1000, data = res_data$Residual[res_data$Dataset == "Dataset_1"], max.dist = 50, dir=pi/2, tol=pi/8)
vario.135 <- variog(coords = res_data[res_data$Dataset == "Dataset_1", c("Easting", "Northing")]/1000, data = res_data$Residual[res_data$Dataset == "Dataset_1"], max.dist = 50, dir=3*pi/4, tol=pi/8)

png(file = paste(graphics.path, loss_metric, "_", vgm_model, "_directional.png", sep = ""),
    height = 2400,
    width = 3000,
    res = 300)
par(mfrow = c(2, 2))
plot(vario.0, ylab = expression(gamma), xlab = "Distance (km)", main = expression(0 * degree))
lines(ml_fit)
plot(vario.45, ylab = expression(gamma), xlab = "Distance (km)", main = expression(45 * degree))
lines(ml_fit)
plot(vario.90, ylab = expression(gamma), xlab = "Distance (km)", main = expression(90 * degree))
lines(ml_fit)
plot(vario.135, ylab = expression(gamma), xlab = "Distance (km)", main = expression(135 * degree))
lines(ml_fit)
dev.off()

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

### x1_dist = observation locations/coordinates from TEPHRA2/OpenQuake.

#x1_dist <- res_data[, c("Easting", "Northing")] # If we are computing residual error.
x1_dist <- coordinates(r.pred)
colnames(x1_dist) <- c("Easting", "Northing")

#----------
# Kriging
#----------

kriging_pred <- krige.conv(coords = res_data[, c("Easting", "Northing")]/1000, 
                           data = res_data$Residual, locations = x1_dist/1000,
                           krige = krige.control(type.krige = "sk", 
                                                 obj.model = ml_fit))

kriging_df <- cbind(x1_dist, kriging_pred$predict, kriging_pred$krige.var)
kriging_df <- as.data.frame(kriging_df); colnames(kriging_df) <- c("Easting", "Northing", "Kriging_Mean", "Kriging_Variance")

if(loss_metric == "MSE"){
    
    # Subtract kriged residual from forward prediction.
    kriging_df$Kriging_Mean <- as.vector(r.pred) - kriging_df$Kriging_Mean
    
}

if(loss_metric == "Chi"){
    
    kriging_df$Kriging_Mean <- kriging_df$Kriging_Mean*sqrt(as.vector(r.pred))
    kriging_df$Kriging_Variance <- kriging_df$Kriging_Variance*as.vector(r.pred)
    
    # Subtract kriged residual from forward prediction.
    kriging_df$Kriging_Mean <- as.vector(r.pred) - kriging_df$Kriging_Mean  
    
}

if(loss_metric == "MSL"){
    
    # Mean and variance of log_e(Actual + 1)
    lognormal_mu <- (log10(as.vector(r.pred)+1) - kriging_df$Kriging_Mean)/log10(exp(1)) 
    lognormal_sig <- kriging_df$Kriging_Variance/((log10(exp(1)))^2)
    
    # Mean and variance of Actual using characteristics of a lognormal rv
    kriging_df$Kriging_Mean <- exp(lognormal_mu + 0.5*lognormal_sig) - 1
    kriging_df$Kriging_Variance <- (exp(lognormal_sig) - 1)*exp(2*lognormal_mu + lognormal_sig)
    
}

kriging_df$Kriging_Mean[kriging_df$Kriging_Mean<0] <- 0

plot(kriging_df$Kriging_Mean, as.vector(r.pred))
abline(a = 0, b = 1)

#----------
# Worden et al. method
#----------

ngrid <- nrow(x1_dist)

### mu_im1 = mean MMI-equivalent from TEPHRA2/OpenQuake.

mu_im1 <- rep(0, ngrid)

phi_im1 <- rep(sqrt(ml_fit$cov.par[1] + ml_fit$nugget), ngrid)
tau_im1 <- rep(0, ngrid)

# Observations
im2 <- res_data$Residual

# Additional sigma for Dataset 2 observations
obs_extra_sig <- rep(0, length(res_data$Residual))
obs_extra_sig[res_data$Dataset == "Dataset_2"] <- sqrt(ml_fit$nugget + ml_fit$sigmasq)

# Observation locations
x2_dist <- res_data[, c("Easting", "Northing")]

# Observation mu
mu_im2 <- rep(0, length(res_data$Residual)) 
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
    worden_df$Worden_Mean <- as.vector(r.pred) - Res_Pred
    worden_df$Worden_Variance <- sig_im1im1_im2_diag^2
    
}

if(loss_metric == "Chi"){
    
    worden_df$Worden_Mean <- Res_Pred*sqrt(as.vector(r.pred))
    worden_df$Worden_Variance <- (sig_im1im1_im2_diag^2)*as.vector(r.pred)
    
    # Subtract kriged residual from forward prediction.
    worden_df$Worden_Mean <- as.vector(r.pred) - worden_df$Worden_Mean 
    
}

if(loss_metric == "MSL"){
    
    # Mean and variance of log_e(Actual + 1)
    lognormal_mu <- (log10(as.vector(r.pred)+1) - Res_Pred)/log10(exp(1)) 
    lognormal_sig <- (sig_im1im1_im2_diag^2)/((log10(exp(1)))^2)
    
    # Mean and variance of Actual using characteristics of a lognormal rv
    worden_df$Worden_Mean <- exp(lognormal_mu + 0.5*lognormal_sig) - 1
    worden_df$Worden_Variance <- (exp(lognormal_sig) - 1)*exp(2*lognormal_mu + lognormal_sig)
    
}

worden_df$Worden_Mean[worden_df$Worden_Mean < 0] <- 0

#----------
# Residual error and variogram after kriging/Worden
#----------

# Residual error before and after kriging (run with x1_dist for points only, above code with r.pred will give errors)

# lognormal_mu <- (log10(as.vector(res_data$Prediction)+1) - kriging_df$Kriging_Mean)/log10(exp(1)) 
# lognormal_sig <- kriging_df$Kriging_Variance/((log10(exp(1)))^2)
# 
# kriging_df$Kriging_Mean <- exp(lognormal_mu + 0.5*lognormal_sig) - 1
# kriging_df$Kriging_Variance <- (exp(lognormal_sig) - 1)*exp(2*lognormal_mu + lognormal_sig)
# 
# mean((kriging_df$Kriging_Mean-res_data$Load)^2)
# mean((res_data$Prediction-res_data$Load)^2)
# 
# lognormal_mu <- (log10(as.vector(res_data$Prediction)+1) - Res_Pred)/log10(exp(1)) 
# lognormal_sig <- (sig_im1im1_im2_diag^2)/((log10(exp(1)))^2)
# 
# # Mean and variance of Actual using characteristics of a lognormal rv
# worden_df$Worden_Mean <- exp(lognormal_mu + 0.5*lognormal_sig) - 1
# worden_df$Worden_Variance <- (exp(lognormal_sig) - 1)*exp(2*lognormal_mu + lognormal_sig)
# 
# mean((worden_df$Worden_Mean-res_data$Load)^2)

# Variogram plot after kriging/worden method

# vgm_after_kriging <- variog(coords = res_data[res_data$Dataset == "Dataset_1", c("Easting", "Northing")]/1000, 
#                             #data = kriging_df$Kriging_Mean[res_data$Dataset == "Dataset_1"] - res_data$Load[res_data$Dataset == "Dataset_1"], max.dist = 50)
#                             
#                             data = log10((kriging_df$Kriging_Mean[res_data$Dataset == "Dataset_1"] + 1)/(res_data$Load[res_data$Dataset == "Dataset_1"]+1)), max.dist = 50)

# png(file = paste(graphics.path, loss_metric, "_", vgm_model, "_vgm_after_kriging.png", sep = ""),
#     height = 1200,
#     width = 1500,
#     res = 300)
# plot(vgm_after_kriging, ylab = expression(gamma), xlab = "Distance (km)", main = "After kriging (MSLE)")
# dev.off()
# 
# vgm_after_worden <- variog(coords = res_data[res_data$Dataset == "Dataset_1", c("Easting", "Northing")]/1000, 
#data = worden_df$Worden_Mean[res_data$Dataset == "Dataset_1"] - res_data$Load[res_data$Dataset == "Dataset_1"], max.dist = 50)

#                            data = log10((worden_df$Worden_Mean[res_data$Dataset == "Dataset_1"] + 1)/(res_data$Load[res_data$Dataset == "Dataset_1"]+1)), max.dist = 50)
# 
# png(file = paste(graphics.path, loss_metric, "_", vgm_model, "_vgm_after_worden.png", sep = ""),
#     height = 1200,
#     width = 1500,
#     res = 300)
# plot(vgm_after_worden, ylab = expression(gamma), xlab = "Distance (km)", main = "After Worden (MSLE)")
# dev.off()

#--- Comparative plots
##########-----

# Implement method for plot

# Input coordinates of volcano
kelud.df <- data.frame(easting=644183.3,northing=9123213.6) 

# Observation locations
keep_cols <- c("Easting", "Northing", "Load")
obs <- full_data[keep_cols] # Untransformed scale for plotting later.
names(obs) <- (c("x","y","z")) # Set z as the load of ash (our response variable)
obs.df <- obs #  save a dataframe of data points
coordinates(obs) <- c("x", "y") # save as SPDF (Spatial Points Data Frame)
proj4string(obs) <- CRS(sr_utm49) # define projection of data points as UTM Zone 49

# Remove observations outside raster extent:
obs <- obs[obs@coords[, 2] < r.pred@extent@ymax, ]


# Create rasters from kriged and worden output:
kriged_raster <- rasterFromXYZ(kriging_df[, c("Easting", "Northing", "Kriging_Mean")])
worden_raster <- rasterFromXYZ(worden_df[, c("Easting", "Northing", "Worden_Mean")])

# Write raster
writeRaster(kriged_raster, file='kriged_raster.tif', overwrite=T)
writeRaster(worden_raster, file='worden_raster.tif')

kriged_sd_raster <- rasterFromXYZ(kriging_df[, c("Easting", "Northing", "Kriging_Variance")])
kriged_sd_raster <- sqrt(kriged_sd_raster)

worden_sd_raster <- rasterFromXYZ(worden_df[, c("Easting", "Northing", "Worden_Variance")])
worden_sd_raster <- sqrt(worden_sd_raster)

# Uniform uncertainty map for forward prediction:
uniform_sd <- rasterFromXYZ(data.frame("Easting" = kriging_df[, "Easting"],
                                       "Northing" = kriging_df[, "Northing"], 
                                       "RMSE" = sqrt(mean((res_data$Load - res_data$Prediction)^2))))

# Set color pallette
my.palette <- brewer.pal(n = 9, name = "OrRd")
my.palette.2 <- brewer.pal(n = 7, name = "OrRd")
div.palette <- rev(brewer.pal(n = 11, name = "RdYlBu"))
div.palette.2 <- colorRampPalette(c("red", "white", "blue"))
# Plot residual map for areas of over/underprediction

png(file = paste(graphics.path, loss_metric, "_", vgm_model, "_spatial_pred_highres.png", sep = ""),
    height = 1700,
    width = 3000,
    res = 300)

par(mfrow = c(2, 3), mai=c(0.62,0.62,0.42,0.62))

plot(r.pred, 
     col = my.palette,
     breaks= seq(0, 240, by = 30),
     main = "Forward prediction")
contour(r.pred, 
        #nlevels = 8, #number of lines
        #levels = c(10,30, 50, 70, 90), #specify which elevations to put contour lines
        labcex=.5,
        add=TRUE)
points(obs, col="blue", pch=1, cex=obs@data[["z"]]/55)

plot(kriged_raster, 
     col = my.palette,
     breaks= seq(0, 240, by = 30),
     main = "Kriging prediction")
contour(kriged_raster, 
        #nlevels = 8, #number of lines
        #levels = c(10,30, 50, 70, 90), #specify which elevations to put contour lines
        labcex=.5,
        add=TRUE)
points(obs, col="blue", pch=1, cex=obs@data[["z"]]/55)

plot(worden_raster, 
     col = my.palette,
     breaks= seq(0, 240, by = 30),
     main = "Worden method")
contour(worden_raster, 
        #nlevels = 8, #number of lines
        #levels = c(10,30, 50, 70, 90), #specify which elevations to put contour lines
        labcex=.5,
        add=TRUE)
points(obs, col="blue", pch=1, cex=obs@data[["z"]]/55)

plot(uniform_sd, 
     col = my.palette,
     breaks=seq(14, 22, by = 1), 
     main = "Forward MSE")
contour(uniform_sd, 
        #nlevels = 8, #number of lines
        #levels = c(10,30, 50, 70, 90), #specify which elevations to put contour lines
        labcex=.5,
        add=TRUE)
points(obs, col="blue", pch=1, cex=obs@data[["z"]]/55)

plot(kriged_sd_raster, 
     col = my.palette.2, 
     breaks=seq(0, 120, by = 20), #seq(14, 22, by = 1),
     main = "Kriging standard deviation")
contour(kriged_sd_raster, 
        #nlevels = 8, #number of lines
        #levels = c(10,30, 50, 70, 90), #specify which elevations to put contour lines
        labcex=.5,
        add=TRUE)
points(obs, col="blue", pch=1, cex=obs@data[["z"]]/55)

plot(worden_sd_raster, 
     col = my.palette.2, 
     breaks=seq(0, 120, by = 20), # seq(14, 22, by = 1)
     main = "Worden standard deviation")
contour(worden_sd_raster, 
        #nlevels = 8, #number of lines
        #levels = c(10,30, 50, 70, 90), #specify which elevations to put contour lines
        labcex=.5,
        add=TRUE)
points(obs, col="blue", pch=1, cex=obs@data[["z"]]/55)

dev.off()

png(file = paste(graphics.path, loss_metric, "_", vgm_model, "_kriged_resid.png", sep = ""),
    height = 1525,
    width = 1700,
    res = 300)

plot(r.pred - kriged_raster, 
     col = div.palette,
     zlim = c(-110, 110),
     # breaks= seq(-110, 110, by = 20),
     main = "Kriged residuals")
points(obs, col="blue", pch=1, cex=obs@data[["z"]]/55)

dev.off()

png(file = paste(graphics.path, loss_metric, "_", vgm_model, "_worden_resid.png", sep = ""),
    height = 1525,
    width = 1700,
    res = 300)

plot(r.pred - worden_raster, 
     col = div.palette,
     zlim = c(-110, 110),
     # breaks= seq(-110, 110, by = 20),
     main = "Worden residuals")
points(obs, col="blue", pch=1, cex=obs@data[["z"]]/55)

dev.off()

arg <- list(at = c(-1, 0, 1), labels = c("0.1", "1", "10"))
arg.2 <- list(at = c(0.1, 1, 10), labels = c("0.1", "1", "10"))


png(file = paste(graphics.path, loss_metric, "_", vgm_model, "_kriged_resid_ratio.png", sep = ""),
    height = 1525,
    width = 1700,
    res = 300)

plot(log10(kriged_raster/r.pred),  
     zlim= c(-10, 2),
     col = rev(div.palette.2(100)),
     axis.args = arg,
     main = "Ratio between kriged predictions and model predictions")
points(obs, col="blue", pch=1, cex=obs@data[["z"]]/55)

dev.off()

png(file = paste(graphics.path, loss_metric, "_", vgm_model, "_worden_resid_ratio.png", sep = ""),
    height = 1525,
    width = 1700,
    res = 300)

plot(log10(worden_raster/r.pred), 
     #col = div.palette,
     zlim= c(-1.25, 1.25),
     col = rev(div.palette.2(13)),
     axis.args = arg,
     main = "Ratio between Worden predictions and model predictions")
points(obs, col="blue", pch=1, cex=obs@data[["z"]]/55)

dev.off()

plot(worden_raster/r.pred, 
     #col = div.palette,
     zlim= c(10^(-1.25), 10^1.25),
     col = rev(div.palette.2(100)),
     axis.args = arg.2,
     main = "Ratio between Worden predictions and model predictions")
points(obs, col="blue", pch=1, cex=obs@data[["z"]]/55)


save(worden_raster, worden_sd_raster, file = "G:/Tephra_Spatial/Data/msl_matern_wolden_rasters.RData")