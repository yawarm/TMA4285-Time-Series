# Packages
install.packages("ggplot2")
# Load the ggplot2 package
library(ggplot2)
library(reshape2)

#-------
#Reading data

CPI_data <- read.csv('/home/yawar/Downloads/monthly percent-change in CPI - Ark 1.csv') # start from 2010 M01 to 2023 M08
Unemployment_data <- read.csv('/home/yawar/Downloads/Unemployment rate (LFS) - Ark 1.csv') # start from 2010 M01 to 2023 M08
Commodity_consum_data <- read.csv('/home/yawar/Downloads/Totalt varekonsum (Varekonsumindeksen) - Ark 1.csv') # start from 2010 M01 to 2023 M08
Producer_price_index_data <- read.csv('/home/yawar/Downloads/ProducerPriceIndex - Ark 1.csv') # start from 2010 M01 to 2023 M08
#-------

# Function to calculate and plot the Autocorrelation Function (ACF) for a given time series data
# Parameters:
#   data: Time series data vector
#   max_lag: Maximum lag for calculating ACF (default is 40)
#   ci_level: Confidence level for confidence intervals (default is 0.95)
#   main_title: Title for the plot (default is "ACF Plot")
#   xlab: Label for the x-axis (default is "Lag")
#   ylab: Label for the y-axis (default is "ACF")
calculate_and_plot_acf <- function(data, max_lag = 40, ci_level = 0.95, main_title = "ACF Plot", xlab = "Lag", ylab = "ACF") {
  num_obs <- length(data)  # Get the number of observations in the input data
  mean_val <- mean(data)   # Calculate the mean of the input data
  acf_values <- numeric(max_lag + 1)  # Initialize a numeric vector to store ACF values for each lag
  
  # Calculate ACF for each lag up to max_lag
  for (lag in 0:max_lag) {
    numerator_sum <- 0
    denominator1_sum <- 0
    denominator2_sum <- 0
    
    # Calculate the numerator and denominators for the autocorrelation formula
    for (t in (lag + 1):num_obs) {
      numerator_sum <- numerator_sum + (data[t] - mean_val) * (data[t - lag] - mean_val)
      denominator1_sum <- denominator1_sum + (data[t] - mean_val)^2
      denominator2_sum <- denominator2_sum + (data[t - lag] - mean_val)^2
    }
    
    # Calculate and store the autocorrelation value for the current lag
    acf_values[lag + 1] <- numerator_sum / sqrt(denominator1_sum * denominator2_sum)
  }
  
  # Plot ACF
  plot(acf_values, type = "h", col = "blue", lwd = 2, xlab = xlab, ylab = ylab, main = main_title)
  abline(h = 0, col = "black", lty = 2)  # Add a horizontal line at y = 0
  
  # Add confidence intervals
  ci_upper <- qnorm((1 + ci_level) / 2) / sqrt(num_obs)
  ci_lower <- -ci_upper
  abline(h = c(ci_upper, ci_lower), col = "red", lty = 2)  # Add red dashed lines for upper and lower confidence intervals
}


# Function to calculate the Autocorrelation Function (ACF) for a given time series data
# Parameters:
#   data: Time series data vector
#   max_lag: Maximum lag for calculating ACF (default is 40)
# Returns:
#   acf_values: Vector of ACF values for each lag up to max_lag
calculate_acf <- function(data, max_lag = 40) {
  num_obs <- length(data)          # Get the number of observations in the input data
  mean_value <- mean(data)         # Calculate the mean of the input data
  acf_values <- numeric(max_lag + 1)  # Initialize a numeric vector to store ACF values for each lag
  
  # Calculate ACF for each lag up to max_lag
  for (lag in 0:max_lag) {
    numerator_sum <- 0
    denominator1_sum <- 0
    denominator2_sum <- 0
    
    # Calculate the numerator and denominators for the autocorrelation formula
    for (t in (lag + 1):num_obs) {
      numerator_sum <- numerator_sum + (data[t] - mean_value) * (data[t - lag] - mean_value)
      denominator1_sum <- denominator1_sum + (data[t] - mean_value)^2
      denominator2_sum <- denominator2_sum + (data[t - lag] - mean_value)^2
    }
    
    # Calculate and store the autocorrelation value for the current lag
    acf_values[lag + 1] <- numerator_sum / sqrt(denominator1_sum * denominator2_sum)
  }
  
  return(acf_values)
}


# Function to compute the Durbin-Levinson algorithm for partial autocorrelations
# Parameters:
#   acf: Autocorrelation function values
# Returns:
#   Partial autocorrelation coefficients
durbin_levinson <- function(acf){
  order <- length(acf) - 1  # Order of the autoregressive model
  phi <- matrix(nrow = order + 1, ncol = order + 1)  # Matrix to store Durbin-Levinson coefficients
  
  phi[1, 1] <- 0
  phi[2, 2] <- acf[2]
  
  # Iterate through each lag to compute the coefficients
  for (n in 2:order){
    numerator <- acf[n + 1]
    denominator <- 1
    
    # Iterate through previous coefficients to compute new ones
    for (k in 1:(n - 1)){
      if ((n - 1) != k){
        phi[n, k + 1] <- phi[n - 1, k + 1] - phi[n, n] * phi[n - 1, n - k]
        phi[k + 1, n] <- phi[n, k + 1]
      }
      
      numerator <- numerator - phi[n, k + 1] * acf[n - k + 1]
      denominator <- denominator - phi[n, k + 1] * acf[k + 1]
    }
    
    # Compute and store the new coefficient
    phi[n + 1, n + 1] <- as.numeric(numerator / denominator)
  }
  
  # Return only the diagonal values, excluding the first element
  return(diag(phi)[-1])
}

# Function to calculate and plot the Partial Autocorrelation Function (PACF)
# Parameters:
#   vals: Time series values
#   acf: Autocorrelation function values
#   ci_level: Confidence level for confidence intervals (default is 0.95)
#   main_title: Title for the plot (default is "PACF Plot")
#   xlab: Label for the x-axis (default is "Lag")
#   ylab: Label for the y-axis (default is "PACF")
calculate_and_plot_pacf <- function(vals, acf, ci_level = 0.95, main_title = "PACF Plot", xlab = "Lag", ylab = "PACF"){
  order <- length(acf) - 1  # Order of the autoregressive model
  T <- length(vals)  # Length of the time series
  
  # Compute the Durbin-Levinson coefficients
  coeffs <- durbin_levinson(acf)
  
  # Calculate confidence intervals
  ci_upper <- qnorm((1 + ci_level) / 2) / sqrt(T)
  ci_lower <- -ci_upper
  
  # Plot PACF with confidence intervals and a horizontal line at 0
  plot(coeffs, type = "h", col = "blue", lwd = 2, xlab = xlab, ylab = ylab, main = main_title)
  abline(h = c(ci_upper, ci_lower), col = "red", lty = 2)
  abline(h = 0, col = "black", lty = 2)
}

# Function to calculate the Partial Autocorrelation Function (PACF)
# Parameters:
#   vals: Time series values
#   acf: Autocorrelation function values
# Returns:
#   Partial autocorrelation coefficients
calculate_pacf <- function(vals, acf){
  pacf <- durbin_levinson(acf)
  return(pacf)
}



#---------------

#Elementary data exploration of time series data

# Extracting relevant columns and converting to numeric values
cpi_values <- as.numeric(CPI_data$X0.1)
unemployment <- as.numeric(Unemployment_data$X3.9)
consum <- as.numeric(Commodity_consum_data$X111.9)
producer <- as.numeric(Producer_price_index_data$X84.3)

# Create a time vector based on the length of CPI values
time <- 1:length(cpi_values)

# Display summary statistics for each dataset
summary(cpi_values)
summary(unemployment)
summary(consum)
summary(producer)

# Plotting the original time series data
plot(time, cpi_values, type = "l", xlab = "Time", ylab = "CPI", main = "Original CPI Time Series")
plot(time, unemployment, type = "l", xlab = "Time", ylab = "Original Data", main = "Original Unemployment Time Series")
plot(time, consum, type = "l", xlab = "Time", ylab = "Original Data", main = "Original Commodity Consum Time Series")
plot(time, producer, type = "l", xlab = "Time", ylab = "Original Data", main = "Original Producer Price Index Time Series")

# Create Histogram and Density Plot for CPI values
ggplot(data = CPI_data, aes(x = cpi_values)) +
  geom_histogram(binwidth = 0.2, fill = "blue", color = "black") +
  geom_density(color = "red") +
  labs(title = "Histogram and Density Plot for CPI", x = "CPI Value") +
  theme_minimal()

# Create Histogram and Density Plot for Unemployment values
ggplot(data = Unemployment_data, aes(x = unemployment)) +
  geom_histogram(binwidth = 0.2, fill = "blue", color = "black") +
  geom_density(color = "red") +
  labs(title = "Histogram and Density Plot for Unemployment", x = "Unemployment Value") +
  theme_minimal()


# ACF (Autocorrelation Function) Analysis

# Calculate and plot ACF for CPI data
calculate_and_plot_acf(cpi_values, max_lag = 20, ci_level = 0.95, main_title = "ACF Plot for CPI Data")
# Calculate ACF values for CPI data
acf_cpi_untrans <- calculate_acf(cpi_values, max_lag = 20)

# Calculate and plot ACF for Unemployment data
calculate_and_plot_acf(unemployment, max_lag = 20, ci_level = 0.95, main_title = "ACF Plot for Unemployment Data")
# Calculate ACF values for Unemployment data
acf_unemployment_untrans <- calculate_acf(unemployment, max_lag = 20)

# Calculate and plot ACF for Commodity Consumption data
calculate_and_plot_acf(consum, max_lag = 20, ci_level = 0.95, main_title = "ACF Plot for Commodity Consum Data")
# Calculate ACF values for Commodity Consumption data
acf_consum_untrans <- calculate_acf(consum, max_lag = 20)

# Calculate and plot ACF for Producer Price Index data
calculate_and_plot_acf(producer, max_lag = 20, ci_level = 0.95, main_title = "ACF Plot for Producer Price Index Data")
# Calculate ACF values for Producer Price Index data
acf_producer_untrans <- calculate_acf(producer, max_lag = 20)


# PACF (Partial Autocorrelation Function) Analysis

# Calculate and plot PACF for CPI data
calculate_and_plot_pacf(cpi_values, acf_cpi_untrans, ci_level = 0.95, main_title = "PACF Plot for CPI Data")
# Calculate PACF values for CPI data
calculate_pacf(cpi_values, acf_cpi_untrans)

# Calculate and plot PACF for Unemployment data
calculate_and_plot_pacf(unemployment, acf_unemployment_untrans, ci_level = 0.95, main_title = "PACF Plot for Unemployment Data")
# Calculate PACF values for Unemployment data
calculate_pacf(unemployment, acf_unemployment_untrans)

# Calculate and plot PACF for Commodity Consumption data
calculate_and_plot_pacf(consum, acf_consum_untrans, ci_level = 0.95, main_title = "PACF Plot for Commodity Consum Data")
# Calculate PACF values for Commodity Consumption data
calculate_pacf(consum, acf_consum_untrans)

# Calculate and plot PACF for Producer Price Index data
calculate_and_plot_pacf(producer, acf_producer_untrans, ci_level = 0.95, main_title = "PACF Plot for Producer Price Index Data")
# Calculate PACF values for Producer Price Index data
calculate_pacf(producer, acf_producer_untrans)


# Heat Map Analysis

# Combine data into a single data frame
combined_data <- data.frame(
  CPI = CPI_data$X0.1,
  Unemployment = Unemployment_data$X3.9,
  CommodityConsumption = as.numeric(Commodity_consum_data$X111.9),
  ProducerPriceIndex = as.numeric(Producer_price_index_data$X84.3)
)

# Calculate the covariance matrix
cov_matrix <- cov(combined_data)

# Melt the matrix for ggplot
melted_cov_matrix <- melt(cov_matrix)

# Create a ggplot heatmap
ggplot(melted_cov_matrix, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white", size = 0.5) +  # Add border to tiles
  scale_fill_gradientn(colors = colorRampPalette(c("blue", "white", "red"))(100)) +  # Define color gradient
  labs(title = "Covariance Heatmap",
       x = "Variables",
       y = "Variables") +  # Add titles
  theme_minimal() +  # Use minimal theme
  theme(
    axis.text = element_text(size = 10),  # Adjust text size for axis
    axis.title = element_text(size = 12, face = "bold")  # Adjust title size and style for axis
  )




#---------------
# Data Transformation

# Calculate the natural logarithm (ln) of the data
log_cpi <- (cpi_values)
log_unemployment <- log(unemployment)
log_consum <- log(consum)
log_producer <- log(producer)

# Initialize a vector for log differences
log_diff_cpi <- numeric(length(log_cpi) - 1)
log_diff_unemployment <- numeric(length(log_unemployment) - 1)
log_diff_consum <- numeric(length(log_consum) - 1)
log_diff_producer <- numeric(length(log_producer) - 1)

# Calculate log differences
for (i in 1:(length(log_cpi) - 1)) {
  log_diff_cpi[i] <- log_cpi[i + 1] - log_cpi[i]
}
for (i in 1:(length(log_unemployment) - 1)) {
  log_diff_unemployment[i] <- log_unemployment[i + 1] - log_unemployment[i]
}
for (i in 1:(length(log_consum) - 1)) {
  log_diff_consum[i] <- log_consum[i + 1] - log_consum[i]
}
for (i in 1:(length(log_producer) - 1)) {
  log_diff_producer[i] <- log_producer[i + 1] - log_producer[i]
}

# Create a time vector for log differences
time_log_diff <- 1:length(log_diff_cpi)


# Plot Transformed Data

# Plot Log Differences for CPI
plot(time_log_diff, log_diff_cpi, type = "l", col = "blue", xlab = "Time Period", ylab = "Log Differences CPI", main = "Log Differences CPI")

# Plot Log Differences for Unemployment
plot(time_log_diff, log_diff_unemployment, type = "l", col = "blue", xlab = "Time Period", ylab = "Log Differences Unemployment", main = "Log Differences Unemployment")

# Plot Transformed Commodity Consumption
plot(time_log_diff, log_diff_consum, type = "l", col = "blue", xlab = "Time Period", ylab = "Log Differences Commodity Consumption", main = "Log Differences Commodity Consumption")

# Plot Transformed Producer Price Index
plot(time_log_diff, log_diff_producer, type = "l", col = "blue", xlab = "Time Period", ylab = "Log Differences Producer Price Index", main = "Log Differences Producer Price Index")

# Histograms and Density Plots

# Create data frames for transformed data
CPI_data_trans <- data.frame(time_log_diff, log_diff_cpi)
Unemployment_data_trans <- data.frame(time_log_diff, log_diff_unemployment)
Commodity_consum_data_trans <- data.frame(time_log_diff, log_diff_consum)
Producer_Price_Index_trans <- data.frame(time_log_diff, log_diff_producer)

# Plot Histogram and Density for Transformed CPI
ggplot(data = CPI_data_trans, aes(x = log_diff_cpi)) +
  geom_histogram(binwidth = 0.2, fill = "blue", color = "black") +
  geom_density(color = "red") +
  labs(title = "Histogram and Density Plot", x = "Transformed CPI Value") +
  theme_minimal()

# Plot Histogram and Density for Transformed Unemployment
ggplot(data = Unemployment_data_trans, aes(x = log_diff_unemployment)) +
  geom_histogram(binwidth = 0.2, fill = "blue", color = "black") +
  geom_density(color = "red") +
  labs(title = "Histogram and Density Plot", x = "Transformed Unemployment Value") +
  theme_minimal()

# Plot Histogram and Density for Transformed Commodity Consumption
ggplot(data = Commodity_consum_data_trans, aes(x = log_diff_consum)) +
  geom_histogram(binwidth = 0.2, fill = "blue", color = "black") +
  geom_density(color = "red") +
  labs(title = "Histogram and Density Plot", x = "Transformed Commodity Consum Value") +
  theme_minimal()

# Plot Histogram and Density for Transformed Producer Price Index
ggplot(data = Producer_Price_Index_trans, aes(x = log_diff_producer)) +
  geom_histogram(binwidth = 0.2, fill = "blue", color = "black") +
  geom_density(color = "red") +
  labs(title = "Histogram and Density Plot", x = "Transformed Producer Price Index Value") +
  theme_minimal()


# Summary Statistics

# Summary for Log Differences of CPI
summary(log_diff_cpi)

# Summary for Log Differences of Unemployment
summary(log_diff_unemployment)

# Summary for Log Differences of Commodity Consumption
summary(log_diff_consum)

# Summary for Log Differences of Producer Price Index
summary(log_diff_producer)


# Autocorrelation Function (ACF) Analysis

# ACF Plot for Transformed CPI
calculate_and_plot_acf(log_diff_cpi, max_lag = 20, ci_level = 0.95, main_title = "ACF Plot for CPI Transformed")
acf_cpi_trans <- calculate_acf(log_diff_cpi, max_lag = 20)

# ACF Plot for Transformed Unemployment
calculate_and_plot_acf(log_diff_unemployment, max_lag = 20, ci_level = 0.95, main_title = "ACF Plot for Unemployment Transformed")
acf_unemployment_trans <- calculate_acf(log_diff_unemployment, max_lag = 20)

# ACF Plot for Transformed Commodity Consumption
calculate_and_plot_acf(log_diff_consum, max_lag = 20, ci_level = 0.95, main_title = "ACF Plot for Commodity Consum Transformed")
acf_consum_trans <- calculate_acf(log_diff_consum, max_lag = 20)

# ACF Plot for Transformed Producer Price Index
calculate_and_plot_acf(log_diff_producer, max_lag = 20, ci_level = 0.95, main_title = "ACF Plot for Producer Price Index Transformed")
acf_producer_trans <- calculate_acf(log_diff_producer, max_lag = 20)


# Partial Autocorrelation Function (PACF) Analysis

# PACF Plot for Transformed CPI
calculate_and_plot_pacf(log_diff_cpi, acf_cpi_trans, main_title = "PACF Plot for CPI Transformed")
calculate_pacf(log_diff_cpi, acf_cpi_trans)

# PACF Plot for Transformed Unemployment
calculate_and_plot_pacf(log_diff_unemployment, acf_unemployment_trans, ci_level = 0.95, main_title = "PACF Plot for Unemployment Transformed")
calculate_pacf(log_diff_unemployment, acf_unemployment_trans)

# PACF Plot for Transformed Commodity Consumption
calculate_and_plot_pacf(log_diff_consum, acf_consum_trans, ci_level = 0.95, main_title = "PACF Plot for Commodity Consum Transformed")
calculate_pacf(log_diff_consum, acf_consum_trans)

# PACF Plot for Transformed Producer Price Index
calculate_and_plot_pacf(log_diff_producer, acf_producer_trans, ci_level = 0.95, main_title = "PACF Plot for Producer Price Index Transformed")
calculate_pacf(log_diff_producer, acf_producer_trans)


# Heat Map

# Combine transformed data into a single data frame
combined_trans_data <- data.frame(
  CPI = log_diff_cpi,
  Unemployment = log_diff_unemployment,
  CommodityConsumption = log_diff_consum,
  Producer_Price_Index = log_diff_producer
)

# Calculate the covariance matrix for transformed data
cov_matrix_trans <- cov(combined_trans_data)

# Melt the covariance matrix for ggplot
melted_cov_matrix_trans <- melt(cov_matrix_trans)

# Create a ggplot heatmap for the covariance matrix
ggplot(melted_cov_matrix_trans, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white", size = 0.5) +  # Add border to tiles
  scale_fill_gradientn(colors = colorRampPalette(c("blue", "white", "red"))(100)) +
  labs(title = "Covariance Heatmap",
       x = "Variables",
       y = "Variables") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10),  # Adjust text size for axis
    axis.title = element_text(size = 12, face = "bold")  # Adjust title size and style for axis
  )



# Parameter Estimation Including Uncertainty

# Selecting data for parameter estimation
data <- log_diff_cpi  # Dependent variable
cov1 <- log_diff_unemployment  # First covariate
cov2 <- log_diff_consum  # Second covariate
cov3 <- log_diff_producer  # Third covariate


# ARIMA/ARMAX(1,1,1) Likelihood Function

# ARIMA/ARMAX(1,1,1) Likelihood Function

# Function to compute the log-likelihood for an ARIMA/ARMAX(1,1,1) model
# Parameters: 
# - parameters: Model parameters (phi, theta, sigma2, beta1, beta2, beta3)
# - data: Time series data
# - cov1, cov2, cov3: Covariates for enhancing the model

likelihood_ARIMA_ARMAX <- function(parameters) {
  # Extract ARIMA model parameters
  ar_coefficient <- parameters["phi"]
  ma_coefficient <- parameters["theta"]
  innovation_variance <- parameters["sigma2"]
  
  # Extract covariate coefficients for ARMAX
  covar_coefficient1 <- parameters["beta1"]
  covar_coefficient2 <- parameters["beta2"]
  covar_coefficient3 <- parameters["beta3"]
  
  # Calculate the state space representation
  transition_matrix <- ar_coefficient
  observation_matrix <- 1
  observation_variance <- innovation_variance
  state_innovation_matrix <- ma_coefficient
  
  initial_state_mean <- 0
  initial_state_covariance <- innovation_variance / (1 - ar_coefficient^2)
  
  # Initialize variables for the Kalman filter
  num_obs <- length(data)
  state_estimates <- matrix(0, num_obs, 1)
  state_covariances <- matrix(0, num_obs, 1)
  log_likelihood <- 0
  
  for (time_point in 1:num_obs) {
    # Prediction step
    prediction_error <- data[time_point] - observation_matrix * initial_state_mean - covar_coefficient1 * cov1[time_point] - covar_coefficient2 * cov2[time_point] - covar_coefficient3 * cov3[time_point]
    prediction_error_variance <- observation_matrix * initial_state_covariance * observation_matrix + observation_variance + covar_coefficient1^2 * var(cov1) + covar_coefficient2^2 * var(cov2) + covar_coefficient3^2 * var(cov3)
    kalman_gain <- initial_state_covariance * observation_matrix / prediction_error_variance
    
    # Update state estimates
    updated_state_mean <- initial_state_mean + kalman_gain * prediction_error
    updated_state_covariance <- initial_state_covariance - kalman_gain * prediction_error_variance * kalman_gain
    
    # Update log-likelihood
    log_likelihood <- log_likelihood - 0.5 * (log(2 * pi) + log(prediction_error_variance) + prediction_error^2 / prediction_error_variance)
    
    # Store state estimates
    state_estimates[time_point] <- updated_state_mean
    state_covariances[time_point] <- updated_state_covariance
    
    # Time update for the next iteration
    initial_state_mean <- transition_matrix * updated_state_mean
    initial_state_covariance <- transition_matrix * updated_state_covariance * transition_matrix + state_innovation_matrix * observation_variance * state_innovation_matrix
  }
  
  # Return negative log-likelihood to be minimized
  return(-log_likelihood)
}



# GARCH(1,1) Likelihood Function

# Function to compute the log-likelihood for a GARCH(1,1) model with covariates
# Parameters: 
# - params: GARCH model coefficients (omega, alpha, beta, beta1, beta2, beta3)
# - data: Time series data
# - cov1, cov2, cov3: Covariates for enhancing the model

likelihood_GARCH <- function(params, data, cov1, cov2, cov3) {
  # Extract GARCH model parameters
  omega <- params[1]
  alpha <- params[2]
  beta <- params[3]
  beta1 <- params[4]
  beta2 <- params[5]
  beta3 <- params[6]
  
  # Initialize variables
  T <- length(data)
  sigma_sq <- numeric(T)
  loglik <- 0
  
  # Set initial value for sigma_sq[1]
  sigma_sq[1] <- var(data) / (1 - alpha - beta)
  
  # Calculate the log-likelihood
  for (t in 2:T) {
    sigma_sq[t] <- omega + alpha * data[t - 1]^2 + beta * sigma_sq[t - 1] + beta1 * data[t - 1] * cov1[t] + beta2 * data[t - 1] * cov2[t] + beta3 * data[t - 1] * cov3[t]
    loglik <- loglik - 0.5 * (log(sigma_sq[t]) + (data[t]^2) / sigma_sq[t])
  }
  
  return(-loglik)
}


#-------------

# ARIMA parameter estimation

# Initial parameters for ARIMA model
initial_parameters_ARIMA <- c(phi = 0.1, theta = 0.1, sigma2 = 0.1, beta1 = 0.1, beta2 = 0.1, beta3 = 0.1)

# Optimize ARIMA parameters using likelihood function
optim_results_ARIMA <- optim(par = initial_parameters_ARIMA, fn = likelihood_ARIMA_ARMAX, method = "L-BFGS-B", hessian = TRUE)
ARIMA_estimates <- optim_results_ARIMA$par

# Display ARIMA parameter estimates
ARIMA_estimates

hessian_matrix_ARIMA <- optim_results_ARIMA$hessian

# Calculate covariance matrix and standard errors for ARIMA estimates
covariance_matrix_ARIMA <- solve(hessian_matrix_ARIMA)
std_errors_ARIMA <- sqrt(diag(covariance_matrix_ARIMA))
std_errors_ARIMA

# GARCH parameter estimation

# Initial parameters for GARCH model
initial_parameters_GARCH <- c(omega = 0.1, alpha = 0.5, beta = 0.1, beta1 = 0.1, beta2 = 0.1, beta3 = 0.1)

# Optimize GARCH parameters using likelihood function
optim_results_GARCH <- optim(par = initial_parameters_GARCH, fn = likelihood_GARCH, data = data, cov1 = cov1, cov2 = cov2, cov3 = cov3, method = "BFGS", hessian = TRUE)
GARCH_estimates <- optim_results_GARCH$par

# Display GARCH parameter estimates
GARCH_estimates

hessian_matrix_GARCH <- optim_results_GARCH$hessian

# Calculate covariance matrix and standard errors for GARCH estimates
covariance_matrix_GARCH <- solve(hessian_matrix_GARCH)
std_errors_GARCH <- sqrt(diag(covariance_matrix_GARCH))
std_errors_GARCH



#-----------
# Set the significance level (e.g., 0.05 for a 5% significance level)
alpha <- 0.05

# Determination of significance for the ARIMA(1,1,1) and GARCH(1,1) estimations

# Calculate the t-statistics for each parameter using the bootstrap estimates
t_statistics_ARIMA <- (ARIMA_estimates - initial_parameters_ARIMA) / std_errors_ARIMA
t_statistics_GARCH <- (GARCH_estimates - initial_parameters_GARCH) / std_errors_GARCH

# Calculate the degrees of freedom for the t-distribution
df_ARIMA <- length(data) - length(initial_parameters_ARIMA)  # Assuming that's the correct degrees of freedom
df_GARCH <- length(data) - length(initial_parameters_GARCH)  # Assuming that's the correct degrees of freedom

# Calculate the two-tailed p-values
p_values_ARIMA <- 2 * (1 - pt(abs(t_statistics_ARIMA), df_ARIMA))
p_values_GARCH <- 2 * (1 - pt(abs(t_statistics_GARCH), df_GARCH))

# Check which parameters are statistically significant for ARIMA(1,1,1)
significant_parameters_ARIMA <- data.frame(
  Parameter = names(initial_parameters_ARIMA),
  Coefficient = ARIMA_estimates,
  T_Statistic = t_statistics_ARIMA,
  P_Value = p_values_ARIMA,
  Significant = p_values_ARIMA < alpha
)

# Check which parameters are statistically significant for GARCH(1,1)
significant_parameters_GARCH <- data.frame(
  Parameter = names(initial_parameters_GARCH),
  Coefficient = GARCH_estimates,
  T_Statistic = t_statistics_GARCH,
  P_Value = p_values_GARCH,
  Significant = p_values_GARCH < alpha
)

# Display the results for ARIMA(1,1,1)
significant_parameters_ARIMA

# Display the results for GARCH(1,1)
significant_parameters_GARCH




#-------------

# Calculate AIC for ARIMA(1,1,1)
n_ARIMA <- length(data)
loglik_ARIMA <- -optim_results_ARIMA$value  # Negative log-likelihood from optimization
k_ARIMA <- length(initial_parameters_ARIMA)  # Number of parameters in ARIMA(1,1,1)
AIC_ARIMA <- 2 * k_ARIMA - 2 * loglik_ARIMA

# Calculate AIC for GARCH(1,1)
loglik_GARCH <- -optim_results_GARCH$value  # Negative log-likelihood from optimization
k_GARCH <- length(initial_parameters_GARCH)  # Number of parameters in GARCH(1,1)
AIC_GARCH <- 2 * k_GARCH - 2 * loglik_GARCH

# Compare AIC values
cat("AIC for ARIMA(1,1,1):", AIC_ARIMA, "\n")
cat("AIC for GARCH(1,1):", AIC_GARCH, "\n")

# Model selection
if (AIC_ARIMA < AIC_GARCH) {
  cat("ARIMA(1,1,1) is preferred.\n")
} else if (AIC_GARCH < AIC_ARIMA) {
  cat("GARCH(1,1) is preferred.\n")
} else {
  cat("Both models are equally good.\n")
}




#------------------------

# Forecasting of the inflation rate using GARCH(1,1)

data <- log_diff_cpi
length(data)
# Number of runs
n_runs <- 100

# Initialize a matrix to store forecast_values for each run
all_forecast_values <- matrix(NA, nrow = length(data), ncol = n_runs)

# Move the declaration of forecast_values outside the loop
forecast_values <- numeric(length(data))

for (run in 1:n_runs) {
  # Simulate GARCH(1,1) conditional variances
  sigma_sq <- numeric(length(data))
  sigma_sq[1] <- GARCH_estimates["omega"] / (1 - GARCH_estimates["alpha"] - GARCH_estimates["beta"])
  
  for (t in 2:length(data)) {
    sigma_sq[t] <- GARCH_estimates[1] +
      GARCH_estimates[2] * data[t-1]^2 +
      GARCH_estimates[3] * sigma_sq[t-1] +
      GARCH_estimates[5] * data[t-1] * cov2[length(cov2)] +
      GARCH_estimates[6] * data[t-1] * cov3[length(cov3)]
  }

  # Simulate innovations
  innovations <- rnorm(length(data) + length(forecast_values), mean = 0, sd = sqrt(abs(sigma_sq)))
  
  # Forecast future values
  forecast_values <- innovations[length(data) + 1:length(forecast_values)]
  
  # Store forecast_values in the matrix
  all_forecast_values[, run] <- forecast_values
}

n_sim <- 50

# Create a sequence of time periods for your data and forecasts
time_periods <- 1:(length(data) + n_sim)

# Plot the seasonal differenced data in black
plot(time_periods, c(data, rep(NA, n_sim)), type = "l", col = "black", xlab = "Time", ylab = "CPI Value", main = "Forecasted montly change in CPI using a GARCH(1,1) model")

# Add the forecasted inflation in blue
lines(time_periods[length(data)+1:length(all_forecast_values[, 1])], all_forecast_values[, 1], col = "blue")
lines(time_periods[length(data)+1:length(all_forecast_values[, 2])], all_forecast_values[, 2], col = "red")
lines(time_periods[length(data)+1:length(all_forecast_values[, 3])], all_forecast_values[, 3], col = "orange")
lines(time_periods[length(data)+1:length(all_forecast_values[, 4])], all_forecast_values[, 4], col = "green")



# Calculate the average forecast
average_forecast <- rowMeans(all_forecast_values)

# Calculate standard deviation
forecast_std <- apply(all_forecast_values, 1, sd)

# Calculate other measures like minimum, maximum, etc.
forecast_min <- apply(all_forecast_values, 1, min)
forecast_max <- apply(all_forecast_values, 1, max)

# Display the results
result_table <- data.frame(
  Average_Forecast = average_forecast,
  Forecast_StdDev = forecast_std,
  Forecast_Min = forecast_min,
  Forecast_Max = forecast_max
)

confidence_level_95 <- 0.95
confidence_level_50 <- 0.5
lower_bound_95 <- apply(all_forecast_values, 1, function(x) quantile(x, (1 - confidence_level_95) / 2))
upper_bound_95 <- apply(all_forecast_values, 1, function(x) quantile(x, 1 - (1 - confidence_level_95) / 2))
lower_bound_50 <- apply(all_forecast_values, 1, function(x) quantile(x, (1 - confidence_level_50) / 2))
upper_bound_50 <- apply(all_forecast_values, 1, function(x) quantile(x, 1 - (1 - confidence_level_50) / 2))

# Assuming you have already calculated result_table as mentioned in the previous code

# Plot average forecast with confidence intervals
plot(time_periods, c(data, rep(NA, n_sim)), type = "l", col = "black", xlab = "Time", ylab = "CPI Value", main = "Forecasted monthly change in CPI using GARCH(1,1) model")

# Plot average forecast
lines(time_periods[length(data) + 1:length(average_forecast)], average_forecast, col = "blue")

# Plot 95% confidence intervals
lines(time_periods[length(data) + 1:length(average_forecast)], lower_bound_95, col = "red", lty = 2)
lines(time_periods[length(data) + 1:length(average_forecast)], upper_bound_95, col = "red", lty = 2)
lines(time_periods[length(data) + 1:length(average_forecast)], lower_bound_50, col = "orange", lty = 2)
lines(time_periods[length(data) + 1:length(average_forecast)], upper_bound_50, col = "orange", lty = 2)

# Add legend
legend("topleft", legend = c("Observed", "Average Forecast", "95% Confidence Intervals", "50% Confidence Intervals"), col = c("black", "blue", "red", 'orange'), lty = c(1, 1, 2, 2))

