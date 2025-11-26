############################
## ASPEN PARAMETERISATION ##
############################

## Author: Willem Bonnaffé (w.bonnaffe@gmail.com)

############################
## ONTOGENY - FORMAT DATA ##
############################

## Load data
data = read.csv('bella-1980-table-4.csv', sep=',', header=T)
head(data)

## Rename columns
colnames(data) = c('age', 'height', 'density', 'wood.weight', 'leaves.weight', 'total.weight')

## Height ~ Age
plot(data$age, data$height, xlab='Age (years)', ylab='Height (m)', xlim=c(0,6), ylim=c(0,4), type='b', bty='l')

## Biomass ~ Age
plot(data$age, data$total.weight/data$density * 1.2, xlab='Age (years)', ylab='Total Biomass (kg/tree)', xlim=c(0,6), ylim=c(0,0.3), bty='l')

#
###

###############################
## ONTOGENY - AGE TO BIOMASS ##
###############################

## Response and explanatory variable
Y = log(data$total.weight/data$density * 1.2)
X = data$age

## Quadratic interpolation
interpolation = lm(Y ~ X + I(X^2))
print(interpolation$coefficients)

## Age-to-biomass model
age_to_biomass = function(A, bbeta=interpolation$coefficients) exp(bbeta[1] + bbeta[2] * A + bbeta[3] * A^2)

## Visualise
x = seq(0, 10, .1)
plot(data$age, data$total.weight/data$density * 1.2, xlim=c(0,10), ylim=c(0,1), xlab='Age (years)', ylab='Total Biomass (kg/tree)')
lines(x, age_to_biomass(x), col="red")

print(age_to_biomass(0.5))

#
###

####################################
## ONTOGENY - DIAMETER TO BIOMASS ##
####################################

rho = 1318/12.7 # (Bella et al. 1980, p. 6(12), Table 1, ratio between average height and Dbhob)

## Diameter to biomass (Bella et al. 1980)
diameter_to_biomass = function(D, a, b, rho) exp(a + b * log(D^3 * rho)) * 1.2 # With scaling by 1.2 for below ground biomass
a = -0.80319; b = 0.936736; rho = 103.7795

#
###

#########################
## SURVIVAL - ANALYSIS ##
#########################

## Mortalities and diameter at breast height (Zegler et al. 2012)
m = c(0.16, 0.7, 0.8, 0.5, 1/83) # last point is asymptotic survival
s = 1 - m
D = c(1.0, 2.5, 5.0, 10.0, 2774/103.7795) # last point is asymptotic diameter

## Converting to biomass
B = diameter_to_biomass(D,a,b,rho)/1000

## Survival model
sigmoid = function(x, a, b, c) a/(1+exp(-(x - b)*c))

## Visualisation
par(mfrow=c(3,1), bty='l')
plot(B, D, type='b', xlab='Biomass (kg/tree)', ylab='Diameter (cm)')
plot(B, D * rho, type='b', xlab='Biomass (kg/tree)', ylab='Height (cm)')
plot(B, s, type='b', xlab='Biomass (kg/tree)', ylab='Survival Probability')
x = seq(0,max(B),0.1)
lines(x, sigmoid(x, 0.988, -75.0, 0.025), col='blue') # Without Canker disease
lines(x, sigmoid(x, 0.988, 20.0, 0.05), col='red') # With Canker disease
par(mfrow=c(1,1))

#
###

#################
## CONSUMPTION ##
#################

## Nitrogen uptake rate
root_shoot_ratio = 0.2 # (Liepins et al. 2017, p. 63, Fig. 3)
mass_root_system = 2778 # (g root / m^2) (Steele et al. 1997, p. 581(5) Table 4.)
uptake_rate_per_root_mass = 150 * 1e-6 # (g N / g root / h) (Coleman et al. 1998, p. 522(10), Fig. 6)
mass_total = mass_root_system + mass_root_system/root_shoot_ratio # (g)
print(mass_total)
uptake_rate_nitrogen = mass_root_system * uptake_rate_per_root_mass / mass_total # (g N / g plant / h)
uptake_rate_nitrogen = uptake_rate_nitrogen * 24 * 365 * 0.5 # (kg N / kg plant / yr)
print(uptake_rate_nitrogen)

## Nitrogen efficiency
nitrogen_content = (10+25)/2 # (mg N / g plant) (Coleman et al. 1998, p. 522(10), Fig. 5)
nitrogen_content = nitrogen_content * 1e-3 # (kg N / kg plant)
print(nitrogen_content)
relative_growth = 0.03 # (1 / day) (Coleman et al. 1998, p. 522(10), Fig. 5)
nitrogen_efficiency =  (exp(relative_growth)-1)/nitrogen_content * 365 * 0.5 # (kg N / kg / day) (Coleman et al. 1998, p. 517(5))
print(nitrogen_efficiency)

## Biomass uptake rate
uptake_rate_biomass = uptake_rate_nitrogen * nitrogen_efficiency
print(uptake_rate_biomass)

#
###

################################
## REPRODUCTION - FORMAT DATA ##
################################

## Goal: digitise Fig. 6 from Cole et al. (2021)
## DOI: 10.1093/aob/mcaa070

## Load digitised data
scale = read.csv('cole-2021-fig-6-y-axis-pixel-values.csv', header=T, sep=',')
pixel_values = read.csv('cole-2021-fig-6-pixel-values.csv', header=T, sep=',')

## Compute pixels per reproduction
pixels_per_reproduction = abs(mean(diff(scale$y - max(scale$y))))/abs(mean(diff(scale$reproduction)))
# 112

## Convert pixel values to reproduction
reproduction = abs(pixel_values$y - max(pixel_values$y))/pixels_per_reproduction 

## Combine
data_age_reproduction = data.frame(cbind(pixel_values$age, reproduction))
colnames(data_age_reproduction) = c('age', 'reproduction')

## Visualise
plot(data_age_reproduction$age, data_age_reproduction$reproduction, type='b')

## Write csv
write.table(x=data_age_reproduction, file='cole-2021-fig-6.csv', sep=',')

#
###

###################################
## REPRODUCTION - ANALYSIS DATA  ##
###################################

## Load data
data = read.csv('cole-2021-fig-6.csv', header=T, sep=',')
head(data)

## Format data
data$reproduction = data$reproduction/max(data$reproduction) * 0.75 # Scale by max reproduction (Worrell et al. 1999)

## Change colnames
colnames(data) = c('a','r')

## Define function
reproduction = function(x, a, b, c) a/(1+exp(-(x - b)*c))

## Visualise
par(mfrow=c(2,1))
A = seq(0,10,0.1)
W = age_to_biomass(A)
#
plot(data$a, data$r, type='b', xlim=c(0,10), ylim=c(0,1), xlab='Age (years)', ylab='Reproduction Probability')
lines(A, reproduction(W, 0.75, 0.5, 10.0), col='red')
#
plot(age_to_biomass(data$a), data$r, type='b', xlim=c(0,10), ylim=c(0,1), xlab='Total Biomass (kg/tree)', ylab='Reproduction Probability')
lines(W, reproduction(W, 0.75, 0.5, 10.0), col='red')
#
par(mfrow=c(1,1))

#
###

#####################
## METABOLIC COSTS ##
#####################



#
###

#############################
## PREDATION - FORMAT DATA ##
#############################

## Goal: digitise Fig. 7 from Brice et al. (2024)
## DOI:

## Load digitised data
scale_x = read.csv('brice-2024-fig-7-x-axis-pixel-values.csv', header=T, sep=',')
scale_y = read.csv('brice-2024-fig-7-y-axis-pixel-values.csv', header=T, sep=',')
pixel_values = read.csv('brice-2024-fig-7-pixel-values.csv', header=T, sep=',')

## Compute pixels per reproduction
pixels_per_value_x = abs(mean(diff(scale_x$x - max(scale_x$x))))/abs(mean(diff(scale_x$value)))
print(pixels_per_value_x)
pixels_per_value_y = abs(mean(diff(scale_y$y - max(scale_y$y))))/abs(mean(diff(scale_y$value)))
print(pixels_per_value_y)

## Convert pixel values to reproduction
values_x = rev(abs(pixel_values$x - max(pixel_values$x))/pixels_per_value_x)
values_y = abs(pixel_values$y - max(pixel_values$y))/pixels_per_value_y

## Combine
data_combined = data.frame(cbind(values_x, values_y))
colnames(data_combined) = c('height', 'browse.probability')

## Visualise
plot(data_combined$height, data_combined$browse.probability, type='b')

## Write csv
write.table(x=data_combined, file='brice-2024-fig-7.csv', sep=',')

#
###

##########################
## PREDATION - ANALYSIS ##
##########################

## Load data
data = read.csv('brice-2024-fig-7.csv', sep=',', header=T)

## Convert to biomass
H = data$height
D = H / rho
B = diameter_to_biomass(D, a, b, rho) / 1000

## Define function
browse_probability = function(x, a, b, c) a/(1+exp(-(x - b)*c))

## Visualise
par(mfrow=c(2,1))
#
plot(B, data$browse.probability, type='b', xlab='Biomass (kg)', ylab='Browse Probability')
x = seq(0.0, 3.0, 0.01)
lines(x, browse_probability(x, 0.8, 0.3, -6), col='red')
plot(B * 100, data$browse.probability, type='b', xlab='Stand Biomass (kg)', ylab='Browse Probability')
lines(x * 100, browse_probability(x * 100, 0.8, 0.3 * 100, -6 / 100), col='red')
#
par(mfrow=c(1,1))

#
###

i = 1
num_iterations = 10
avg_time = 0.0
time_factor = 1.0 # 60 for minutes
for(i in 1:num_iterations){
  
  ## Start
  t0 = Sys.time()
  
  ## Code
  Sys.sleep(1.0)
  print(i)
  
  ## End
  t1 = Sys.time()
  delta_t = as.numeric(difftime(t1,t0,"seconds"))
  avg_time = mean(c(avg_time, delta_t))
  exp_time = avg_time * (num_iterations - i)
  message(paste("tti:", round(avg_time/time_factor,2), "s"))
  message(paste("ttc:", round(exp_time/time_factor,2), "s"))
}