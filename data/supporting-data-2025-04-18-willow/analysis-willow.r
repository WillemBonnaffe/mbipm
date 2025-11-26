###############
## LOAD DATA ##
###############

#
###

#######################
## HEIGHT TO BIOMASS ##
#######################

## Load data (Atchley 1989)
dat = read.csv('atchley-1989-table-8.csv', sep=',', header=T)
head(dat)

## Height to biomass (cm to grams) (Atchley 1989)
xx = dat$Height..cm.
yy = log(dat$Shoot.wt...g. + dat$Root.wt...g.)
plot(xx, exp(yy), bty='l', xlim=c(0,300), ylim=c(0,5000))
interpolation = lm(yy ~ xx)
height_to_biomass = function(x, params = interpolation$coefficient) exp(params[1] + params[2] * x)
lines(xx, height_to_biomass(xx), col='red')
print(height_to_biomass(300))

## Shoot to root 
root_to_shoot_mean = mean(1/dat$Shoot.root)
root_to_shoot_sd = sd(1/dat$Shoot.root)
print(root_to_shoot_mean)
print(root_to_shoot_sd)

## Height to biomass (m to kg) (Mosseler et al. 2016)
xx = seq(0, 5, 0.1)
height_to_biomass = function(x) (-0.459 + 0.782 * x ) * 1.3 # Adding below-ground biomass
yy = height_to_biomass(xx)
plot(xx, yy, type='l', xlab='Height (m)', ylab='Biomass (kg)')
lines(xx, rep(0.3, length(xx)), lty=2)
lines(c(-1,1)*100, c(-1,1)*100, lty=3)
print(height_to_biomass(3.0))

#
###

##############
## SURVIVAL ##
##############

## Survival
# 0.3 kg at 2 years old ()
# In first year should be seedlings and experience high mortality
xx = seq(0, 5.0, .01)
yy = function(x, a, b, c) a/(1+exp(-(x-b)*c))
plot(xx,yy(xx, 0.983, 0.1, 5), xlab='Biomass (kg)', ylab='Survival Rate (prop.)', type='l')
mean(yy(seq(0.0,5.0,0.1), 0.983, 0.1, 5))
lines(c(0.65, 0.65), c(0, 1.0), lty=2)

#
###

##################
## REPRODUCTION ##
##################

## Survival
# 0.350 kg in first year
# In first year should be seedlings and experience high mortality
xx = seq(0, 5.0, .01)
yy = function(x, a, b, c) a/(1+exp(-(x-b)*c))
plot(xx,yy(xx, 0.4, 0.75, 10), xlab='Biomass (kg)', ylab='Probability of reproduction (prop.)', type='l')
mean(yy(seq(0.0,5.0,0.1), 0.983, 0.1, 5))
lines(c(0.65, 0.65), c(0, 1.0), lty=2)

#
###

######################
## BROWSING EVASION ##
######################

xx = seq(0, 5.0, 0.01)
yy = function(x, a, b, c) a/(1+exp(-(x-b)*c))
plot(xx, yy(xx, 1.0, 2.0, -10.0), xlab='Biomass (kg)', ylab='Probability of browsing (prop.)', type='l')

#
###