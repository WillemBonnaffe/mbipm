##############
## INITIATE ##
##############

## Load packages
library(readxl)

## Script label
SCRIPT_LABEL = "A1"
FIG_RES = 300
FIG_WID = 6
FIG_HEI = 6
SAVE_FIG = T
START_YEAR = 1996

## Graphical parameters
par(bty='l')

#
###

######################################
## FORMAT GRASS BIOMASS TIME SERIES ##
######################################

## LOAD AND CHECK

## Load a specific sheet by name
df = data.frame(read_excel("data/master-v2025-05-15-17-06.xlsx", sheet = "grass-biomass"))

# Ensure correct data types
df$year <- as.factor(df$year)
df$biomass <- as.numeric(df$biomass)

## FIGURE
file_name = paste("A1-outputs/", SCRIPT_LABEL, "-fig-grass-biomass-time-series.png", sep="")
if (SAVE_FIG) png(file_name, width=FIG_WID, height=FIG_HEI, res=FIG_RES, units="in")
boxplot(biomass ~ year, data = df,
        xlab = "Year",
        ylab = "Grass Biomass (??)",
        las = 2,         # Rotate x-axis labels
        col = "lightgreen")
if (SAVE_FIG) dev.off()

## SUM FOR EACH YEAR

## Information
length(unique(df$cell)) # 1978 pixels
length(unique(df$year)) # 16 years

## Ensure correct data types
df$year <- as.numeric(as.character(df$year))
head(df$year)

## Compute biomass time series across all cells for each year
biomass_sum_by_year = aggregate(biomass ~ year, data = df, sum)

## Plot
plot(biomass_sum_by_year$year, biomass_sum_by_year$biomass, type='b', xlab='Year', ylab='Grass Biomass (??)', bty='l')

#
###

######################################
## FORMAT WILLOW HEIGHT TIME SERIES ##
######################################

## LOAD AND CHECK DATA

## Load a specific sheet by name
df = data.frame(read_excel("data/master-v2025-05-15-17-06.xlsx", sheet = "willow-heights"))

## Check data
head(df)

## Visualise distribution for single year
s = (df$year == "2001")
length(s)
hist(as.numeric(df$fall_height[s]))

# Ensure year is treated as a factor or sorted
df$year <- as.factor(df$year)

# Convert fall_height to numeric in case it's not already
df$fall_height <- as.numeric(df$fall_height)

## FIGURE
file_name = paste("A1-outputs/", SCRIPT_LABEL, "-fig-willow-height-time-series.png", sep="")
if (SAVE_FIG) png(file_name, width=FIG_WID, height=FIG_HEI, res=FIG_RES, units="in")
boxplot(fall_height ~ year, data = df,
        xlab = "Year",
        ylab = "Fall Height (cm)",
        # main = "Time Series of Willow Fall Heights",
        las = 2,         # Rotate x-axis labels for readability
        col = "lightblue")
if (SAVE_FIG) dev.off()

## COMPUTE DISTRIBUTION TIME SERIES

# Ensure numeric
df$fall_height <- as.numeric(df$fall_height)
df$year <- as.character(df$year)

# Define height bins
bins <- seq(0, ceiling(max(df$fall_height, na.rm = TRUE)), by = 5)
bin_labels <- paste(head(bins, -1), bins[-1], sep = "-")
years <- sort(unique(df$year))

# Build count matrix
mat <- sapply(years, function(y) {
  h <- df$fall_height[df$year == y]
  tab <- table(cut(h, breaks = bins, labels = bin_labels, include.lowest = TRUE))
  as.numeric(tab[bin_labels])  # ensure consistent order
})
rownames(mat) <- bin_labels

# Plot
image(1:ncol(mat), 1:nrow(mat), t(mat), axes = FALSE, col = topo.colors(100),
      xlab = "Year", ylab = "Height Bin (cm)")

# Axes
axis(1, at = 1:ncol(mat), labels = years, las = 2, cex.axis = 0.7)
axis(2, at = 1:nrow(mat), labels = rownames(mat), las = 2, cex.axis = 0.7)

#
###

#####################################
## FORMAT ASPEN HEIGHT TIME SERIES ##
#####################################

## LOAD AND CHECK DATA

## Load a specific sheet by name
df = data.frame(read_excel("data/master-v2025-05-15-17-06.xlsx", sheet = "aspen-heights"))

## Check data
head(df)

## Visualise distribution for single year
s = (df$year == "2001")
length(s)
hist(as.numeric(df$height_cm[s]))

# Ensure year is treated as a factor or sorted
df$year <- as.factor(df$year)

# Convert fall_height to numeric in case it's not already
df$height_cm <- as.numeric(df$height_cm)

## FIGURE
file_name = paste("A1-outputs/", SCRIPT_LABEL, "-fig-aspen-height-time-series.png", sep="")
if (SAVE_FIG) png(file_name, width=FIG_WID, height=FIG_HEI, res=FIG_RES, units="in")
boxplot(height_cm ~ year, data = df,
        xlab = "Year",
        ylab = "Height (cm)",
        # main = "Time Series of Willow Fall Heights",
        las = 2,         # Rotate x-axis labels for readability
        col = "lightblue")
if (SAVE_FIG) dev.off()

#
###

############################
## FORMAT ELK TIME SERIES ##
############################

## LOAD AND CHECK DATA

## Load a specific sheet by name
df = data.frame(read_excel("data/master-v2025-05-15-17-06.xlsx", sheet = "elk-counts-uncorrected"))

## Check data
head(df)

## FIGURE
file_name = paste("A1-outputs/", SCRIPT_LABEL, "-fig-elk-counts-time-series.png", sep="")
if (SAVE_FIG) png(file_name, width=FIG_WID, height=FIG_HEI, res=FIG_RES, units="in")
plot(df$year, df$northern_YNP_count, type="b", xlab="Year", ylab="Counts", bty='l')
if (SAVE_FIG) dev.off()

#
###

##############################
## FORMAT BISON TIME SERIES ##
##############################

## LOAD AND CHECK DATA

## Load a specific sheet by name
df = data.frame(read_excel("data/master-v2025-05-15-17-06.xlsx", sheet = "bison-counts"))

## Check data
head(df)

## FIGURE
file_name = paste("A1-outputs/", SCRIPT_LABEL, "-fig-bison-counts-time-series.png", sep="")
if (SAVE_FIG) png(file_name, width=FIG_WID, height=FIG_HEI, res=FIG_RES, units="in")
plot(df$year, df$northern_YNP_count, type="b", xlab="Year", ylab="Counts", bty='l')
if (SAVE_FIG) dev.off()

#
###

#############################
## FORMAT WOLF TIME SERIES ##
#############################

## LOAD AND CHECK DATA

## Load a specific sheet by name
df = data.frame(read_excel("data/master-v2025-05-15-17-06.xlsx", sheet = "wolf-counts"))

## Check data
head(df)

## FIGURE
file_name = paste("A1-outputs/", SCRIPT_LABEL, "-fig-wolf-counts-time-series.png", sep="")
if (SAVE_FIG) png(file_name, width=FIG_WID, height=FIG_HEI, res=FIG_RES, units="in")
plot(df$year, df$northern_YNP_count, type="b", xlab="Year", ylab="Counts", bty='l')
if (SAVE_FIG) dev.off()

#
###

###############################
## FORMAT COUGAR TIME SERIES ##
###############################

## LOAD AND CHECK DATA

## Load a specific sheet by name
df = data.frame(read_excel("data/master-v2025-05-15-17-06.xlsx", sheet = "cougar-counts"))

## Check data
head(df)

## FIGURE
file_name = paste("A1-outputs/", SCRIPT_LABEL, "-fig-cougar-counts-time-series.png", sep="")
if (SAVE_FIG) png(file_name, width=FIG_WID, height=FIG_HEI, res=FIG_RES, units="in")
plot(df$year, df$northern_YNP_count, type="b", xlab="Year", ylab="Counts")
if (SAVE_FIG) dev.off()

#
###

####################################
## WILLOW HEIGHT TO BIOMASS MODEL ##
####################################

## Height (m) to biomass (kg)
height_to_biomass_willow = function(x, a = -0.459, b = 0.782) {
  height_to_above_ground_biomass_ = (a + b * x ) # (Mosseler et al. 2016, p. 103(7), Table 5)
  height_to_biomass_ = height_to_above_ground_biomass_ * 1.3 # Correcting for below ground biomass (Atchley 1989, p. 41(53), Table 8)
  return(height_to_biomass_)
}

## Convert height (m) to biomass (kg) in willow
df = data.frame(read_excel("data/master-v2025-05-15-17-06.xlsx", sheet = "willow-heights"))
df$calendar_year <- as.factor(df$calendar_year)
df$fall_height <- as.numeric(df$fall_height)
df$biomass <- height_to_biomass_willow(df$fall_height * 0.01) # kg(cm) -> kg(m)
s = which(df$biomass < 0)
length(s)/length(df$biomass) # 181/6043 too small for model
df$biomass[s] = 0
hist(df$biomass)

#
###

#############################
## ASPEN HEIGHT TO BIOMASS ##
#############################

## Ratio height-to-diameter (cm to cm)
rho = mean(c(1318, 1712, 1328)/c(12.7, 16.5, 11.5)) # H/D ratio (Bella et al. 1980, p. 6(12), Table 1)
# [1] 107.6718

## Diameter (cm) to biomass (g) (Bella et al. 1980)
diameter_to_biomass = function(D, a = -0.80319, b = 0.936736, rho = 107.6718) {
  height_to_above_ground_biomass_ = exp(a + b * log(D^3 * rho)) # Regression model (Bella et al. 1980, p. 8, Table 3, row "Y5")
  height_to_biomass_ = height_to_above_ground_biomass_ * 1.2 # With scaling by 1.2 for below ground biomass (Liepins et al. 2017, p. 63(7), Fig. 3)
  return(height_to_biomass_)
}

## Height (cm) to biomass (g) (Bella et al. 1980)
height_to_biomass_aspen = function(H, a = -0.80319, b = 0.936736, rho = 107.6718) {
  height_to_above_ground_biomass_ = exp(a + b * (3 * log(H) - 2 * log(rho))) # (Bella et al. 1980, p. 8, Table 3, row "Y5")
  height_to_biomass_ = height_to_above_ground_biomass_ * 1.2 # With scaling by 1.2 for below ground biomass (Liepins et al. 2017, p. 63(7), Fig. 3)
  return(height_to_biomass_)
}

## Convert height (cm) to biomass (kg) in aspen
df = data.frame(read_excel("data/master-v2025-05-15-17-06.xlsx", sheet = "aspen-heights"))
df$calendar_year <- as.factor(df$calendar_year)
df$height <- as.numeric(df$height)
df$biomass <- height_to_biomass_aspen(df$height)/1000 # g(cm) --> kg(cm)
hist(df$biomass)

#
###

###################
## COMBINED PLOT ##
###################

## Initiate time series
time_series_inference = NULL
ref_years = seq(START_YEAR, 2025, 1)
time_series_inference = cbind(time_series_inference, ref_years)
time_series_inference = cbind(0:(nrow(time_series_inference)-1), time_series_inference)

## FIGURE
file_name = paste("A1-outputs/", SCRIPT_LABEL, "-fig-all-time-series.png", sep="")
if (SAVE_FIG) png(file_name, width=FIG_WID, height=FIG_HEI*2, res=FIG_RES, units="in")
par(mfrow=c(7,1), mar=c(4,4,1,1), bty='l')

## Grass
df = data.frame(read_excel("data/master-v2025-05-15-17-06.xlsx", sheet = "grass-biomass"))
sum_by_year = aggregate(biomass ~ calendar_year, data = df, sum)
s = (sum_by_year$calendar_year >= START_YEAR)
year = sum_by_year$calendar_year[s]
count = sum_by_year$biomass[s] * 1.12085 * 1900 * 100 # lbs/acre --> kg/ha --> kg for total study area (Smith et al. 2023, p. 246(2))
plot(year, count, type='b', xlim=c(START_YEAR, 2025), xlab="Year", ylab='Total Biomass (kg)')

## Collect
s = match(ref_years, year)
time_series_inference = cbind(time_series_inference, count[s])

## Willow
df = data.frame(read_excel("data/master-v2025-05-15-17-06.xlsx", sheet = "willow-heights"))
df$calendar_year <- as.factor(df$calendar_year)
df$fall_height <- as.numeric(df$fall_height)
height_to_biomass = function(x) (-0.459 + 0.782 * x ) * 1.3
df$biomass <- height_to_biomass_willow(df$fall_height * 0.01)
s = which(df$biomass < 0); df$biomass[s] = 0
sum_by_year = aggregate(biomass ~ calendar_year, data = df, mean) # E[Z_willow]
sum_by_year$calendar_year <- as.numeric(as.character(sum_by_year$calendar_year))
s = (sum_by_year$calendar_year >= START_YEAR)
year = sum_by_year$calendar_year[s]
count = sum_by_year$biomass[s]
plot(year, count, type='b', xlim=c(START_YEAR, 2025), xlab="Year", ylab='Mean p.c. Biomass (kg)')

## Collect
s = match(ref_years, year)
time_series_inference = cbind(time_series_inference, count[s])

## Aspen
df = data.frame(read_excel("data/master-v2025-05-15-17-06.xlsx", sheet = "aspen-heights"))
df$calendar_year <- as.factor(df$calendar_year)
df$height_cm <- as.numeric(df$height_cm)
df$biomass <- height_to_biomass_aspen(df$height_cm)/1000 # g(cm) --> kg(cm)
sum_by_year = aggregate(biomass ~ calendar_year, data = df, mean) # E[Z_aspen]
sum_by_year$calendar_year <- as.numeric(as.character(sum_by_year$calendar_year))
s = (sum_by_year$calendar_year >= START_YEAR)
year = sum_by_year$calendar_year[s]
count = sum_by_year$biomass[s]
plot(year, count, type='b', xlim=c(START_YEAR, 2025), xlab="Year", ylab='Mean p.c. Biomass (kg)')

## Collect
s = match(ref_years, year)
time_series_inference = cbind(time_series_inference, count[s])

## Load a specific sheet by name
df = data.frame(read_excel("data/master-v2025-05-15-17-06.xlsx", sheet = "elk-counts-uncorrected"))
s = which(df$calendar_year >= START_YEAR)
year = df$calendar_year[s]
count = df$northern_YNP_count[s]
plot(year, count, type='b', xlim=c(START_YEAR, 2025), xlab="Year", ylab='Count')

## Collect
s = match(ref_years, year)
time_series_inference = cbind(time_series_inference, count[s])

## Load a specific sheet by name
df = data.frame(read_excel("data/master-v2025-05-15-17-06.xlsx", sheet = "bison-counts"))
s = which(df$calendar_year >= START_YEAR)
year = df$calendar_year[s]
count = df$northern_YNP_count[s]
plot(year, count, type='b', xlim=c(START_YEAR, 2025), xlab="Year", ylab='Count')

## Collect
s = match(ref_years, year)
time_series_inference = cbind(time_series_inference, count[s])

## Load a specific sheet by name
df = data.frame(read_excel("data/master-v2025-05-15-17-06.xlsx", sheet = "wolf-counts"))
s = which(df$calendar_year >= START_YEAR)
year = df$calendar_year[s]
count = df$northern_YNP_count[s]
plot(year, count, type='b', xlim=c(START_YEAR, 2025), xlab="Year", ylab='Count')

## Collect
s = match(ref_years, year)
time_series_inference = cbind(time_series_inference, count[s])

## Load a specific sheet by name
df = data.frame(read_excel("data/master-v2025-05-15-17-06.xlsx", sheet = "cougar-counts"))
s = which(df$calendar_year >= START_YEAR)
year = df$calendar_year[s]
count = df$northern_YNP_count[s]
plot(year, count, type='b', xlim=c(START_YEAR, 2025), xlab="Year", ylab='Count')

## Collect
s = match(ref_years, year)
time_series_inference = cbind(time_series_inference, count[s])

par(mfrow=c(1,1))
dev.off()

## Save table
colnames(time_series_inference) = c("time_step","year","biomass_grass", "mean_biomass_willow", "mean_biomass_aspen", "count_elk", "count_bison", "count_wolf", "count_cougar")
file_name = paste("A1-outputs/", SCRIPT_LABEL, "-tab-all-time-series.csv", sep="")
write.csv(x = time_series_inference, file = file_name, row.names = F)

#
###

###########################
## COMPACT VISUALISATION ##
###########################

## FIGURE
file_name = paste("A1-outputs/", SCRIPT_LABEL, "-fig-all-time-series-pairwise.png", sep="")
if (SAVE_FIG) png(file_name, width=FIG_WID, height=FIG_HEI, res=FIG_RES, units="in")
#
labels = c("G", "Wi", "A", "E", "B", "W", "C")
colvect = rev(rainbow(31, start=0.0, end=0.5))[time_series_inference[,1]+1]
par(mfrow=c(7,7), mar=c(1,1,1,1)*.2, bty="o")
#
for (i in 1:7){
  for (j in 1:7){
    if (i == j)
    {
      plot(-1:1, -1:1, xaxt="n", yaxt="n", cex=0)
      text(x=0, y=0, labels=labels[i], cex=2)
    } else {
      plot(time_series_inference[,j+2], time_series_inference[,i+2], type="l", xaxt="n", yaxt="n")
      points(time_series_inference[,j+2], time_series_inference[,i+2], col=colvect, pch=16)
    }
  }
}
#
par(mfrow=c(1,1))
if (SAVE_FIG) dev.off()
  
#
###