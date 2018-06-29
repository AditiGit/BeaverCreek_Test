# Process Picarro data for Beaver Creek test experiment
# This workhorse script summarizes individual (raw) Picarro observations to 
# summaries of "samples" (groups of consecutive observations made from a given 
# core at a point in time). It computes gas concentration changes, performs 
# some QC, merges the Picarro data with valve map and other ancillary data,
# and writes SUMMARYDATA_FILE.
# 
# Ben Bond-Lamberty and Aditi Sengupta June 2018

source("0-functions.R")

SCRIPTNAME  	<- "3-summarize.R"
PROBLEM       <- FALSE


# ==============================================================================
# Main 

openlog(file.path(outputdir(), paste0(SCRIPTNAME, ".log.txt")), sink = TRUE) # open log

printlog("Welcome to", SCRIPTNAME)

printlog("Reading in raw data...")
read_csv(RAWDATA_FILE) %>%
  # Convert date/time to POSIXct
  mutate(DATETIME = ymd_hms(paste(DATE, TIME))) %>%
  select(-DATE, -TIME) %>%
  arrange(DATETIME) %>%
  print_dims("rawdata") ->
  rawdata
print(summary(rawdata))
printlog("First timestamp:")
print(min(rawdata$DATETIME))
printlog("Last timestamp:")
print(max(rawdata$DATETIME))

# -----------------------------------------------------------------------------
# Prep work: data cleaning, dates, sample numbers, elapsed time

# Assign a different sample number to each sample group 
# (we know we're on a new sample when MPVPosition changes)
printlog("Assigning sample numbers and computing elapsed time...")
rawdata %>%
  mutate(newsample = MPVPosition != lag(MPVPosition)) %>%
  replace_na(list(newsample = FALSE)) %>% 
  mutate(samplenum = cumsum(newsample)) %>%
  select(-newsample) %>%
  group_by(samplenum) %>%
  mutate(elapsed_seconds = as.double(difftime(DATETIME, min(DATETIME), units = "secs"))) ->
  rawdata_samples

printlog("Reading key data...")
library(readxl)
keydata <- readxl::read_excel("Sample_ValveKey.xlsx") %>% 
  mutate(Valve = as.integer(Valve))

printlog("Removing ambient samples amd matching with valve data...")
AMBIENT_VALVE <- 16
rawdata_samples %>%
  ungroup %>% 
  filter(MPVPosition != AMBIENT_VALVE) %>% 
  left_join(keydata, by = c("MPVPosition" = "Valve")) %>% 
  group_by(samplenum) %>%
  #  filter(max(elapsed_seconds) <= MAX_MEASUREMENT_TIME) %>%
  print_dims("rawdata_samples") ->
  rawdata_samples

# Not sure why there are duplicate CO2 concentrations, but remove
rawdata_samples %>% 
  filter(elapsed_seconds < 60) %>% 
  group_by(samplenum, MPVPosition, Sample, elapsed_seconds) %>% 
  summarise(CO2_dry = mean(CO2_dry), CH4_dry = mean(CH4_dry), DATETIME = mean(DATETIME)) ->
  rawdata_samples_trunc

unmatched <- filter(rawdata_samples_trunc, is.na(Sample))
printlog(nrow(unmatched), "unmatched observations; valve numbers", 
         paste(unique(unmatched$MPVPosition), collapse = " "))
rawdata_samples_trunc <- filter(rawdata_samples_trunc, !is.na(Sample))

printlog("Visualizing...")
p <- ggplot(rawdata_samples_trunc, aes(x = elapsed_seconds, color = Sample, group = samplenum)) +
  ggtitle("Concentration by sampling time and valve")
print(p + geom_line(aes(y = CO2_dry)) + xlim(c(0, 60)))
save_plot("co2_by_valve", ptype = ".png")
print(p + geom_line(aes(y = CH4_dry)) + xlim(c(0, 60)))
save_plot("ch4_by_valve", ptype = ".png")

printlog("Computing concentration change rates...")
# We use a robust linear regression here
rawdata_samples_trunc %>% 
  filter(elapsed_seconds >= MIN_FLUXCALC_TIME, elapsed_seconds <= MAX_FLUXCALC_TIME) %>% 
  group_by(samplenum) %>% 
  do(m = lm(CO2_dry ~ elapsed_seconds, data = .)) %>% 
  summarise(samplenum = samplenum,
            intercept = m$coefficients[1], slope = m$coefficients[2], r2 = summary(m)$adj.r.squared, rse = summary(m)$sigma) ->
  slopes


rawdata_samples_trunc %>% 
  group_by(samplenum) %>% 
  summarise(DATETIME = mean(DATETIME), Sample = unique(Sample)) %>% 
  left_join(slopes, by = "samplenum") ->
  summarised_data

# Diagnostic plots

qc1 <- qplot(DATETIME, slope, data = summarised_data, geom = "line") + 
  facet_wrap(~Sample, scales = "free") +
  theme(axis.text.x = element_text(size = 6))
save_plot("qc1-slopes", qc1)
qc2 <- qplot(DATETIME, r2, data = summarised_data, geom = "line") + 
  facet_wrap(~Sample, scales = "free") +
  theme(axis.text.x = element_text(size = 6))
save_plot("qc1-r2", qc2)


save_data(summarised_data)

printlog("All done.")
closelog()
