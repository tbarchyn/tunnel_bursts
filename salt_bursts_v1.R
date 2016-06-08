# salt_bursts: scratch sheet to examine intermittent saltation near threshold
# Copyright 2016 Tom Barchyn
# 08 June 2016
# this program has no warranty whatsoever!

# idea:
# at a point responses of saltation are possibly the result of advecting bursts down
# the tunnel. These are not technically streamers, as the frequency of turbulence in
# the tunnel is much higher and eddy sizes sharply constrained. These are natural
# consequence of instability in fluid/impact threshold (not related to spatially
# and temporally autocorrelated forcing). This mechanism, alone, without reference
# to larger eddies, should produce a characteristic frequency-magnitude distribution.
# the idea is to look into this.

# This is a scratch sheet - meant for experimentation and understanding - lots of to dos

# This could be done analytically, done here numerically to set up for more complicated
# montecarlo ideas (which must be done numerically)


library (lattice)
library (plotly)
library (sfsmisc)

############################################################
# GLOBAL PARAMETERS  
crosswind_max <- 1.5                    # max crosswind seed distance
crosswind_granularity <- 0.01           # granularity of seeds
upwind_max <- 5.0                       # max upwind seed distance
upwind_granularity <- 0.01              # granularity of seeds

# saltation downwind ramp up
# saltation ramps up to some equilibrium magnitude, I don't have a good function
# here for this ramp up - so here is a quick approximation (replace with better
# version in time)

max_salt_magnitude <- 1.0               # some number for maximum saltation
ramp_up_distance <- 10.0                # the distance to reach maximum saltation
ramp_up_exponent <- 2.0                 # exponent on saltation ramp up

# burst cross-wind function
# the cross-wind diffusion is modeled as a gaussian curve with a width that
# increases linearly with downwind distance - replace with better version!
st_dev_width_mult <- 0.05               # this is multiplied by distance

# measurement threshold
# it is clear there is some measurement threshold in all instruments - this is minimum
# value that can be measured
measurement_threshold <- 0.001          # this end measurements are removed
loglog_breaks <- 20                     # the number of loglog breaks

############################################################
# FUNCTIONS
get_centerline_mag <- function (upwind) {
    # function to return the centerline magnitude of a burst as a function of
    # upwind distance to burst seed. Can pass an array of upwind seeds.
    # upwind = the upwind distance to burst seed
    
    mag <- upwind * NA
    mag [upwind >= ramp_up_distance] <- max_salt_magnitude
    mag [upwind < ramp_up_distance] <- max_salt_magnitude - 
        (max_salt_magnitude * ((ramp_up_distance - upwind[upwind < ramp_up_distance]) /
                                    ramp_up_distance)^ramp_up_exponent)
    return (mag)
}

get_mag <- function (upwind, crosswind) {
    # function to return the recorded magnitude as a function of seed location
    # upwind = the upwind location of the seed relative to the receptor
    # crosswind = the crosswind location of the seed relative to receptor
    
    max_mag <- get_centerline_mag (upwind)            # get the centerline magnitude
    st_dev <- upwind * st_dev_width_mult              # get the standard deviation
    zero_dens <- dnorm (0.0, mean = 0.0, sd = st_dev) # get a zero density
    
    # calculate the magnitude as fraction of centerline 'max_mag'
    mag <- max_mag * (dnorm (crosswind, mean = 0.0, sd = st_dev) / zero_dens)
    return (mag)
}

visualize_burst <- function () {
    # function to make a visualization of a canonical burst with specified parameters
    # this just hooks the bins specified above for seed locations as approximate
    # dimensions of interest - can be modified here if desired.
    
    crosswind_bins <- seq (0, crosswind_max, crosswind_granularity)
    upwind_bins <- seq (0, upwind_max, upwind_granularity)
    receptors <- matrix (NA, nrow = length (upwind_bins), ncol = length (crosswind_bins))
    
    # burst seed is at 0.0, 0.0
    for (i in 1:nrow(receptors)) {
        for (j in 1:ncol(receptors)) {
            receptors[i, j] <- get_mag (upwind_bins[i], crosswind_bins[j])
        }
    }
    #rownames (receptors) <- upwind_bins
    #colnames (receptors) <- crosswind_bins
    
    return (receptors)
}

get_burst_dist <- function () {
    # function to get distribution of bursts, based on the specified matrix of seed points
    # returns a frequency output from 'hist'
    # This assumes equal probability of a seed emminating from every spot in the grid.
    crosswind_bins <- seq (0, crosswind_max, crosswind_granularity)
    upwind_bins <- seq (0, upwind_max, upwind_granularity)
    mags <- matrix (NA, nrow = length (upwind_bins), ncol = length (crosswind_bins))
    
    # assign the magnitudes
    for (i in 1:nrow(mags)) {
        for (j in 1:ncol(mags)) {
            mags[i, j] <- get_mag (upwind_bins[i], crosswind_bins[j])
        }
    }
    
    # extract magnitude and frequency from bins
    mags_vec <- as.vector (mags)                              # convert to vector
    summary (mags_vec)                                        # throw summary
    mags_vec <- mags_vec [mags_vec > measurement_threshold]   # remove zeros, which aren't measured
    freq <- hist (mags_vec, breaks = 20, plot = F)            # create bins 
    return (freq)    
}

############################################################
# MAIN
# first, let's visualize the concentration behaviour of the burst (note that axes labels
# depend on the granularities both being 0.01 m!)
receptors <- visualize_burst ()
levelplot (receptors, col.regions = colorRampPalette(c('Blue', 'Red'))(100),
            main = 'canonical burst fluxes', xlab = 'downwind (cm)', ylab = 'crosswind (cm)')

# run some samples with various upwind max distances
upwind_max <- 2.5
fetch_2.5m <- get_burst_dist ()
upwind_max <- 5.0
fetch_5m <- get_burst_dist ()
upwind_max <- 10.0
fetch_10m <- get_burst_dist ()

# plot log log frequency magnitude
plot (fetch_2.5m$mids, fetch_2.5m$density, cex = 0, log = 'xy', main = 'frequency vs. burst flux',
      col = 'red', xlab = 'burst flux', ylab = 'density', xaxt = 'n', yaxt = 'n',
      xlim = c(10^-2, 10^0), ylim = c(10^-1, 10^2))
eaxis (1, at = c(10^-2, 10^-1, 10^0, 10^1, 10^2, 10^3, 10^4, 10^5))
eaxis (2, at = c(10^-2, 10^-1, 10^0, 10^1, 10^2, 10^3, 10^4, 10^5))
abline (h = c(10^-2, 10^-1, 10^0, 10^1, 10^2, 10^3, 10^4, 10^5), col = 'grey')
abline (v = c(10^-2, 10^-1, 10^0, 10^1, 10^2, 10^3, 10^4, 10^5), col = 'grey')

# add the lines and points
points (fetch_2.5m$mids, fetch_2.5m$density, col = 'red')
lines (fetch_2.5m$mids, fetch_2.5m$density, col = 'red')
points (fetch_5m$mids, fetch_5m$density, col = 'blue')
lines (fetch_5m$mids, fetch_5m$density, col = 'blue')
points (fetch_10m$mids, fetch_10m$density, col = 'forestgreen')
lines (fetch_10m$mids, fetch_10m$density, col = 'forestgreen')

legend (x = 0.3, y = 100,  legend = c('2.5 m', '5.0 m', '10.0 m'),
            col = c('red', 'blue', 'forestgreen'), lty = 1, title = 'fetch (m)')



