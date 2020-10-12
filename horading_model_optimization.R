# -----------------------------------------------------------------------------
# Optimization code for:
# Lichti et al. (submitted) Shade, cache-pilferage, and anti-predator behavior 
#   in foragers may drive seed trait evolution in scatter-hoarded plants. 
#   Diversity X:XX-XX.
#
# Copyright (C) 2020 Nathanael Lichti (nlichti@purdue.edu) 
# Edited 10/12/2020
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>
#
# -----------------------------------------------------------------------------

source('hoarding_model_07232020-1635.R')
# See note in hoarding_model_07232020-1635.R header re: unit conversions.

library(parallel)
library(optimx)
library(DEoptim)
library(readxl)
library(dplyr)
library(ggplot2)

# -----------------------------------------------------------------------------
# Some function definitions used in this analysis

# Foraging utility optimization under given conditions
objective <- function(pars, .state, .params, pnames = names(pars)){
  names(pars) <- pnames
 
  # separate behavioral (theta) and non-behavioral state variables
  theta <- grab(c('rate','accept','use','attention','vigilance'), .state)
  # if(theta['vigilance.h1'] <= 0) return(1e12) 
  if(any(abs(theta)>12)) return(1e12) 
  not_theta <- grab(names(theta), .state, invert = TRUE)
  
  # sub in values from optimizer
  replace_theta <- names(pars)[names(pars) %in% names(theta)]
  replace_not_theta <- names(pars)[names(pars) %in% names(not_theta)]
  theta[replace_theta] <- pars[replace_theta]
  not_theta[replace_not_theta] <- pars[replace_not_theta]
  
  # gradient of perceived benefit with regard to behavioral variables
  v <- netBenefit(theta, not_theta = not_theta, args = .params)
  -v 
}

# One-step replacement of selected parameter and state values
reset <- function(x, state, params){
  if(length(state) > 0){
    x$state[names(state)] <- state
  }
  if(length(params) > 0){
    for(i in names(params)) x$parameters[[i]] <- params[[i]]
  }
  x
}

# Function to execute optimization under one replicate of the modeling experiment  
# at specified values for cache density, predation risk, energy content of seeds,
# utiltiy gain of eating 1 kJ of energy (gain), utility gained by recaching 1 kJ 
# of energy (recache), and the foraging mode ('retrieval','pilfer',or 'choice').
# Set useDE = TRUE to force use of differential evolution (DEoptim) instead of 
# L-BFGS-B (DE is always used for mode = 'choice').  nDE scales the number
# of DE agents used (the actual number is nDE * length(par) where par is the
# parameter vector to be optimized).
f = function(density, risk, energy, gain, recache, mode, useDE=FALSE, nDE = 20){
  r1$state['gain.energy.eat'] <- gain
  r1$state['density.sp1.ownCache.h1'] <- ifelse(mode == 'pilfer', 0, density)
  r1$state['density.sp1.otherCache.h1'] <- ifelse(mode == 'retrieval', 0, 5.9 * density)
  if(mode != 'choice'){
    r1$state['attention.pilfer'] <- ifelse(mode == 'pilfer', 12, -12)
    r1$state['attention.retrieve'] <- ifelse(mode == 'retrieval', 12, -12)
  }
  r1$parameters$recache_factor <- recache #log10(recache)
  r1$parameters$habitats[1,1] <- 10^risk
  r1$parameters$fixed_attributes[1,1] <- energy
  r1$parameters$handling_time[1,1] <- r1$parameters$handling_time[1,1] * energy/60
  pars <- switch(
    which(mode == c('retrieval','pilfer','choice')),
    r1$state[c(8:9,11)],
    r1$state[c(7,9,11)],
    r1$state[c(7:9,11:13)]
  )
  ans <- with(r1, 
          if(mode != 'choice' && !useDE){
            aa <- optim(par = pars, fn = objective, .state = state, .params = parameters, method='L-BFGS-B',
                   lower=rep(-12, length(pars)), upper=rep(12, length(pars)))
            data.frame(value=aa$value, t(aa$par))
          } else {
            DEoptim::DEoptim(fn = objective, .state = state, .params = parameters, pnames = names(pars),
                   lower=rep(-12, length(pars)), upper=rep(12, length(pars)),
                   control = list(NP = length(pars)*nDE, reltol = 1e-4, trace = 0))
          }
  )
  if(!is.data.frame(ans)){
    ans <- with(ans$optim, setNames(c(bestval, bestmem), c('value', names(pars))))
    ans <- as.data.frame(t(ans))
  }
  if(!any(grepl('searchRate', names(ans)))) 
    ans['rate.h1.searchRate'] <- r1$state['rate.h1.searchRate']
  if(!any(grepl('relocationRate', names(ans)))) 
    ans['rate.h1.relocationRate'] <- r1$state['rate.h1.relocationRate']
  if(!any(grepl('attention.retrieve', names(ans)))) 
    ans['attention.retrieve'] <- NA
  if(!any(grepl('attention.pilfer', names(ans)))) 
    ans['attention.pilfer'] <- NA
  cbind(
    data.frame(mode = mode, risk = risk, density = density, energy = energy, gain = gain, recache = recache),
    ans
    )
}

# -----------------------------------------------------------------------------
# Read basic parameterization files.  
# Parameters are defined in a set of Excel spreadsheets, along with initial 
# values for a dynamic version of the model that is not used here.

r1 <- reset(
  setup('parameters.xlsx', 'initial_values.xlsx'), 
  state = c(
    density.sp1.free.h1 = 0,
    density.sp1.ownCache.h1 = 200,
    density.sp1.otherCache.h1 = 200 * 5.9,
    density.sp1.abandoned.h1 = 0,
    density.sp1.dead.h1 = 0,
    competition.scca.popden = 0,
    gain.energy.eat = 2,
    attention.retrieve = 0,
    attention.pilfer = 2
  ),
  params = list(
    states = {
      x = r1$parameters$states
      x[1,1] = '1.0'
      x
      }
  )
)

## TEST the optimization for the a given point
f(density = 1200, risk = 0, energy = 35, gain = 1, recache = 0.5, mode = 'pilfer')
f(density = 1200, risk = 0, energy = 35, gain = 1, recache = 0.5, mode = 'retrieval')
f(density = 1200, risk = 0, energy = 35, gain = 1, recache = 0.5, mode = 'choice', nDE = 20)

# -----------------------------------------------------------------------------
# Accumulation of optimization results - this code can be rerun multiple times
# to add new evaluation points until a desired resolution is achieved in graphs.

# Container for all model optimization results
allFits <- list()

# Grid points for the current run
settings_now <- bind_rows(
  # scenarios for Figures 1 and 2
  expand.grid(
    mode = c('retrieval','pilfer','choice'),
    risk = seq(-2, 3, 0.02),
    density = c(200,1200,3200),  
    energy = c(5, 35, 65), 
    gain = c(0.2, 0.5, 1, 2, 5), 
    recache = 1/10 
  ),
  # parameters for figure 3
  expand.grid(
    density = 1200, gain = 5, risk = seq(0,1.5,0.1),
    energy = seq(5, 65, 2.5)[!seq(5, 65, 2.5) %in% c(5,35,65)],
    mode = c('pilfer','retrieval'), recache = 0.1
  )
)

# Save the current run's setting 
allFits$settings <- rbind(allFits$settings, settings_now)

# Execute the current run and append the results to the container slots
# Set up cluster and export necessary funcions
ncore <- detectCores()
clstr <- makeCluster(ncore)
clusterEvalQ(clstr, source('hoarding_model_07232020-1635.R'))
clusterEvalQ(clstr, check_packages('DEoptim'))
clusterExport(clstr, c('f','objective','r1'))

# At least 2 optimization passes may be needed. The L-BFGS-B algorithm is fast 
# but fragile and sensitive to initial parameter values.  Failures in this problem
# occur sporatically and on the lower slopes of the curves for pilferage. They
# are obvious in plots of the results and appear as radical deviations from the 
# overall curve.  Rerunning the problematic points using the DE algorithm (slow
# but reliable for difficult problems) will obtain reasonable results.
system.time({
  allFits$fits <- append(
    allFits$fits,
    with(
      settings_now, 
      clusterMap(clstr, f, density = density, risk = risk, energy = energy, 
                 gain = gain, recache = recache, mode = mode, 
                 useDE=TRUE, SIMPLIFY = FALSE, nDE = 25, 
                 # set useDE to force use of differential evolution
                 .scheduling = 'dynamic')
    ) 
    # For serial processing (slow), use:
    # mapply(f, density = density, risk = risk, energy = energy, gain = gain, 
    #        recache = recache, mode = mode, SIMPLIFY = FALSE)
    # ) 
  )
})

stopCluster(clstr)
save(allFits, file = 'optimization_fits.rda')

results <- do.call(bind_rows, allFits$fits) %>%
  mutate(
    mode = mode,
    risk = round(risk, 2),
    density = density,
    energy = energy,
    gain = gain,
    recache = recache,   
    net_utility = -value,
    pilfer = exp(attention.pilfer)/(exp(0) + exp(attention.pilfer) + exp(attention.retrieve)),
    vigilance = plogis(vigilance.h1),
    cache = plogis(use.sp1.cache),
    search = plogis(rate.h1.searchRate) * r1$parameters$max_search_rate,
    relocation = plogis(rate.h1.relocationRate) * r1$parameters$max_relocation_rate
  ) %>% 
  group_by(mode, risk, density, energy, gain, recache) %>%
  summarize(
    pilfer = pilfer[which.max(net_utility)],
    vigilance = vigilance[which.max(net_utility)],
    cache = cache[which.max(net_utility)],
    search = search[which.max(net_utility)],
    relocation = relocation[which.max(net_utility)],
    net_utility = max(net_utility)
  ) %>% ungroup()

save(results, file = 'optimization_results.rda')

# Calculate contrasts for Figure 2
contrast <- results %>% 
  group_by(risk, density, energy, gain, recache) %>%
  summarize(
    dBenefit = net_utility[mode == 'retrieval'] - net_utility[mode == 'pilfer']
  ) %>% ungroup()

# -----------------------------------------------------------------------------
# Figure 1: Foraging rates by risk, density, seed value, and strategy

# Labelers for facet panels
# Density
dpanel <- function(x){
  x <- lapply(x, function(i) format(as.numeric(i)/10000, digits=2))
  lapply(unname(x), lapply, function(values) {
    values <- paste0('paste(', values, ', ~caches/m^2)')
    c(parse(text = as.character(values)))
  })
} 

# Gain (forager state)
gpanel <- function(x){
  lapply(unname(x), lapply, function(values) {
    values <- paste0('paste(italic(g)[eat] == ', values, ')')
    c(parse(text = as.character(values)))
  })
} 

class(dpanel) <- class(gpanel) <- append(class(dpanel), 'labeller')

# Set up plot
plot1 <- ggplot(
  filter(results, mode != 'choice', recache == 0.1, energy %in% c(5, 35, 65), 
         density %in% c(200, 1200, 3200), gain %in%c(0.2,1,5)), 
  aes(x = risk, y = net_utility, col = mode, lty=factor(energy))
) + 
  geom_hline(yintercept=0.05, col='black') +
  geom_line() + 
  facet_grid(density~gain, labeller = labeller(density = dpanel, gain=gpanel)) +
  labs(
    x = expression(paste('Risk coefficient (',log[10],italic(C),')')), 
    y='Perceived benefit (utility/min)\n'
    ) +
  scale_color_manual(name = 'Strategy', limits = c('retrieval','pilfer'), values = c('steelblue','darkred')) +
  scale_linetype_manual(name = expression(kJ/seed), values = c(4,5,1), labels = c(5,35,65)) +
  theme_bw(base_size=18)

# save and view
ggsave(plot1, file = 'figure1_utility_rates.png')
shell.exec('figure1_utility_rates.png')  # windows only

# -----------------------------------------------------------------------------
# Figure 2: Cache owner's advantage as a function of 
#           seed value, cache density, and risk

# Function to find the risk value at which the max-benefit curve 
# crosses a given threshold value.  
findGUD <- function(x, y, threshold = 0.5){
  approx(x, y, xout = threshold)$y
}

# Locate arrow placement in Figure 2.  The arrowplace variable determines
# wheather the arrow is above or below the zero line.
GURisk <- results %>% group_by(mode, energy, density, gain, recache) %>%
  summarize(
    risk = findGUD(x = net_utility, y = risk, threshold = 0.05)
  ) %>%
  mutate(
    arrowplace = ifelse(gain > 1, -1, 1)
  )

# Set up Figure 2
plot2 <- ggplot(
  filter(contrast, recache == 0.1, gain %in% c(0.2,1,5), energy %in% c(5, 35, 65))#, density %in% c(200, 1200, 3200)) 
) + 
  geom_segment(
    aes(x = risk, xend = risk, y = arrowplace * 0.3, yend = arrowplace * 0.1, col = mode), 
    data = filter(GURisk, energy == 65, recache == 0.1, gain %in% c(0.2,1,5)),
    arrow = arrow(length = unit(0.15, 'cm'), type = 'closed')
    ) + scale_color_manual(name = 'Strategy', limits = c('retrieval','pilfer'), values=c('steelblue','darkred')) +
  geom_line(aes(x=risk, y=dBenefit, lty=factor(energy)), col = 'steelblue') + 
  geom_hline(yintercept = 0) +
  facet_grid(density~gain, labeller = labeller(density = dpanel, gain=gpanel)) +
  labs(
    x = expression(paste('Risk coefficient (',log[10],italic(C),')')), 
    y='Difference in perceived benefit (utility/min)\n'
    ) +
  scale_linetype_manual(name = expression(kJ/seed), values = c(4,5,1), labels = c(5,35,65)) +
  theme_bw(base_size=18)

# Save and view
ggsave(plot2, file = 'figure2_rate_contrast.png')
shell.exec('figure2_rate_contrast.png')  # Windows only

# -----------------------------------------------------------------------------
# Figure 3: Risk-Value phase plane diagram

# Detailed results for 1200 caches/ha, gain = 5, recache = 0.1
energy_risk <- results %>% filter(mode != 'choice') %>%
  group_by(mode, recache, gain, energy, density) %>% summarize(
    risk = findGUD(e=energy, d = density, m=mode, rc = recache, 
                   g = gain, x = net_utility, y = risk, 
                   threshold = 0.05)
  )   

# Function to find the boundaries of the shaded area in Figure 3
ribbon <- function(upper, lower){
  x = c(sort(upper$risk), sort(lower$risk))
  data.frame(
    x = x,
    yupper = ifelse(
      x < max(upper$risk), 
      approx(upper$risk, upper$energy, x)$y,
      max(upper$energy)
    ),
    ylower = ifelse(
      x > min(lower$risk, na.rm = TRUE), 
      approx(lower$risk, lower$energy, x)$y,
      min(lower$energy, na.rm=TRUE)
    )
  )
}

# Text annotations for Figure 3
annotations <- data.frame(
  x = c(0.25, 0.8, 1.3, 0.6, 1.55, 0.7), 
  y = c(55, 35, 15, 68, 68, 68), 
  col = c('black','black','black','red','steelblue','black'),
  lab = c(' Pilferage\n competition',' No pilferage\n competition', 'Not used','A','B','C')
)

# Functions to smooth the peak isopleth (C)
findPeak <- function(risk, diff){
  f <- splinefun(x=risk,y=diff)
  negf <- function(x) -f(x)
  optimize(negf, c(0,1))$minimum
}
findIntersect <- function(risk, diff){
  f <- splinefun(x=risk,y=diff^2)
  negf <- function(x) f(x)
  optimize(negf, c(-2,0.25))$minimum
}

# maximum contrast
peak_contrast <- contrast %>% group_by(energy, density, gain, recache) %>%
  summarize(
    intersect = findIntersect(risk, dBenefit),
    risk = findPeak(risk, dBenefit)
  ) %>% ungroup() %>% filter(gain == 5, density == 1200, recache == 0.1) %>%
  mutate(
    risk_smooth = smooth.spline(energy, risk, spar=0.6)$y
  )

# Set up Figure 3
plot3 <- ggplot(filter(energy_risk, density == 1200, gain==5, recache == 0.1)) + 
  geom_ribbon(data = ribbon(
    filter(energy_risk, density == 1200, gain==5, recache == 0.1, mode=='pilfer'),
    filter(energy_risk, density == 1200, gain==5, recache == 0.1, mode=='retrieval')
  ), mapping = aes(x=x, ymax=yupper, ymin=ylower), 
    fill = 'grey', alpha = 0.2) +
  
  geom_line(
    data = peak_contrast, 
    aes(x = intersect, y = energy)
  ) +
  
  geom_line(
    data = peak_contrast, 
    aes(x = risk_smooth, y = energy)
    ) +
  
  geom_line(aes(x = risk, y = energy, group = mode, col=mode)) +
  geom_text(data = annotations, mapping = aes(x=x, y=y, label = lab),
            inherit.aes = FALSE, size = 6) +
  scale_color_manual(name = 'Strategy', values = c('darkred','steelblue'), guide=FALSE) +
  scale_x_continuous(name = '\n Risk coefficient', breaks = NULL, limits = c(0.1,1.6)) +
  scale_y_continuous(name = 'Seed value\n', breaks = NULL) +
  theme_bw(base_size=18)

# Save and view
ggsave(plot3, file = 'figure3_risk-value_phase_plane.png')
shell.exec('figure3_risk-value_phase_plane.png')  # Windows only

#-------------------------------------------------------------------------------
# Supplement figure S2: retrieval with choice of strategy
plotS2 <- ggplot(
  filter(results, recache == 0.1, energy %in% c(5, 35, 65), 
         density %in% c(200, 1200, 3200), gain %in%c(0.2,1,5)), 
  aes(x = risk, y = net_utility, col = mode, alpha=mode, lwd=mode, lty=factor(energy))
) + 
   geom_hline(yintercept=0.05, col='black') +
   geom_line() + 
   facet_grid(density~gain, labeller = labeller(density = dpanel, gain=gpanel)) +
   labs(
     x = expression(paste('Risk coefficient (',log[10],italic(C),')')), 
     y='Perceived utility rate per minute\n'
   ) +
   scale_color_manual(name = 'Strategy', limits = c('retrieval','pilfer','choice'), values = c('steelblue','darkred','black')) +
   scale_alpha_manual(limits = c('retrieval','pilfer','choice'), values = c(0.3,0.3,1), guide=FALSE) + 
   scale_size_manual(limits = c('retrieval','pilfer','choice'), values = c(1.5,1.5,0.5), guide=FALSE) + 
   scale_linetype_manual(name = expression(kJ/seed), values = c(4,5,1), labels = c(5,35,65)) +
   theme_bw(base_size=18)
 ggsave(plotS2, file = 'figureS2_utility_with_choice.png')
 shell.exec('figureS2_utility_with_choice.png') 

 # Figure S3: Vigilance
 plotS3 <- ggplot(
   filter(results, mode!='choice',recache == 0.1, energy %in% c(5, 35, 65), #net_utility > 1e-6, 
          density %in% c(200, 1200, 3200), gain %in%c(0.2,1,5)), 
   aes(x = risk, y = vigilance, col = mode, lty=factor(energy))
 ) + 
   geom_line() + 
   facet_grid(density~gain, labeller = labeller(density = dpanel, gain=gpanel)) +
   labs(
     x = expression(paste('Risk coefficient (',log[10],italic(C),')')), 
     y='Proportion of attention to vigilance (v)\n'
   ) +
   scale_color_manual(name = 'Strategy', limits = c('retrieval','pilfer'), values = c('steelblue','darkred')) +
   scale_linetype_manual(name = expression(kJ/seed), values = c(4,5,1), labels = c(5,35,65)) +
   theme_bw(base_size=18)
 ggsave(plotS3, file = 'figureS4_vigilance.png')
 shell.exec('figureS3_vigilance.png') 
 
# Figure S4: Choice between caching and eating collected seeds
plotS4 <- ggplot(
  filter(results, mode!='choice',recache == 0.1, energy %in% c(5, 35, 65), net_utility > 1e-6, 
         density %in% c(200, 1200, 3200), gain %in%c(0.2,1,5)), 
  aes(x = risk, y = cache, col = mode, lty=factor(energy))
) + 
  geom_line() + 
  facet_grid(density~gain, labeller = labeller(density = dpanel, gain=gpanel)) +
  labs(
    x = expression(paste('Risk coefficient (',log[10],italic(C),')')), 
    y='Probability that collected seeds will be cached\n'
  ) +
  scale_color_manual(name = 'Strategy', limits = c('retrieval','pilfer','choice'), values = c('steelblue','darkred','black')) +
  scale_linetype_manual(name = expression(kJ/seed), values = c(4,5,1), labels = c(5,35,65)) +
  theme_bw(base_size=18)
 ggsave(plotS3, file = 'figureS4_seed_allocation.png')
 shell.exec('figureS4_seed_allocation.png') 
 






