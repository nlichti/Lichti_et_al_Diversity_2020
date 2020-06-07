# -----------------------------------------------------------------------------
# Optimization code for:
# Lichti et al. (submitted) Shade, cache-pilferage, and anti-predator behavior 
#   in foragers may drive seed trait evolution in scatter-hoarded plants. 
#   Diversity X:XX-XX.
#
# Copyright (C) 2020 Nathanael Lichti (nlichti@purdue.edu) 
# Edited 06/06/2020
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

source('hoarding_model_06062020-1600.R')
# See note in hoarding_model_06062020-1600.R header re: unit conversions.

# -----------------------------------------------------------------------------
# Some function definitions used in this analysis

# Foraging optimization under given conditions
objective <- function(pars, state, params){
  
  # separate behavioral (theta) and non-behavioral state variables
  theta <- grab(c('rate','accept','use','attention','vigilance'), state)
  if(theta['vigilance.h1'] <= 0) return(1e12) 
  not_theta <- grab(names(theta), state, invert = TRUE)
  
  # sub in values from optimizer
  replace_theta <- names(pars)[names(pars) %in% names(theta)]
  replace_not_theta <- names(pars)[names(pars) %in% names(not_theta)]
  theta[replace_theta] <- pars[replace_theta]
  not_theta[replace_not_theta] <- pars[replace_not_theta]
  
  # gradient of perceived benefit with regard to behavioral variables
  v <- netBenefit(theta, not_theta = not_theta, args = params)
  if(is.na(v)) 1e12 else(-v)
}

reset <- function(x, state, params){
  if(length(state) > 0){
    x$state[names(state)] <- state
  }
  if(length(params) > 0){
    for(i in names(params)) x$parameters[[i]] <- params[[i]]
  }
  x
}

f = function(density, risk, energy, mode){
  r2$state['density.sp1.ownCache.h1'] <- ifelse(mode == 'pilfer', 0, density)
  r2$state['density.sp1.otherCache.h1'] <- ifelse(mode == 'retrieval', 0, 5.9 * density)
  r2$parameters$habitats[1,1] <- risk
  r2$parameters$fixed_attributes[1,1] <- energy
  r2$parameters$handling_time[1,1] <- r2$parameters$handling_time[1,1] * energy/60
  pars <- r2$state[7:12]
  with(r2, optim(par = pars, fn = objective, state = state, params = parameters, 
                 method='L-BFGS-B', lower=c(-12,1e-5,-12,-12,-12,-12), upper=c(12,Inf, 12,12,12,12))
             )
}
ff <- function(e, d){
  x <- optimize(peakdiff, c(3.5, 6), energy = e, density = d)
  x
}
peakdiff <- function(risk, energy, density){
  f_retrieve <- f(density = density, risk = risk, energy = energy, mode = 'retrieval')
  f_pilfer <- f(density = density, risk = risk, energy = energy, mode = 'pilfer')
  f_retrieve$value -  f_pilfer$value # note: these values are negative
}
getBenefit <- function(x){-x$value}
getPar <- function(x, which){x$par[which]}
findGUD <- function(e, d, m, x, threshold = 5000){
  x <- filter(x, energy == e, density == d & mode == m) 
  with(x, approx(net_utility, risk, xout = threshold)$y)
}
forageRate <- function(i, density, risk, energy, mode, fits, state, args){
  theta <- fits[[i]]$par
  not_theta <- state[!names(state) %in% names(theta)]
  not_theta['density.sp1.ownCache.h1'] <- density[i] * as.numeric(mode[i] != 'pilfer')
  not_theta['density.sp1.otherCache.h1'] <- density[i] * 5.9 * as.numeric(mode[i] != 'retrieval')
  args$habitats[1,1] <- risk[i]
  args$fixed_attributes[1,1] <- energy[i]
  args$handling_time[1,1] <- r2$parameters$handling_time[1,1] * energy[i]/60
  args <- parseInputs(theta, not_theta, args)   
  args <- resourceUse(args)
  args <- attributeUse(args)
  sum(args$f)/24 # seeds per hour
}


# -----------------------------------------------------------------------------
# Read basic parameterization files.  
# Parameters are defined in a set of Excel spreadsheets, along with initial 
# values for a dynamic version of the model that is not used here.

r1 <- setup('parameters.xlsx', 'initial_values.xlsx')

r2 <- reset(
  r1, 
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
      },
    recache_factor = 0.30
  )
)


# -----------------------------------------------------------------------------
# Accumulation of optimization results - this code can be rerun multiple times
# to add new evaluation points until a desired resolution is achieved in graphs.

# Container for all model optimization results
allFits <- list()

# Grid points for the current run
settings_now <- expand.grid(
  risk = seq(0.2,8,0.1),        
  density = c(200,1200,3200),  
  energy = c(5, 35, 65), 
  mode = c('retrieval','pilfer','choice')
)

# Save the current run's setting 
allFits$settings <- rbind(allFits$settings, settings_now)

# Execute the current run and append the results to the container slots
system.time({
  allFits$fits <- append(
    allFits$fits,
    with(
      settings_now, 
      mapply(f, density = density, risk = risk, energy = energy, mode = mode, SIMPLIFY = FALSE)
      ) 
  )
})

results <- with(
  allFits, 
  mutate(
    settings,
    net_utility = sapply(fits, getBenefit),
    attention_pilfer = sapply(fits, getPar, which='attention.pilfer'),
    foraging_rate = sapply(row_number(), forageRate, density, risk, energy, mode, 
                          fits, r2$state, r2$parameters)
  ) 
)
save(results, file = 'optimization_results.rda')
save(allFits, file = 'optimization_fits.rda')

# Utility (x energy)
epanel <- function(x) paste(x, 'kJ/seed')
ggplot(
  filter(results, mode != 'choice', energy %in% c(5, 35, 65), density %in% c(200, 1200, 3200)), 
  aes(x = risk, y = net_utility/86400, col = mode, lty=factor(density))
  ) + 
  geom_line() + 
  geom_hline(yintercept=5000/86400, col='black') +
  facet_wrap(.~energy, nrow=3, labeller = labeller(energy = epanel)) +
  labs(x = 'Risk coefficient', y='Perceived utility rate per second') +
  scale_color_manual(name = 'Strategy', values = c('steelblue','darkred')) +
  scale_linetype_manual(name = expression(Caches/m^2), values = c(4,5,1), labels = c(200,1200,3200)/10000) +
  theme_bw(base_size=18)


# -----------------------------------------------------------------------------
# Figure 1: Foraging rates by risk, density, seed value, and strategy
dpanel <- function(x)  paste(format(as.numeric(x)/10000, digits=2), 'caches per square meter')
annotations <- data.frame(
  risk = c(3.26, 6.075), 
  utility = 9, 
  arrow_x = c(3.25, 6.05), 
  arrow_y0 = 7.5, arrow_y1 = 4.8,
  lab = c('A','B'),
  density = 200
)
plot1 <- ggplot(
  filter(results, mode != 'choice', energy %in% c(5, 35, 65), density %in% c(200, 1200, 3200)), 
  aes(x = risk, y = net_utility/1440, col = mode, lty=factor(energy))
) + 
  geom_line() + 
  geom_hline(yintercept=5000/1440, col='black') +
  facet_wrap(.~density, nrow=3, labeller = labeller(density = dpanel)) +
  labs(x = 'Risk coefficient', y='Perceived utility rate per minute\n') +
  scale_color_manual(name = 'Strategy', values = c('steelblue','darkred')) +
  scale_linetype_manual(name = expression(kJ/seed), values = c(4,5,1), labels = c(5,35,65)) +
  geom_text(mapping = aes(x = risk, y = utility, label = lab, size = 1.2),
            data = annotations, show.legend = FALSE, inherit.aes = FALSE) +
  geom_segment(mapping = aes(x = arrow_x, y = arrow_y0, xend=arrow_x, yend = arrow_y1),
               arrow=arrow(length = unit(0.03,'npc'), type = 'closed'),
               data = annotations, show.legend = FALSE, inherit.aes = FALSE) + 
  theme_bw(base_size=18)
ggsave(plot1, file = 'figure1_utility_rates.png')
# shell.exec('figure1_utility_rates.png')  # windows only

# Foraging rates on an item per unit time basis
ggplot(
  filter(results, mode != 'choice'), 
  aes(x = risk, y = foraging_rate, col = mode, lty=factor(density))
) + 
  geom_line() + 
  facet_wrap(.~energy, nrow=2) +
  theme_bw() + labs(x = 'Risk coefficient', y='Foraging rate (seeds/hour)')


# -----------------------------------------------------------------------------
# Figure 2: Cache owner's advantage as a function of seed value, cache density, and risk
contrast <- filter(results, mode=='choice') %>% 
  mutate(advantage = filter(results, mode == 'retrieval')$net_utility - filter(results, mode == 'pilfer')$net_utility)
max_contrast <- contrast %>% group_by(energy, density) %>% 
  summarize() %>% ungroup() 
max_contrast2 <- filter(max_contrast, density %in% c(200,1200,3200)) %>%
  with({x <- mapply(ff, e = energy, d = density); as.data.frame(t(x))}) %>%
  bind_cols(filter(max_contrast, density %in% c(200,1200,3200))) %>% 
  select(energy, density, minimum, objective) %>%
  transmute( 
    energy = as.numeric(energy),
    density = as.numeric(density),
    risk = as.numeric(minimum),
    advantage = -as.numeric(objective)
  )
plot2 <- ggplot(
  filter(contrast, energy %in% c(5, 35, 65), density %in% c(200, 1200, 3200)) 
) + 
  # geom_line(data = max_contrast2, mapping = aes(x = risk, y = advantage/1440), color='steelblue') +
  geom_line(aes(x=risk, y=advantage/1440, lty=factor(energy)), col = 'steelblue') + 
  facet_wrap(.~density, nrow=3, labeller = labeller(density = dpanel)) +
  labs(x = 'Risk coefficient', y='Difference in rates (utility/min)\n') +
  scale_linetype_manual(name = expression(kJ/seed), values = c(4,5,1), labels = c(5,35,65)) +
  theme_bw(base_size=18)
ggsave(plot2, file = 'figure2_owner_advantage.png')
# shell.exec('figure2_owner_advantage.png')  # Windows only


# -----------------------------------------------------------------------------
# Figure 3: Risk-Value phase plane diagram

energy_risk <- results %>% filter(mode != 'choice') %>%
  group_by(energy, density, mode) %>% summarize() %>% ungroup() 
energy_risk <- bind_cols(
  energy_risk,
    risk =  with(energy_risk, {
      mapply(findGUD, e = energy, d = density, m = mode, 
             MoreArgs = list(x = results, threshold = 5000))
    })
  )   
peak = contrast %>% group_by(energy) %>% summarize(
  risk = risk[which.max(advantage)]
)

ribbon <- function(upper, lower){
  data.frame(
    x = c(sort(upper$risk), sort(lower$risk)),
    yupper = c(upper$energy, rep(max(upper$energy), nrow(lower))),
    ylower = c(rep(min(lower$energy), nrow(upper)), upper$energy)
  )
}

annotations <- data.frame(
  x = c(3.85, 5.4, 7.1, 4.52, 6.87), 
  y = c(55, 35, 15, 73, 73), 
  lab = c(' Pilferage\n competition',' No pilferage\n competition', 'Not used','A','B')
)

plot3 <- ggplot(filter(energy_risk, density == 1200)) + 
  geom_ribbon(data = ribbon(
    filter(energy_risk, density == 1200, mode=='pilfer'),
    filter(energy_risk, density == 1200, mode=='retrieval')
  ), mapping = aes(x=x, ymax=yupper, ymin=ylower), 
    fill = 'grey', alpha = 0.2) +
  # geom_line( # This will trace the optima w.r.t risk across energy value
  #   data = filter(max_contrast2, density==1200), 
  #   mapping = aes(x = risk, y = energy), col='grey'
  #   ) +
  geom_line(aes(x = risk, y = energy, group = mode, col=mode)) +
  geom_text(data = annotations, mapping = aes(x=x, y=y, label = lab),
            inherit.aes = FALSE, size = 6) +
  scale_color_manual(name = 'Strategy', values = c('steelblue','darkred'), guide=FALSE) +
  scale_x_continuous(name = '\n Risk coefficient', breaks = NULL, limits = c(3.5,7.5)) +
  scale_y_continuous(name = 'Seed value\n', breaks = NULL) +
  theme_bw(base_size=18)
ggsave(plot3, file = 'figure3_risk-value_phase_plane.png')
# shell.exec('figure3_risk-value_phase_plane.png')  # Windows only

# -----------------------------------------------------------------------------
# Figure 4 - net benefit with establishment costs
# proportionality constants u_travel and u_fear were chosen 
# by trial and error for illustrative purposes

radius = function(n, density) sqrt((n/(density/10000))/pi)
weighted_radius <- function(x, k) 2 * pi * x^2 / k
mean_travel <- function(n, density){
  r <- radius(n, density)
  k <- 2 * pi * integrate(I, 0, r)$value
  integrate(weighted_radius, 0, r, k = k)$value
}

peakbenefit <- function(density, risk, energy, n, travel_cost, fear_cost=1){
  density <- exp(density)
  f_retrieve <- f(density = density, risk = risk, energy = energy, mode = 'retrieval')
  f_pilfer <- f(density = density, risk = risk, energy = energy, mode = 'pilfer')
  z <- mean_travel(n, density)
  establishment_cost <- z * 2 * travel_cost +
                        (fear_cost * (z*2/3.5 +  40) * 10^risk / 86400)  #seconds/minute
  (f_retrieve$value -  f_pilfer$value)/86400 + establishment_cost # note: these values are negative
}

delta <- NULL  # in case it needs to be rerun to get more points
risks <- c(seq(1, 4.6, 0.05), seq(4.65, 4.79, 0.01), seq(4.8, 6, 0.05))
pargrid <- expand.grid(energy=c(5,35,65), risk = risks)
system.time({ # this takes a while - parallelize it!
  op <- with(pargrid, mapply(function(e,r){
  x <- optimize(peakbenefit, c(-16,16), energy = e, risk = r, n = 500, 
                travel_cost = 0.0002, fear_cost = 0.0002)
  as.data.frame(x)
}, e = energy, r = risk))
})
op <- as.data.frame(t(op))
delta <- bind_cols(
  pargrid,
  mutate(
    op, 
    density = exp(as.numeric(minimum)),
    benefit=-(as.numeric(objective) - as.numeric(objective)[1:3])
  )
) %>% bind_rows(delta)
optima <- sapply(c(5,35,65), function(i){
  with(
    filter(delta, energy==i), 
    c(risk = risk[which.max(benefit)], benefit = max(benefit))
  )}) %>% t() %>% as.data.frame()

plot4inset <- ggplot(filter(delta, risk > 4.5, risk < 4.8)) +
  geom_line(aes(x=risk, y=100^benefit, linetype=factor(energy)), color = 'steelblue') +
  geom_point(data=optima, aes(x = risk, y = 100^benefit), size = 3, col = 'steelblue') +
  scale_linetype_manual(values = c(4,5,1), guide=FALSE) +
  scale_x_continuous(name = NULL) +
  scale_y_continuous(name = NULL, breaks = c(1.200, 1.205, 1.210, 1.215), 
                     labels = round(log(c(1.200, 1.205, 1.210, 1.215),100), 4)) +
  theme_bw(base_size=12)
plot4 <- ggplot(filter(delta, risk < 5.75), aes(x=risk, y=benefit, linetype=factor(energy))) +
  geom_line(color = 'steelblue') +
  annotation_custom(
    ggplotGrob(plot4inset),
    xmin = 0.8, xmax = 3.2, ymin = 0.02, ymax = 0.045
  ) +
  scale_linetype_manual(name = expression(kJ/seed), values = c(4,5,1), labels = c(5,35,65)) +
  scale_x_continuous(name = expression(paste('\n Risk coefficient (', italic(C), ')'))) +
  scale_y_continuous(name = expression(paste('Change in net benefit from', italic(C)==1))) +
  theme_bw(base_size=18)
ggsave(plot4, file = 'figure4_benefit_gradient.png')
# shell.exec('figure4_benefit_gradient.png')  # windows only




