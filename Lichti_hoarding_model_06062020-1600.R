# -----------------------------------------------------------------------------
# Model code and associated utility functions for:
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
#
# NOTE: This code has been borrowed from another project that aims to model 
# seed-hoarder-predator interactions over time using a differential equation 
# model.  It is designed to track populations on a per hectare, per day basis,
# but many of the parameters are more intuitively defined on a smaller scale
# of square meters per second or minute. Unit conversions are included 
# in the code.
#
# -----------------------------------------------------------------------------
# General utilities
check_packages <- Vectorize(function(x) {
  if (!require(x, character.only = TRUE)) {
    message(paste('Installing required package:', x))
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

mapvalues <- function (x, from, to, warn_missing = TRUE) {
  # this is direct from the plyr package.
  if (length(from) != length(to)) {
    stop("`from` and `to` vectors are not the same length.")
  }
  if (!is.atomic(x)) {
    stop("`x` must be an atomic vector.")
  }
  if (is.factor(x)) {
    levels(x) <- mapvalues(levels(x), from, to, warn_missing)
    return(x)
  }
  mapidx <- match(x, from)
  mapidxNA <- is.na(mapidx)
  from_found <- sort(unique(mapidx))
  if (warn_missing && length(from_found) != length(from)) {
    message("The following `from` values were not present in `x`: ", 
            paste(from[!(1:length(from) %in% from_found)], collapse = ", "))
  }
  x[!mapidxNA] <- to[mapidx[!mapidxNA]]
  x
}

vec_eval <- function(x){
  dims <- dim(x)
  x <- paste(x, collapse = ',')
  x <- paste('c(', x, ')', collapse = '')
  x <- eval.parent(parse(text = x))
  dim(x) <- dims
  x
}

vgsub <- function(pattern, replacement, x, ...){
  for(i in seq_along(pattern)){
    x <- gsub(pattern[i], replacement[i], x, ...)
  }
  x
}

grab <- function(what, from, invert = FALSE, ret.names = FALSE, ...){
  if(is.null(names(from)) && is.character(from)){
    names(from) <- from
    ret.names <- TRUE
  } 
  if(ret.names) from <- setNames(names(from), names(from))
  what <- sapply(what, grepl, x=names(from), ...)
  what <- if(invert){ 
    rowSums(!what) == ncol(what) 
  } else {
    rowSums(what) > 0
  }
  from[what]
}

sumBy <- function(what, from, invert = FALSE, ...){
  values <- lapply(what, grab, from = from, invert = invert, ...)
  setNames(sapply(values, sum), what)
} 

# -----------------------------------------------------------------------------
# utilities specific to this problem
loadSheet <- function(file, sheet, n=0){
  check_packages('readxl')
  x <- as.data.frame(read_excel(file, sheet = sheet, skip = 1))
  if(n > 0){
    rn <- apply(x[,1:n,drop = FALSE], 1, paste, collapse = '.')
    rownames(x) <- rn
    x <- x[,-(1:n), drop = FALSE]
  } 
  x <- as.matrix(x)
  if(is.character(x)){
    x[] <- gsub(' ', '', x)
  }
  x
}

trimNames <- function(x){
  x <- strsplit(x, '.', fixed=TRUE)[[1]]
  if(length(x) == 1) return(x)
  if(x[1] == x[length(x)]) x = x[-length(x)]
  paste(x, collapse = '.')
}

flattenMatrix <- function(x){
  x <- as.matrix(x)
  atrs <- attributes(x)
  if(is.null(rownames(x))){
    pointnames <- colnames(x)
  } else {
    xnames <- expand.grid(rownames(x), colnames(x))
    pointnames <- apply(xnames, 1, paste, collapse = '.')
  }
  x <- c(x)
  attributes(x) <- list(folded = atrs)
  x <- setNames(x, pointnames)
  x[!is.na(x)]
} 

formatRate <- function(x, linearize = TRUE, max = 1, tol = 1e-5){
  if(linearize){
    y <- x[,1]
    y[y > 1 - tol] <- 1 - tol
    y[y < tol] <- tol
    x[,1] <- qlogis(y/max)
    x
  } else {
    plogis(x) * max
  }
}

formatUse <- function(x, linearize = TRUE, dims = NULL, names = NULL){
  if(linearize){
    x <- (log(x) - log(x[,1]))[,-1, drop=FALSE]
    colnames(x) <- paste(colnames(x), 'use', sep='.')
    x
  } else {
    x <- matrix(x, dims[1], dims[2] - 1)
    x <- cbind(0, x)
    dimnames(x) <- names
    exp(x)/rowSums(exp(x))
  }
}

formatAttention <- function(x, linearize = TRUE, names = NULL){
  if(linearize){
    (mgcv::notLog(x) - mgcv::notLog(x[1,]))[-1,,drop = FALSE]
  } else {
    x <- setNames(c(0, x), names)
    mgcv::notExp(x)/sum(mgcv::notExp(x))
  }
}

normalizeAttention <- function(A, args){
  strategies <- args$states[,'strategy']
  densities <- args$n
  nii <- numeric(0)
  for(i in unique(strategies)){
    ni <- grab(names(strategies[strategies == i]), densities)
    if(sum(ni)==0){
      ni[] <- 1/length(ni)
    } else {
      ni <- ni/sum(ni)
    }
    nii <- c(nii, ni)
  }
  A * nii[names(A)] 
}

formatAcceptance <- function(x, linearize = TRUE, tol = 1e-5){
  if(linearize){
    y <- x[,1]
    y[y > 1 - tol] <- 1 - tol
    y[y < tol] <- tol
    x[,1] <- qlogis(y)
    x
  } else {
    plogis(x)
  }
}

fixGains <- function(inits, params, warn = TRUE){
  params$fixed_gains <- with(params, {
    is_fixed <- fixed_gains
    gains <- inits$gain
    toxins <- toxins[, 'cache', drop = FALSE]
    rn <- rownames(toxins)
    cn <- colnames(toxins)
    check <- all(gains[rn,cn] == toxins) 
    if(!check){
      if(warn) 
        message(paste(
          'Initial values for toxins with use =', cn, 
          'overridden by parameter inputs.'
        ))
      gains[rn,cn] <- toxins
    }
    replace <- !colnames(gains) %in% cn
    v <- gains[rn, replace]
    gains[rn, replace] <- ifelse(v <= toxins, v, toxins)
    ifelse(is_fixed, gains, NA)
  })
  params
}

formatGain <- function(values, params, setup = FALSE){
  if(setup){
    values[!is.na(params$fixed_gains)] <- NA
    return(values)
  }
  G <- params$fixed_gains
  G[is.na(G)] <- values
  G
}

formatTransitions <- function(x){
  N <- 0
  v <- vec_eval(x)
  x[] <- ifelse(
    is.na(v),
    paste('ifelse(N==0,0,', x, ')', sep=''),
    x
  )
  x
}

fixTransitions <- function(params){
  params$eat_matrix <- formatTransitions(params$eat_matrix)
  params$cache_matrix <- formatTransitions(params$cache_matrix)
  params
}

matchParameters <- function(values, template, map = NULL){
  if(is.null(dim(values))){
    template_names <- names(template)
    if(!is.null(map)){
      names(template) <- vgsub(names(map), map, names(template))
    }
    n <- strsplit(names(values), '.', fixed = TRUE)
    if(length(n[[1]]) > 1){
      n <- lapply(n, '[', -1)
    }
    n <- sapply(n, paste, collapse = '.')
    for(i in n){
      to_replace <- grep(i, names(template))
      replace_with <- grep(i, names(values))
      template[to_replace] <- values[replace_with]
    }
    names(template) <- template_names
  } else {
    n <- rownames(values)
    template <- matrix(
      NA, ncol(values), length(template), 
      dimnames = list(
        colnames(values),
        names(template)
      )
    )
    for(i in n){
      to_replace <- grep(i, colnames(template))
      replace_with <- grep(i, rownames(values))
      template[,to_replace] <- values[replace_with,]
    }
    template <- t(template)
    if(ncol(template) == 1){
      template <- setNames(c(template), rownames(template))
    }
  }
  template
}

importParameters <- function(file){
  check_packages('readxl')
  list(
    fixed_attributes = loadSheet(file, sheet = 'resource_attributes', n = 1),
    states = loadSheet(file, sheet = 'states', n = 1),
    nutrients = loadSheet(file, sheet = 'nutrients', n = 1),
    toxins = loadSheet(file, sheet = 'toxins', n = 1),
    handling_time = loadSheet(file, sheet = 'handling_time', n = 1),
    general = loadSheet(file, sheet = 'general', n = 1),
    habitats = loadSheet(file, sheet = 'habitats', n = 1),
    gammas = loadSheet(file, sheet = 'gammas', n = 1),
    eat_matrix = loadSheet(file, sheet = 'eat_transform', n = 1),
    cache_matrix = loadSheet(file, sheet = 'cache_transform', n = 1),
    fixed_gains = loadSheet(file, sheet = 'fix_gain', n = 1)
  )
}

importStates <- function(file){
  list(
    density = loadSheet(file, sheet = 'density', n = 3),
    competition = loadSheet(file, sheet = 'competition', n = 1),
    rate = loadSheet(file, sheet = 'rate', n = 1),
    vigilance = loadSheet(file, sheet = 'vigilance', n = 1),
    accept = loadSheet(file, sheet = 'accept', n = 1),
    use = loadSheet(file, sheet = 'use', n = 1),
    attention = loadSheet(file, sheet = 'attention', n = 1),
    gain = loadSheet(file, sheet = 'gain', n = 1),
    resources = loadSheet(file, sheet = 'resources', n = 1)
  )
}

getSetupInfo <- function(state, params){
  # params$gain_mask <- is.na(state$gain)
  params$names <- list(
    resources = {resources = rownames(params$fixed_attributes)},
    states = {states = rownames(params$states)},
    habitats = {habitats = rownames(params$habitats)},
    types = apply(
      expand.grid(resources, states, habitats), 1, paste, collapse = '.'
    ),
    fixed_attributes = {fa = dimnames(params$fixed_attributes)},
    variable_attributes = {va = dimnames(state$resources)},
    attributes = c(fa[[2]], va[[2]]),
    nutrients = rownames(params$nutrients),
    toxins = rownames(params$toxins),
    use = dimnames(state$use),
    attention = rownames(state$attention),
    acceptance = rownames(state$accept),
    gain = dimnames(state$gain)
  )
  params$dims <- with(
    params$names, 
    list(
      N = length(resources),
      M = length(states),
      L = length(habitats),
      I = length(types),
      J = length(use[2]),
      K = length(attributes),
      Kvar = length(variable_attributes),
      gain = sapply(gain, length),
      variable_attributes = sapply(variable_attributes, length),
      fixed_attributes = sapply(fixed_attributes, length),
      use = sapply(use, length),
      attention = length(attention),
      acceptance = length(acceptance)
    )
  )
  params
}

setup <- function(parameter_file, initial_values_file, description = NULL){
  if(is.null(description)) 
    description <- format(
      paste('set up ', Sys.time(),
            '\n parameters = ', parameter_file, 
            '\n initial state = ', initial_values_file, 
            sep=''
      )
    )
  params <- importParameters(parameter_file)
  inits <- importStates(initial_values_file)
  params <- getSetupInfo(inits, params)
  params <- fixGains(inits, params)
  params <- fixTransitions(params)
  params <- c(
    grab('general', params, invert = TRUE),
    as.list(params$general[,'value'])
  )
  
  inits$rate <- formatRate(inits$rate, max = params$max_rate)
  inits$accept <- formatAcceptance(inits$accept)
  inits$use <- formatUse(inits$use)
  inits$attention <- formatAttention(inits$attention)
  inits$gain <- formatGain(inits$gain, params, setup = TRUE)
  inits <- unlist(lapply(inits, flattenMatrix))
  names(inits) <- sapply(names(inits), trimNames)
  
  list(
    description = description, 
    parameters = params, 
    state = inits
  )
}

# -----------------------------------------------------------------------------
# Parsers - these functions extract and format the information in the  
# parameters the current (usually linearized) state-values
parseInputs <- function(theta, not_theta, args){
  args$n <- parseDensities(not_theta, args)
  args$r <- parseRates(theta, args)
  args$a <- parseAcceptance(theta, args)
  args$p <- parseAttention(theta, args)
  args$W <- parseUsage(theta, args)
  args$G <- parseGain(not_theta, args)
  args$X <- parseAttributes(not_theta, args)
  args$k <- parseConspicuousness(not_theta, args)
  args$c <- parseBurial(not_theta, args)
  args$h <- parseHandlingTime(not_theta, args)
  args$C <- parseRisk(not_theta, args)
  args$e <- args$excavation_time
  args$m <- args$nutrients
  args$R <- args$max_rate
  args$b <- args$boldness
  args$v <- parseVigilance(theta, args)
  args
}

parseDensities <- function(values, args){
  pmax(grab('density', values), 0)
}

parseRates <- function(values, args){
  rates <- formatRate(
    grab('rate', values), 
    linearize = FALSE, 
    max = args$max_rate
  )
  rates <- matchParameters(rates, args$n)
  rates
}

parseVigilance <- function(values, args){
  matchParameters((grab('vigilance', values)), args$n) #exp
}

parseGain <- function(values, args){
  G <- G2 <- formatGain(grab('gain', values), args)
  to_modify <- with(args$names, grab(toxins, gain[[1]], invert=TRUE))
  G2[to_modify, 'cache'] <- G2[to_modify, 'cache'] / 10^args$recache_factor
  rbind(G, G2)
}

parseAttributes <- function(values, args){
  va <- array(
    grab('resources', values), 
    args$dims$variable_attributes, 
    dimnames = args$names$variable_attributes
  )
  X <- matchParameters(cbind(args$fixed_attributes, va), args$n)
  from_own_cache <- args$states[,'own']
  storage.mode(from_own_cache) <- 'logical'
  cbind(X * !from_own_cache, X * from_own_cache)
}

parseUsage <- function(values, args){
  U <- formatUse(
    grab('use', values), 
    linearize = FALSE, 
    args$dims$use, 
    args$names$use
  )
  matchParameters(U, args$n)
}

parseAttention <- function(values, args){
  A <- formatAttention(
    grab('attention', values), 
    linearize = FALSE, 
    args$names$attention
  )
  A['ignore'] <- 0
  A <- matchParameters(A, args$n, map = args$states[,'strategy'])
  normalizeAttention(A, args)
}

parseAcceptance <- function(values, args){
  A <- formatAcceptance(
    grab('accept', values), 
    linearize = FALSE
  )
  matchParameters(A, args$n)
}

parseConspicuousness <- function(values, args){
  k <- args$states[,'conspicuousness']
  storage.mode(k) <- 'numeric'
  matchParameters(k, args$n)
}

parseBurial <- function(values, args){
  b <- args$states[,'buried']
  storage.mode(b) <- 'logical'
  matchParameters(b, args$n)
}

parseHandlingTime <- function(values, args){
  matchParameters(args$handling_time, args$n)
}

parseRisk <- function(values, args){
  matchParameters(10^args$habitats, args$n)
}

# -----------------------------------------------------------------------------
# model definition
detection <- function(args){
  value <- with(as.list(args),{
    p <- round(p, 8)
    ifelse(
      r < R,
      ifelse(k > 0, p^(1/k) * (1 - (r/R)^k), 0),
      0
    )
  })
  value
}

totalExcavationTime <- function(args){
  value <- with(as.list(args),{
    n <- n/10000 # convert to per square meter basis
    retrieve <- grab('ownCache', n, ret.names=TRUE)
    n[retrieve] <- n[retrieve] * 10^args$retrieval
    e * sum(r * n * d * c)
    # e * sum(r * n * d * c)/sum(r * n * d)
  })
  if(is.na(value)) 0 else(value)
}

totalHandlingTime <- function(args){
  value <- with(as.list(args),{
    n <- n/10000 # convert to per square meter basis
    retrieve <- grab('ownCache', n, ret.names=TRUE)
    n[retrieve] <- n[retrieve] * 10^args$retrieval
    sum(r * n * a * d * W * h)
    # sum(a * d * W * h)
  })
  value
}

baseForagingRate <- function(args){
  value <- with(as.list(args),{
    # convert density to per square meter basis
    n <- n/10000 
    # modify the experienced deensity of own caches to account for 
    # differences increased retreval rates due to either:
    #   (1) concentration of caches in space
    #   (2) memory advantages in cache retieval
    # (or both)
    retrieve <- grab('ownCache', n, ret.names=TRUE)
    n[retrieve] <- n[retrieve] * 10^args$retrieval
    # calculate foraging rate
    f <- (r * a * d * n)/(1 + v + E + H)
    # convert units from items/second -> items/day
    f * 86400 * active   
  })
  value
}

resourceUse <- function(args){
  args$d <- detection(args)
  args$E <- totalExcavationTime(args)
  args$H <- totalHandlingTime(args)
  args$f <- baseForagingRate(args)
  args$U <- with(as.list(args), f * W) 
  args
}

attributeUse <- function(args){
  args$A <- with(as.list(args),{
    t(t(U) %*% X)
  })
  args
}

constrainTheta <- function(theta){
  a <- grab('attention', theta)
  sum(-mgcv::notExp(a^2 - 25))
}

netBenefit <- function(theta, not_theta, args){
  args <- parseInputs(theta, not_theta, args)   
  args <- resourceUse(args)
  args <- attributeUse(args)
  value <- with(as.list(args),{
    sum(A * G) - sum(p * C/(b * v))
  })
  value <- value + constrainTheta(theta) 
  if(is.finite(value)) return(value) else(browser())
}

