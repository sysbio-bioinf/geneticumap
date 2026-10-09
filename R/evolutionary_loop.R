evolution_loop <- function (
  fitness.fun,
  fitness.fun.name,
  minimize = NULL,
  n.objectives = 1L,
  n.dim = NULL, lower = NULL, upper = NULL,
  representation, mu, lambda,
  p.recomb = 0.7, p.mut = 0.3,
  survival.strategy = "plus", n.elite = 0L,
  log.stats = NULL,
  log.pop = FALSE,
  monitor = NULL,
  initial.solutions = NULL,
  parent.selector = NULL,
  survival.selector = NULL,
  mutator = NULL,
  max.iter = 100L,
  ...) {
  checkmate::assertChoice(representation, c("binary", "float", "permutation", "custom"))
  checkmate::assertChoice(survival.strategy, c("comma", "plus"))
  checkmate::assertNumber(p.recomb, lower = 0, upper = 1)
  checkmate::assertNumber(p.mut, lower = 0, upper = 1)
  checkmate::assertFlag(log.pop)
  mu <- checkmate::asInt(mu, lower = 1L)
  lambda.lower <- if (survival.strategy == "plus") 1L else mu
  lambda <- checkmate::asInt(lambda, lower = lambda.lower)

  #single-objective logging setup
  log_stats_list <- list(fitness = list())
  for (stat in log.stats) {
    if (stat == "mean") {
      log_stats_list$fitness[[paste(fitness.fun.name, stat)]] <- list(fun = function(fitness, ...) mean(as.numeric(fitness)))
    } else if (xor(stat == "worst", minimize)) {
      log_stats_list$fitness[[paste(fitness.fun.name, stat)]] <- list(fun = function(fitness, ...) min(as.numeric(fitness)))
    } else {
      log_stats_list$fitness[[paste(fitness.fun.name, stat)]] <- list(fun = function(fitness, ...) max(as.numeric(fitness)))
    }
  }

  control <- ecr::initECRControl(fitness.fun, minimize = minimize, n.objectives = n.objectives)
  control$type <- representation

  control <- ecr::registerECROperator(control, "mutate", mutator)
  control <- ecr::registerECROperator(control, "selectForMating", parent.selector)
  control <- ecr::registerECROperator(control, "selectForSurvival", survival.selector)

  log <- ecr::initLogger(control,
    log.stats = log_stats_list,
    log.pop = log.pop, init.size = 1000L
  )

  gen.fun <- ecr::genReal
  gen.pars <- list(n.dim = n.dim, lower = lower, upper = upper)

  population <- initial.solutions
  if (representation != "custom") {
    population <- do.call(ecr::initPopulation, c(list(mu = mu, gen.fun = gen.fun, initial.solutions = initial.solutions), gen.pars))
  }
  fitness <- ecr::evaluateFitness(control, population, ...)

  for (i in seq_along(population)) {
    attr(population[[i]], "fitness") <- fitness[, i]
  }

  ecr::updateLogger(log, population, fitness = fitness, n.evals = mu)

  n.iter <- 1L
  repeat {
    # Documentation says that generateOffspring passes down further arguments to the mutator, but it does not.
    # To access n.iter in the mutator, we therefore call it directly.
    offspring <- ecr::mutate(control, population[ecr::selectForMating(control, fitness, n.select = lambda)], p.mut = p.mut, n.iter = n.iter, lower = lower, upper = upper)
    fitness.offspring <- ecr::evaluateFitness(control, offspring, ...) #this runs in parallel batches per parallelization unit
    for (i in seq_along(offspring)) {
      attr(offspring[[i]], "fitness") <- fitness.offspring[, i]
    }

    sel <- if (survival.strategy == "plus") {
      ecr::replaceMuPlusLambda(control, population, offspring, fitness, fitness.offspring)
    } else {
      ecr::replaceMuCommaLambda(control, population, offspring, fitness, fitness.offspring, n.elite = n.elite)
    }

    #explicitly remove references to old generations
    rm(population, fitness, offspring, fitness.offspring)

    population <- sel$population
    fitness <- sel$fitness

    #remove duplicate pointers
    rm(sel)

    ecr::updateLogger(log, population, fitness, n.evals = lambda)
    if (is.function(monitor)) monitor()

    message(paste0("Finished iteration ", n.iter, " of ", max.iter))

    if (n.iter >= max.iter) {
      break
    }
    n.iter <- n.iter + 1L
  }
  return(make_ecr_result(control, log, population, fitness))
}

make_ecr_result <- function(control, log, population, fitness, ...) {
  # Single-objective result object
  return(BBmisc::makeS3Obj(
    task = control$task,
    best.x = log$env$best.x,
    best.y = log$env$best.y,
    log = log,
    last.population = population,
    last.fitness = as.numeric(fitness),
    classes = c("ecr_single_objective_result", "ecr_result")
  ))
}