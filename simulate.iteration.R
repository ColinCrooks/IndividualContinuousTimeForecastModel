#' Function to sample one simulation of the cohort over the forecast horizon.
#'
#' All patient tables require the following variables:
#' \describe{
#' \item {patient_id}{integer representing the unique patients}
#' \item {covariates}{all variables needed for prediction/MI plus forescalation}
#' \item {day_diagnosis}{day from diagnosis - this is updated for each transition time so can be used as a predictor}
#' \item {day_inpatient}{day from admission - this is updated for each transition time so can be used as a predictor}
#' \item {state}{integer label of state of patient at that time point}
#' \item {postICU}{binary flag set once patient enters ICU and updated in transitions - can be used as a predictor}
#' \item {day_in_state}{numeric marking time already spent in starting state at beginning of transition for initial truncated sampling}
#'
#' @param sample_i integer for labelling this iteration in the returned matrix
#' @param patients data table containing all patients at date time at start of simulation
#' @param new.patient.sample data table of patients to sample for new admissions with same columns as patients table
#' @param ICU.patient.sample data table of patients to sample for new ICU out of area admissions with same columns as patients table
#' @param tmat matrix of eligible transitions with NA for those not possible with format as defined in mstate or hesim
#' @param states data table of states that are not absorbing with state_id - numeric number of state, and state_name - name of each state
#' @param future.admissions numeric vector of daily predicted admission numbers
#' @param cov.tv  string vector of continuous time varying covariates to be updated at each new transition
#' @param coef.grad matrix of slope coefficients (column 1) and slope standard error (column 2) for longitudinal modelling of covariates listed in cov.tv
#' @param mi.omit.predictors string vector of covariate names to not use in predicting imputations (Important to avoid collinearity)
#' @param mi.omit.imputed string vector of covariate names to not impute missing values for (Hopefully automatically omitted as shouldn't be missing).
#' @param diffICU integer of additional patients to be added from day 1 from ICU sample at start of cohort
#' @param ICU.transfers numeric vector of daily ICU transfers as out of area admissions to be included
#' @param flexsurv_list list of flexsurvreg or flexsurvspline objects for each transition from the flexsurv package.
#'
#'
#' @return  all.trans - data.table containing
#' \describe{
#' \item{sample}{Given label for this simulation}
#' \item{patient_id}{patient_id as given in patients}
#' \item{to}{state transitioned into at that time}
#' }
#'
#' @import data.table
#' @import mice
#'
simulate.iteration <- function(sample_i,
															 patients,
															 new.patient.sample,
															 ICU.patient.sample,
															 tmat,
															 states,
															 future.admissions,
															 cov.tv,
															 coef.grad,
															 mi.omit.predictors,
															 mi.omit.imputed,
															 diffICU,
															 ICU.transfers,
															 flexsurv_list) {
	library(mice)
	library(data.table)

	forecast.length <- length(future.admissions)
	ICU.state <- which(colnames(tmat) %in% "ICU")
	died.state <- which(colnames(tmat) %in% "Died")
	ED.state <- which(colnames(tmat) %in% "ED")
	Ward.state <- which(colnames(tmat) %in% "Ward")
	start_states <- patients$state

	all.trans <-
		data.table::data.table(
			sample = numeric(),
			patient_id = numeric(),
			time_stop = numeric(),
			to = numeric()
		)

	# Helper function adapted from flexsurv package to add covariates to survival model parameters.
	add.covs <- function(surv.model, X) {
		nres <- nrow(X)
		pars <-
			matrix(
				surv.model$res.t[surv.model$dlist$pars, "est"],
				nrow = nres,
				ncol = length(surv.model$dlist$pars) ,
				byrow = T
			)
		xb <-
			if (surv.model$ncovs == 0)
				0
		else
			surv.model$res.t[surv.model$covpars, "est"] %*% t(X)
		for (j in seq(along = surv.model$dlist$pars)) {
			covinds <- surv.model$mx[[surv.model$dlist$pars[j]]]
			if (length(covinds) > 0) {
				pars[, j] <- pars[, j] + xb
			}
			pars[, j] <- surv.model$dlist$inv.transforms[[j]](pars[, j])
		}
		colnames(pars) <- surv.model$dlist$pars
		return(pars)
	}

	# Additional helper function adapted from flexsurv package for constructing model matrix from data table
	compress.model.matrices <- function(mml) {
		cbind.drop.intercept <-
			function(...)
				do.call("cbind", lapply(list(...), function(x)
					x[, -1, drop = FALSE]))
		X <- do.call("cbind.drop.intercept", mml)
		loc.cnames <- colnames(mml[[1]])[-1]
		anc.cnames <-
			unlist(mapply(
				function(x, y)
					sprintf("%s(%s)", x, y),
				names(mml[-1]),
				lapply(mml[-1], function(x)
					colnames(x)[-1])
			))
		cnames <- c(loc.cnames, anc.cnames)
		colnames(X) <- cnames
		X
	}

	# Helper function adapted from flexsurv package to construct model matrix from data.table and expand dummy variables for factors
	form.model.matrix <-
		function(object,
						 newdata,
						 na.action = na.pass,
						 forms = NULL) {
			mfo <- model.frame(object)

			covnames <- attr(mfo, "covnames")
			tt <- attr(mfo, "terms")
			Terms <- delete.response(tt)
			mf <-
				model.frame(Terms,
										newdata,
										xlev = .getXlevels(tt, mfo),
										na.action = na.action)
			if (!is.null(cl <- attr(Terms, "dataClasses")))
				.checkMFClasses(cl, mf)
			if (is.null(forms))
				forms <- object$all.formulae
			mml <- vector(mode = "list", length = length(forms))
			names(mml) <- names(forms)
			forms[[1]] <- delete.response(terms(forms[[1]]))
			for (i in seq_along(forms)) {
				mml[[i]] <- model.matrix(forms[[i]], mf)
			}
			X <- compress.model.matrices(mml)

			attr(X, "newdata") <-
				mf # newdata with any extra variables stripped.  Used to name components of summary list
			X
		}



	# Sample additional patients for ICU if needed
	newpatients <-
		patients[start_states == ICU.state][sample(1:nrow(patients[start_states == ICU.state]) , size = diffICU)]
	newpatients[, patient_id := max(patients$patient_id) + .I]
	patients.added <-
		rbind(patients[start_states %in% c(ED.state, Ward.state, ICU.state)], newpatients)
	start_states.added <-
		c(start_states[start_states %in% c(ED.state, Ward.state, ICU.state)], rep(ICU.state, diffICU))[order(patients.added$patient_id)]
	patients.added <- patients.added[order(patient_id)]
	patients.added[, postICU := postICU == 1]
	patients.added[state == ICU.state, postICU := T] # needed to be logical
	npatients <- nrow(patients.added)
	new.npatients <- nrow(new.patient.sample)

	patients.comb <- data.table::rbindlist(list(patients.added,new.patient.sample,ICU.patient.sample))


	# Multiple imputation of missing variables - one sample for each iteration.
	init <- mice::mice(patients.comb,  maxit = 0)
	meth <- init$method
	predM <-
		mice::quickpred(patients.comb)


	patients.imputed <-
		data.table::data.table(mice::complete(
			mice::mice(
				patients.comb,
				method = meth,
				predictorMatrix = predM,
				m = 1,
				maxit = 50,
				print = FALSE
			)
		))
	rm(patients.comb)
	new.patient.sample.imput <- patients.imputed[(npatients + 1):((npatients  + new.npatients)),]
	new.patient.sample.imput[, postICU := postICU == 1]
	new.patient.sample.imput[state == ICU.state, postICU := T] # needed to be logical


	ICU.patient.sample.imput <- patients.imputed[((npatients + 1 + new.npatients)):(nrow(patients.imputed)),]


	patients.imputed <- patients.imputed[1:npatients,]

	# Set up vectors to accumulate information during iterations
	all.patients <-
		current.patients <-
		data.table::setDT(patients.imputed)[order(patient_id),]
	current.day <- rep(0, nrow(all.patients))
	current.state <- current.to <- start_states.added #

	all.trans <- data.table::rbindlist(list(
		all.trans,
		data.table::data.table(
			sample = rep(sample_i, nrow(patients.imputed)),
			patient_id = patients.imputed$patient_id,
			time_stop = rep(0, nrow(patients.imputed)),
			to = start_states.added
		)
	))

	for (day in 0:(forecast.length - 1)) {
		# To sample from the truncated parametric survival distribution, over sample a number of projected times and only select those samples with a first transition from the future time period
		if (day == 0) {
			next.trans <-
				data.table::data.table(patient_id = numeric(),
															 to = numeric(),
															 time_stop = numeric())
			for (from.state in unique(current.to)) {
				# for each of the states in the current cohort
				if (from.state %in% states$state_id) {
					# Do not simulate further transitions for death state

					## simulate next time and states for people whose current state is from.state
					# Plausible transition states from from.state
					transi <-
						tmat[from.state, is.finite(tmat[from.state, ]), drop = F]

					# Number of patients who could make this transition
					ni <- sum(current.to == from.state)

					# Matrix to store the sampled time for each transition
					t.trans1 <- matrix(-Inf, nrow = ni, ncol = length(transi))

					## simulate times to all potential destination states
					# For each potential transition j index
					for (j in 1:length(transi)) {
						# label of potential to state indexed j
						trans.to <- transi[j]

						# Store model matrix
						X <-
							matrix(rep(0, length(current.to == from.state)), nrow = 1)

						if (length(flexsurv_list[[trans.to]]$covpars) > 0)	{
							# Patients eligible for transition selected (current.to == from.state)
							# covariates in transition model indexed as were included in the model (flexsurv_list[[trans.to]]$covpars)
							X <- form.model.matrix(
								object = flexsurv_list[[trans.to]],
								newdata = current.patients[current.to == from.state],
								na.action = na.omit
							)
						}


						# Baseline parameters for survival function from model combined with the covariates and coefficients (flexsurv_list[[trans.to]]$covpars)
						basepars <- as.list(as.data.frame(add.covs(
							surv.model =
								flexsurv_list[[trans.to]],
							X = X
						)))
						# Identify truncated uniform distribution to sample from the cumulative survival distribution
						fncall <-
							c(list(q = current.patients[current.to == from.state, day_in_state]),
								basepars,
								flexsurv_list[[trans.to]]$aux)
						if (is.null(flexsurv_list[[trans.to]]$dfns$p))
							stop("No random sampling function found for this model")

						p.sample <- do.call(flexsurv_list[[trans.to]]$dfns$p, fncall)


						# function call from flexsurv to sample transition times from the parametric cumulative survival distribution for transition within truncated distribution
						fncall <-
							c(list(p = runif(
								n = length(p.sample),
								min = p.sample,
								max = 1
							)), basepars, flexsurv_list[[trans.to]]$aux)
						if (is.null(flexsurv_list[[trans.to]]$dfns$q))
							stop("No random sampling function found for this model")

						t.trans1[, j] <-
							do.call(flexsurv_list[[trans.to]]$dfns$q, fncall) # Sample a time for each patient

						# If not for escalation then ignore transitions to ICU
						if (match(transi[j], tmat[from.state, ])  == ICU.state)
							t.trans1[current.patients[current.to == from.state, forescalation == 0], j] <-
							Inf

					}



					mc <-
						max.col(-t.trans1)																			# Identify the earliest transition time for each patient
					next.state <-
						match(transi[mc], tmat[from.state, ]) 						# Identify the transition state for the selected transition
					next.time <-
						t.trans1[cbind(seq_along(next.state), mc)]				# Identify the transition time for the selected transition
					next.trans <- rbind(
						next.trans[, .(patient_id, to, time_stop)],
						data.table::data.table(
							patient_id = current.patients[current.to == from.state, patient_id],
							to = next.state,
							time_stop = next.time
						)
					)   # Store the transitions for patients in from.state
				}
			}
			next.trans <- current.patients[next.trans,
																		 on = 'patient_id']
			next.trans[, `:=`(time_stop = time_stop - day_in_state)] # Set time to be from beginning of simulation
			next.trans[time_stop < 0, time_stop := 0.1]  # Shouldn't be possible but catch all to prevent crashing out in big simulation
		}

		# Sample from full untruncated distribution for subsequent transitions
		if (day > 0) {
			# Container of patients who have not yet transitioned out of current day
			patients.withinday <- current.patients
			patients.withinday[, `:=`(last.day = 0, last.state = current.to)]
			step.state <-
				current.to # Remaining states to update within the day
			next.trans <-
				data.table::data.table(patient_id = numeric(),
															 to = numeric(),
															 time_stop = numeric()) # cumulative container of transitions
			while (length(step.state) > 0) {
				for (from.state in unique(step.state)) {
					# for each of the states in the current cohort
					if (from.state %in% states$state_id) {
						# Do not simulate further transitions for death state

						## simulate next time and states for people whose current state is from.state
						# Plausible transition states from from.state
						transi <-
							tmat[from.state, is.finite(tmat[from.state, ]), drop = F]

						# Number of patients who could make this transition
						ni <- sum(step.state == from.state)

						# Matrix to store the sampled time for each transition
						t.trans1 <- matrix(-Inf, nrow = ni, ncol = length(transi))

						## simulate times to all potential destination states for patients in from.state
						patients.todo <- step.state == from.state

						# For each potential transition j index
						for (j in 1:length(transi)) {
							# label of potential to state indexed j
							trans.to <- transi[j]

							# Store model matrix
							X <-
								matrix(rep(0, ni), nrow = 1)

							if (length(flexsurv_list[[trans.to]]$covpars) > 0)	{
								# Patients eligible for transition selected (current.to == from.state)
								# covariates in transition model indexed as were included in the model (flexsurv_list[[trans.to]]$covpars)
								X <- form.model.matrix(
									object = flexsurv_list[[trans.to]],
									newdata = patients.withinday[step.state == from.state],
									na.action = na.omit
								)
							}


							# Baseline parameters for survival function from model combined with the covariates and coefficients (flexsurv_list[[trans.to]]$covpars)
							basepars <- as.list(as.data.frame(add.covs(
								surv.model =
									flexsurv_list[[trans.to]],
								X = X
							)))
							# function call from flexsurv to sample transition times from the parametric survival distribution for transition
							fncall <-
								c(list(n = ni), basepars, flexsurv_list[[trans.to]]$aux)

							if (is.null(flexsurv_list[[trans.to]]$dfns$r))
								stop("No random sampling function found for this model") # Catch for me if I use ineligible distribution

							t.trans1[, j] <-
								do.call(flexsurv_list[[trans.to]]$dfns$r, fncall) # Sample a time for each patient

							# Ignore ICU transition if not for escalation
							if (match(transi[j], tmat[from.state, ])  == ICU.state)
								t.trans1[patients.withinday[step.state == from.state, forescalation == 0], j] <-
								Inf


						}

						mc <-
							max.col(-t.trans1)																			# Identify the earliest transition time for each patient
						next.state <-
							match(transi[mc], tmat[from.state, ]) 						# Identify the transition state for the selected transition
						next.time <-
							t.trans1[cbind(seq_along(next.state), mc)]	+ patients.withinday[patients.todo,	last.day]		# Identify the transition time for the selected transition
						patients.withinday[patients.todo, `:=`(last.day = next.time, last.state = next.state)]
						next.trans <- rbind(
							next.trans,
							data.table::data.table(
								patient_id = patients.withinday[patients.todo, patient_id],
								to = next.state,
								time_stop = next.time
							)
						)   # Store the transitions for patients in from.state

					}
				}
				# Make sure everyone has moved on at least one day or died by repeating sampling for those who haven't
				step.state <-
					patients.withinday[last.day <= 1 &
														 	last.state != died.state, last.state]
				patients.withinday <-
					patients.withinday[last.day <= 1 & last.state != died.state, ]
			}
			next.trans <- current.patients[next.trans,
																		 on = 'patient_id']
			next.trans[time_stop < 0, time_stop := 0.1]  # Shouldn't be possible but catch all to prevent crashing out in big simulation			}
		}


		next.trans[, `:=`(sample = sample_i,
											incr = time_stop)]

		# Update covariates for post ICU
		next.trans[(data.table::shift(to, n = 1L, type = 'lag') == ICU.state) &
							 	patient_id == data.table::shift(patient_id, n = 1L, type = 'lag') ,
							 postICU := 1]
		next.trans[to == ICU.state, postICU := 1]

		# returned stop time updated to be from start of simulation
		next.trans$time_stop <-
			next.trans$time_stop +	current.day[match(next.trans$patient_id, all.patients$patient_id)]

		#select first transition that ends in next day or later
		current.trans <- next.trans[(time_stop) > day , ][order(sample, patient_id, time_stop),
																											head(.SD, 1),
																											by = .(sample, patient_id)]

		#Store update including all transitions up to next day or later, so that state rolled forward reflects the last change even if within the day
		# Otherwise lose ED transitions and patients remain on initial starting state in ED
		all.trans <- data.table::rbindlist(list(all.trans,
																						next.trans[(time_stop) <= day ,
																											 .(sample, patient_id, time_stop, to)],
																						current.trans[, .(sample, patient_id, time_stop, to)]))

		# Update times for patients transitioned on this day
		current.day[match(current.trans$patient_id, all.patients$patient_id)] <-
			current.trans$time_stop

		# Stop transitions when reach accumulating state
		current.state[match(current.trans$patient_id, all.patients$patient_id)] <-
			current.trans$to
		current.day[current.state == died.state] <- Inf

		#update covariates by increment from last update (time_stop)

		all.patients[current.state == ICU.state , postICU := T]

		all.patients[match(current.trans$patient_id, patient_id),
								 `:=`(
								 	day_diagnosis = day_diagnosis + current.trans$incr,
								 	day_inpatient = day_inpatient +  current.trans$incr
								 )]
		all.patients[match(current.trans$patient_id, patient_id),
								 (cov.tv) := lapply(cov.tv,  function(x)
								 	all.patients[match(current.trans$patient_id, patient_id), get(x)] +
								 		(
								 			log(current.trans$incr) *
								 				rnorm(nrow(current.trans),
								 							mean = coef.grad[1, x],
								 							sd = coef.grad[2, x])
								 		))]

		# Identify patients for sampling on next iteration whose last transition ends within the next day
		current.to <-
			current.state[(current.day) > day & (current.day) <= day + 1]
		current.patients <-
			all.patients[(current.day) > day &
									 	(current.day) <= day + 1, ][order(patient_id), ]

		# Sample the number of next day admissions from the predictions allowing expected error around count prediction using Poisson error
		sampleN <-
			round(rpois(n = 1, lambda = future.admissions[day + 1]))

		# Sample new patients from recent admissions in last fortnight
		newpatients <-
			new.patient.sample.imput[sample(
				x = 1:nrow(new.patient.sample.imput),
				size = sampleN,
				replace = T
			),]

		newpatients[, `:=`(sample = sample_i,
											 time_stop = day,
											 to = state)]
		if (sampleN > 0)
			newpatients[, `:=`(patient_id = .I + max(all.patients$patient_id))]

		# Sample new ICU transfers from other providers if specified
		new.ICUtransfers <-
			ICU.patient.sample.imput[sample(
				x = 1:nrow(ICU.patient.sample.imput),
				size = ICU.transfers[day + 1],
				replace = T
			),]


		new.ICUtransfers[, `:=`(
			sample = sample_i,
			time_stop = day,
			to = state,
			postICU = 1
		)]
		if (ICU.transfers[day + 1] > 0)
			new.ICUtransfers[, `:=`(patient_id = .I + max(newpatients$patient_id))]

		# Add beginning states of additional patients into current data table of transitions
		all.trans <- data.table::rbindlist(list(all.trans,
																						newpatients[, .(sample, patient_id, time_stop, to)],
																						new.ICUtransfers[, .(sample, patient_id, time_stop, to)]))

		# Combine next day patients, new admissions and new ICU patients for next day sampling
		current.patients <- data.table::rbindlist(list(current.patients,
																									 newpatients[, !c('sample', 'time_stop', 'to')],
																									 new.ICUtransfers[, !c('sample', 'time_stop', 'to')]))

		if (nrow(current.patients) < 1)
			# Nothing to do!
		{
			next
		}

		# Combine all current patients, and newly sampled new admissions and new ICU patients to overall cohort
		all.patients <- data.table::rbindlist(list(all.patients,
																							 newpatients[, !c('sample', 'time_stop', 'to')],
																							 new.ICUtransfers[, !c('sample', 'time_stop', 'to')]))


		current.day <-
			as.vector(c(current.day, rep(day, sampleN + ICU.transfers[day + 1])))
		current.state <-
			as.vector(c(current.state, newpatients$state, new.ICUtransfers$state))
		current.to <-
			as.vector(c(current.to, newpatients$state, new.ICUtransfers$state))

	}

	# Return data.table output of all transitions for binding in rbindlist
	all.trans

}
