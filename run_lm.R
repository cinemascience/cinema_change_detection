run_lm <- function(x, y, B, alpha, deltasq) {
# Begin checking errors ----------------------------------------------------------------

# Checking for existence.
	if (exists("x") == FALSE) {
		stop("Input the time steps.")
	}
	if (exists("y") == FALSE) {
		stop("Input the data.")
	}
	if (exists("B") == FALSE) {
		stop("Input the desired size of buffer.")
	}
	if (exists("alpha") == FALSE) {
		stop("Input the desired level of significance (alpha).")
	}
	if (exists("deltasq") == FALSE) {
		stop("Input the desired amount of variance (deltasq).")
	}
	
# Checking buffer size.
	if (B < 3) {
		stop("Buffer size must be greater than 3.")
	}

# Checking if buffer is an integer.
	if (all(B == floor(B)) == FALSE) {
		stop("Buffer must be an integer.")
	}
	
# Checking if values are numeric.
	if (all(sapply(x, is.numeric)) == FALSE) {
		stop("Input time steps must be numeric.")
	}
	if (all(sapply(y, is.numeric)) == FALSE) {
		stop("Input data must be numeric.")
	}
	if (is.numeric(deltasq) == FALSE) {
		stop("Input deltasq must be numeric.")
	}
	if (is.numeric(alpha) == FALSE) {
		stop("Input alpha must be numeric.")
	}
	
# Checking if values are positive.
	if ((alpha <= 0) == TRUE) {
		stop("Input alpha must be positive.")
	}
	if ((deltasq < 0) == TRUE) {
		stop("Input deltasq must be positive.")
	}

# Checking if data is a vector.
	if (is.vector(x) == FALSE) {
		stop("Input data must be a vector.")
	}
# End checking errors ----------------------------------------------------------------

# Start smaller functions ----------------------------------------------------------------	
	# Using biglm to compute a linear fit using sufficient statistics
	F_fit <- function(x, y) {
		# Make a data frame, then generate and return the biglm object.
		data <- data.frame(x = x, y = y)
		return(biglm(y ~ x, data = data))
	}
	
	# Calculating the p-value for single-line versus two-line.
	F_pvalue <- function(RSS1, RSS2, p1, p2, deltasq, Tdot) {
		# If statement in the case that both RSS1 and RSS2 are 0.
		if(abs(RSS1 - RSS2) < .Machine$double.eps) { 
			pval <- 1
		} else {
			# Compute the modified F statistic
			num <- (RSS1 - RSS2) / (p2 - p1)				# Numerator
			denom <- (RSS2 + Tdot * deltasq) / (Tdot - p2)	# Denominator
			modF <- num / denom								# The modified F statistic
			
			# Return the critical value.
			pval <- pf(modF, p2 - p1, Tdot - p2, lower.tail = FALSE)
		}
		return(pval)
	}
	
	# Updating the biglm object with new data.
	F_updatefit <- function(obj, x, y) {
		# Make a data frame, then update and return the biglm object.
		data <- data.frame(x = x, y = y)
		return(update(obj, data))
	}	
# End smaller functions ----------------------------------------------------------------	
	
# Initial set-up -----------------------------------------------------------------------

	p1 <- 2							# Number of parameters in the single-line model
	p2 <- 4							# Number of parameters in the two-line model
	T <- length(y)					# Total number of time steps
	curr.ind <- 1:B					# curr indices for the entire simulation
	buff.ind <- (B + 1):(2 * B)		# buff indices for the entire simulation
	part.curr.ind <- curr.ind 		# curr indices for the partition
		
	# Use biglm objects to calculate curr, buff, and curr U buff.
	# The RSSs are in <object>$qr$ss. Tdot is in <object>$n.
	curr.lm <- F_fit(x[curr.ind], y[curr.ind])
	buff.lm <- F_fit(x[buff.ind], y[buff.ind])
	curr.buff.lm <- F_fit(x[c(curr.ind, buff.ind)], y[c(curr.ind, buff.ind)])
	
	# Data Frame to store results.
	results <- data.frame(matrix(vector(), 0, 5, dimnames = list(c(), c("RSS", "Time", "pvalue", "beta0", "beta1"))), stringsAsFactors = F)
	
# Start main algorithm -----------------------------------------------------------------
	while(buff.ind[B] <= T) {
		# Calculating RSS1 and RSS2
		RSS1 <- curr.buff.lm$qr$ss 				# single-line model
		RSS2 <- curr.lm$qr$ss + buff.lm$qr$ss 	# two-line model
	
		# Find the p-value for the modified F statistic.
		p.value <- F_pvalue(RSS1, RSS2, p1, p2, deltasq, curr.buff.lm$n)
		
		# If Else Statement on whether or not to partition the data.
		if(p.value < alpha) { # Rejected the null hypothesis (two-line is better than single-line).
			# Store the important information: RSS, Part, pvalue, beta0, beta1
			curr.lm.coef <- coef(F_fit(x[part.curr.ind], y[curr.ind]))			
			part.results <- c(RSS = curr.lm$qr$ss, Time = curr.ind[length(curr.ind)], pvalue = p.value, beta0 = curr.lm.coef[1], beta1 = curr.lm.coef[2])
			
			results[dim(results)[1] + 1,] <- part.results
			
			# Update the indices
			curr.ind <- buff.ind
			buff.ind <- (buff.ind[B] + 1):(buff.ind[B] + B)
			part.curr.ind <- 1:B	
			
			# Must stop if buff reaches end of simulation.
			if(buff.ind[B] > T) {
				break
			}
			
			# Compute biglm objects (curr, buff, and curr U buff) for the new regions
			curr.lm <- F_fit(x[curr.ind], y[curr.ind])
			buff.lm <- F_fit(x[buff.ind], y[buff.ind])
			curr.buff.lm <- F_fit(x[c(curr.ind, buff.ind)], y[c(curr.ind, buff.ind)])
		} else { # Fail to reject the null hypothesis (single-line is better than two-line). 
			
			# Update the indices of the oldest and newest time step in buff.
			old.ind <- buff.ind[1]
			new.ind <- buff.ind[B] + 1

			# Must stop if new set of time steps reaches end of simulation.
			if(new.ind > T) {
				break
			}

			# Update curr and buff indices.
			curr.ind <- c(curr.ind, old.ind)
			buff.ind <- c(buff.ind[2:B], new.ind)
			
			part.curr.ind <- 1:(length(part.curr.ind) + 1)
			
			# Update the biglm objects for curr and curr U buff.
			curr.lm <- F_updatefit(curr.lm, x[old.ind], y[old.ind])
			curr.buff.lm <- F_updatefit(curr.buff.lm, x[new.ind], y[new.ind])	
			
			# Recompute the biglm object for the updated buff.
			buff.lm <- F_fit(x[buff.ind], y[buff.ind])
		}
	}
	
	# Record results from last partition, because the algorithm is within B steps of the end of the simulation.
	part.curr.ind <- 1:(length(part.curr.ind) + B)									# curr indices within the last partition.
	curr.buff.lm.coef <- coef(F_fit(x[part.curr.ind], y[c(curr.ind, buff.ind)]))	# Calculating the regression coefficients for the last partition.
	
	part.results <- c(RSS = curr.buff.lm$qr$ss, Time = T, pvalue = NA, beta0 = curr.buff.lm.coef[1], beta1 = curr.buff.lm.coef[2])
	results[dim(results)[1] + 1,] <- part.results
	
# End main algorithm -----------------------------------------------------------------	

	# Creating list to return results. The following output are:
	# Cumulative RSS for the entire simulation, the location of the partition, p-value of the modified F Statistic, the intercept coefficient of the regression model, and the slope intercept of the regression model.
		part.list <- list("RSS" = results[, 1],
			"RSSTotal" = sum(results[, 1]),
			"Time" = results[, 2],
			"Partition" = length(results[, 2]),
			"pvalue" = results[, 3],
			"beta0" = results[, 4],
			"beta1" = results[, 5]
		)
	return(part.list)

}