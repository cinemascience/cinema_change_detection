rss_sweep <- function(x, y, B, rss, an, dn) {
# Errors ----------------------------------------------------------------

# Checking for existence.
	if (exists("rss") == FALSE) {
		stop("Input the desired RSS.")
	}
	if (exists("an") == FALSE) {
		stop("Input the desired power value for alpha e.g. 10^(-an).")
	}
	if (exists("dn") == FALSE) {
		stop("Input the desired power value for deltasq e.g. 10^(-dn).")
	}
	
# Checking if values are integers.
	if (all(an == floor(an)) == FALSE) {
		stop("Power for alpha must be an integer.")
	}
	if (all(dn == floor(dn)) == FALSE) {
		stop("Power for deltasq must be an integer.")
	}
	
# Checking if values are numeric.
	if (all(sapply(x, is.numeric)) == FALSE) {
		stop("Input data must be numeric.")
	}
	if (is.numeric(rss) == FALSE) {
		stop("Input RSS must be numeric.")
	}
	if (is.numeric(an) == FALSE) {
		stop("Input an must be numeric.")
	}
	if (is.numeric(dn) == FALSE) {
		stop("Input an must be numeric.")
	}
	
# Checking if values are positive.
	if (any(an <= 0) == TRUE) {
		stop("Power for alpha must be positive.")
	}
	if (any(dn <= 0) == TRUE) {
		stop("Power for deltasq must be positive.")
	}

# Start main function ----------------------------------------------------------------	
	
# Making the range of values for alpha and deltasq based on inputs.
	alpha <- c(10^(-seq(1:an)))
	deltasq <- c(10^(-seq(1:dn)), 0)

	alpha.name <- paste(alpha)
	deltasq.name <- paste(deltasq)

	a.n <- length(alpha)
	d.n <- length(deltasq)

# Creating matrices that will hold the results of Total RSS and Number of Partitions.
	rss.tab <- matrix(, ncol = a.n, nrow = d.n)
	colnames(rss.tab) <- alpha.name
	rownames(rss.tab) <- deltasq.name
	part.tab <- matrix(, ncol = a.n, nrow = d.n)
	colnames(part.tab) <- alpha.name
	rownames(part.tab) <- deltasq.name
	
# Checking all partition values and their respective RSS values.
	for (i in 1:d.n) {
		for (j in 1:a.n) {
			z <- run_lm(x, y, B, alpha[j], deltasq[i])
			rss.tab[i, j] <- z$RSSTotal
			part.tab[i, j] <- z$Partition
		}
	}
	
# Finding the partition requested number closest to our results.
	rss.min <- abs(rss - rss.tab)
	ind <- which(rss.min == min(rss.min), arr.ind = TRUE)
	
	if (nrow(ind) > 1) {
	# If there are ties, it is broken which partition has the lowest RSS.
		part.min <- vector()
		ind.row <- dim(ind)[1]
		
		for (i in 1:ind.row) {
			part.min[i] <- part.tab[ind[i, 1], ind[i, 2]]
		}
		
		part.ind <- which(part.min == min(part.min))
		ind.min <- ind[part.ind, ]
	} else {
		ind.min <- ind
	}
	
	sweep.list <- list("alpha" = alpha[ind.min[2]],
		"deltasq" = deltasq[ind.min[1]]
	)
	
	return(sweep.list)
}