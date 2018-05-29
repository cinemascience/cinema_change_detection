part_sweep <- function(x, y, B, partition, an, dn) {
# Errors ----------------------------------------------------------------
	if (exists("partition") == FALSE) {
		stop("Input the desired number of partitions.")
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
	if (all(partition == floor(partition)) == FALSE) {
		stop("The number of partitions must be an integer.")
	}
	
# Checking if values are numeric.
	if (is.numeric(partition) == FALSE) {
		stop("Input number of partitions must be numeric.")
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
	if (any(partition <= 0) == TRUE) {
		stop("The number of partitions must be positive.")
	}
	
# Checking if data is a vector.
	if (is.vector(x) == FALSE) {
		stop("Input data must be a vector.")
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
	part.min <- abs(partition - part.tab)
	ind <- which(part.min == min(part.min), arr.ind = TRUE)
	
	if (nrow(ind) > 1) {
	# If there are ties, it is broken by the partition that has the lowest RSS.
		rss.min <- vector()
		ind.row <- dim(ind)[1]
		
		for (i in 1:ind.row) {
			rss.min[i] <- rss.tab[ind[i, 1], ind[i, 2]]
		}
		
		rss.ind <- which(rss.min == min(rss.min))
		ind.min <- ind[rss.ind, ]
	} else {
		ind.min <- ind
	}
	
	sweep.list <- list("alpha" = alpha[ind.min[2]],
		"deltasq" = deltasq[ind.min[1]]
	)
	
	return(sweep.list)
}