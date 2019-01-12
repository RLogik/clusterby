#' Data-clustering | summarise
#'
#' Provides customisable summary method tailored for clustering.
#'
#' \code{summarise(tib, ...)}
#' @param tib Tibble data which has been grouped.
#' @examples tib %>% summarise(concentration=mean, names=c('set',';'), attributes=c('list',';'));
#' @export summarise
#' @keywords cluster summarise




summarise <- function(tib, ...) {
	summaries <- list(...);
	summarycols <- names(summaries);
	cols <- names(tib);
	summarycols <- summarycols[which(summarycols %in% cols)];

	method <- list();
	for(col in summarycols) {
		s <- summaries[[col]];
		if(mode(s) == 'function') {
			f <- s;
		} else if(mode(s) == 'character') {
			if(s[1] == 'set') {
				sep = ';'; if(length(s) > 1) sep = s[2];
				f <- function(x) {return(paste(unique(x),collapse=sep));};
			} else if(s[1] == 'list') {
				sep = ';'; if(length(s) > 1) sep = s[2];
				f <- function(x) {return(paste(x,collapse=sep));};
			} else if(s[1] == 'length') {
				f <- length;
			} else if(s[1] == 'min') {
				f <- min;
			} else if(s[1] == 'max') {
				f <- max;
			} else if(s[1] == 'range') {
				f <- function(x) {return(paste(min(x),max(x),sep='-'));};
			} else if(s[1] == 'mean') {
				f <- mean;
			} else if(s[1] == 'var') {
				f <- var;
			} else if(s[1] == 'sd') {
				f <- sd;
			} else {
				f <- function(x) {return(NA);};
			}
		}
		method[[col]] <- f;
	}

	groupname <- 'cluster';
	i <- 0;
	while(groupname %in% cols) {
		groupname <- paste0('cluster', i);
		i <- i+1;
	}

	tib <- tib %>% nest(.key=groupname);
	n <- nrow(tib);
	for(i in c(1:n)) {
		tib_ <- tib[i, groupname][[1]][[1]];
		for(col in summarycols) {
			f <- method[[col]];
			w <- f(as.vector(tib_[, col]));
			tib_[, col] <- w;
		}
		tib[i, groupname][[1]][[1]] <- tib_;
	}
	tib <- tib %>% unnest();

	return(tib);
};