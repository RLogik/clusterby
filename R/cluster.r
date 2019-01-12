#' Data-clustering
#'
#' This package contains methods, which enables clustering in dataframes. Particularly useful for bio-mathematics, cognitive sciences, etc.
#'
#' \code{cluster(df, ...)}
#' @param df Dataframe to be clustered. Method also possible with vectors.
#' @param group Defaults to \code{c()}. Specificies columns, by which data is to be preliminarily divided into groups, within which the clusters are to be built.
#' @param by Specifies the column(s) for geometric data, according to which the clusters are to be built.
#' @param isnear A function. This function operates pairs of entries in the columns with geometric data and returns \code{TRUE}/\code{FALSE} if entries are near. Defaults to a Euclidian metric.
#' @param dist Defaults to \code{Inf}. If the default euclidean metric is used for \code{isnear}, this is the maximum tolerated distance between geometric data.
#' @param strict Defaults to \code{TRUE}. If the default euclidean metric is used for \code{isnear}, this sets the proximity to be a strict \code{< dist} or else \code{<= dist}.
#' @param clustername Defaults to \code{'cluster'}. Running \code{df \%>\% clusterby(...)} returns a data frame, which extends \code{df} by 1 column with this name. This column tags the clusters by a unique index.
#' @param min_cluster_size Defaults to \code{0}. If a cluster is smaller than this, it will not be viewed as a cluster.
#' @param max_cluster_size Defaults to \code{Inf}. If a cluster is larger than this, it will be broken up into smaller pieces.
#' @param split Defaults to \code{FALSE}. If set to \code{TRUE}, then the output will be group the tibble data by cluster (equivalent to performing \code{\%>\% group_by(...)}).
#' @keywords cluster clustering gene
#' @export cluster
#' @examples gene %>% cluster(by='position', group=c('gene','actve'), dist=400, strict=FALSE, clustername='tag');
#' @examples protein3d %>% cluster(by=c('x','y','z'), group=c('gene','actve'), dist=2.5, max=5000);
#' \code{summarise(tib, ...)}
#' @param tib Tibble data which has been grouped.
#' @examples tib %>% summarise(concentration=mean, names=c('set',';'), attributes=c('list',';'));




cluster <- function(data, ...) {
	INPUTVARS <- list(...);
	VARNAMES <- names(INPUTVARS);

	tagcols <- INPUTVARS[['group']];
	by <- INPUTVARS[['by']];
	clustername <- INPUTVARS[['clustername']];
	min_cluster_size <- INPUTVARS[['min']];
	max_cluster_size <- INPUTVARS[['max']];
	splqit <- INPUTVARS[['split']];
	isnear <- INPUTVARS[['near']];
	d <- INPUTVARS[['dist']];
	strict <- INPUTVARS[['strict']];

	if(!is.logical(strict)) strict <- FALSE;
	if(!is.logical(split)) split <- FALSE;
	if(!is.numeric(d)) d <- Inf;
	if(!is.numeric(min_cluster_size)) min_cluster_size <- 0;
	if(!is.numeric(max_cluster_size)) max_cluster_size <- Inf;
	if(!is.function(isnear)) {
		if(strict) {
			isnear <- function(x, y) {
				dx <- sqrt(sum((x-y)^2));
				# dx <- max(abs(x-y));
				return(dx < d);
			};
		} else {
			isnear <- function(x, y) {
				dx <- sqrt(sum((x-y)^2));
				# dx <- max(abs(x-y));
				return(dx <= d);
			};
		}
	}
	if(!is.character(clustername)) clustername <- 'cluster';


	## Erstellung von Spaltennamen (Klusterspalte + Pufferspalte):
	cols <- names(data);
	data <- as.tibble(data);
	chunkname <- 'chunk';
	i <- 0;
	while(chunkname %in% c(cols, clustername)) {chunkname <- paste0('chunk', i); i <- i+1;}

	## Prägruppierung der Daten:
	if(!('group' %in% VARNAMES)) tagcols <- c();
	leer <- list();
	n <- nrow(data);
	for(i in c(1:n)) leer[[i]] <- c();
	data[[clustername]] <- leer;
	tib <- data %>% group_by_at(tagcols) %>% nest(.key=chunkname);
	n <- nrow(tib);

	for(i in c(1:n)) {
		df <- as.data.frame(tib[i, chunkname][[1]][[1]]);
		df_ <- df[, c(by, clustername)];
		n <- nrow(df_);
		edges <- list();
		for(i in c(1:n)) {
			r <- 0;
			e <- c();
			x <- df_[i, by];
			for(j in c((i+1):n)) {
				y <- df_[j, by];
				if(near(x, y)) {
					r <- r +  1;
					e[r] <- j;
				}
			}
			edges[[i]] <- e;
		}
		clusters <- generateclasses(edges, min_cluster_size, max_cluster_size);
		tib[i, chunkname][[1]][[1]][, clustername] <- clusters;
	}

	df <- tib %>% unnest();
	if(split) df <- df %>% group_by_at(c(tagcols, clustername)); #%>% nest(.key='data');

	return(NULL);
};


generateclasses <- function(edges, min_sz, max_sz) {
	if(!is.numeric(min_sz)) min_sz <- 0;
	if(!is.numeric(max_sz)) max_sz <- Inf;
	n <- length(edges);
	if(n == 0) return(c());

	getclass <- function(i, ind, sz) {
		nodes = c();
		children <- c(i);
		bool <- TRUE;
		while(sz > 0) {
			children_ <- c();
			for(i in children) children_ <- c(children_, edges[[i]]);
			children <- children_;
			if(length(children) > 0) break;
			filt <- which(children %in% ind);
			if(length(filt) >= sz) filt <- filt[c(1:sz)];
			children <- children[filt];
			nodes <- c(nodes, children);
			ind <- ind[-c(children)];
			sz <- sz - length(children);
		}
		return(list(indices=ind, nodes=nodes))
	};

	key <- 0;
	ind <- c(1:n);
	classes <- rep(NA,n);
	while(length(ind) > 0) {
		i <- ind[1];
		ind <- ind[-1];
		obj <- getclass(i, ind, max_sz);
		ind <- obj$indices;
		nodes <- obj$nodes;
		if(length(nodes) >= min_sz) {
			classes[nodes] <- key;
			key <- key + 1;
		}
	}
};


# summarise <- function(tib, ...) {
# 	summaries <- list(...);
# 	summarycols <- names(summaries);
# 	cols <- names(tib);
# 	summarycols <- summarycols[which(summarycols %in% cols)];

# 	method <- list();
# 	for(col in summarycols) {
# 		s <- summaries[[col]];
# 		if(mode(s) == 'function') {
# 			f <- s;
# 		} else if(mode(s) == 'character') {
# 			if(s[1] == 'set') {
# 				sep = ';'; if(length(s) > 1) sep = s[2];
# 				f <- function(x) {return(paste(unique(x),collapse=sep));};
# 			} else if(s[1] == 'list') {
# 				sep = ';'; if(length(s) > 1) sep = s[2];
# 				f <- function(x) {return(paste(x,collapse=sep));};
# 			} else if(s[1] == 'length') {
# 				f <- length;
# 			} else if(s[1] == 'min') {
# 				f <- min;
# 			} else if(s[1] == 'max') {
# 				f <- max;
# 			} else if(s[1] == 'range') {
# 				f <- function(x) {return(paste(min(x),max(x),sep='-'));};
# 			} else if(s[1] == 'mean') {
# 				f <- mean;
# 			} else if(s[1] == 'var') {
# 				f <- var;
# 			} else if(s[1] == 'sd') {
# 				f <- sd;
# 			} else {
# 				f <- function(x) {return(NA);};
# 			}
# 		}
# 		method[[col]] <- f;
# 	}

# 	groupname <- 'cluster';
# 	i <- 0;
# 	while(groupname %in% cols) {
# 		groupname <- paste0('cluster', i);
# 		i <- i+1;
# 	}

# 	tib <- tib %>% nest(.key=groupname);
# 	n <- nrow(tib);
# 	for(i in c(1:n)) {
# 		tib_ <- tib[i, groupname][[1]][[1]];
# 		for(col in summarycols) {
# 			f <- method[[col]];
# 			w <- f(as.vector(tib_[, col]));
# 			tib_[, col] <- w;
# 		}
# 		tib[i, groupname][[1]][[1]] <- tib_;
# 	}
# 	tib <- tib %>% unnest();

# 	return(tib);
# };