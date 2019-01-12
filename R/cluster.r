#' Data-clustering
#'
#' This package contains methods, which enables clustering in dataframes. Particularly useful for bio-mathematics, cognitive sciences, etc.
#'
#' \code{cluster(df, ...)}
#' @param df Tibble/Dataframe to be clustered. Method also possible with vectors.
#' @param groupby Defaults to \code{c()}. Specificies columns, by which data is to be preliminarily divided into groups, within which the clusters are to be built.
#' @param by Specifies the column(s) for geometric data, according to which the clusters are to be built.
#' @param isnear A function. This function operates pairs of entries in the columns with geometric data and returns \code{TRUE}/\code{FALSE} if entries are near. Defaults to a Euclidian metric.
#' @param dist Defaults to \code{Inf}. If the default euclidean metric is used for \code{isnear}, this is the maximum tolerated distance between geometric data.
#' @param strict Defaults to \code{TRUE}. If the default euclidean metric is used for \code{isnear}, this sets the proximity to be a strict \code{< dist} or else \code{<= dist}.
#' @param clustername Defaults to \code{'cluster'}. Running \code{df \%>\% clusterby(...)} returns a data frame, which extends \code{df} by 1 column with this name. This column tags the clusters by a unique index.
#' @param min_cluster_size Defaults to \code{0}. If a cluster is smaller than this, it will not be viewed as a cluster.
#' @param max_cluster_size Defaults to \code{Inf}. If a cluster is larger than this, it will be broken up into smaller pieces.
#' @param split Defaults to \code{FALSE}. If set to \code{TRUE}, then the output will be group the tibble data by cluster (equivalent to performing \code{\%>\% group_by(...)}).
#' @export cluster
#' @examples gene %>% cluster(groupby=c('gene','actve'), by='position', dist=400, strict=FALSE, clustername='tag');
#' @examples protein3d %>% cluster(groupby='celltype', by=c('x','y','z'), dist=2.5, max=5000);
#' @keywords cluster clustering gene




cluster <- function(data, ...) {
	INPUTVARS <- list(...);
	VARNAMES <- names(INPUTVARS);

	groupby <- INPUTVARS[['groupby']];
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
	data <- as_tibble(data);
	chunkname <- 'chunk';
	i <- 0;
	while(chunkname %in% c(cols, clustername)) {chunkname <- paste0('chunk', i); i <- i+1;}

	## Prägruppierung der Daten:
	if(!is.vector(groupby)) groupby <- c();
	leer <- list();
	n <- nrow(data);
	for(i in c(1:n)) leer[[i]] <- NA;
	data[, clustername] <- leer;
	# data <- add_column(data, !!(clustername) := leer);
	tib <- data %>% group_by_at(groupby) %>% nest(.key=!!(chunkname));
	m <- nrow(tib);

	for(i in c(1:m)) {
		chunk <- tib[i, chunkname][[1]][[1]];
		edges <- list();
		n <- nrow(chunk);
		hasedges <- FALSE;
		for(i in c(1:n)) {
			r <- 0;
			e <- c();
			if(i < n) {
				x <- chunk[i, by];
				for(j in c((i+1):n)) {
					y <- chunk[j, by];
					if(isnear(x, y)) {
						r <- r +  1;
						e[r] <- j;
					}
				}
			}
			if(length(e) == 0) {
				e <- NA;
			} else {
				hasedges <- TRUE;
			}
			edges[[i]] <- e;
		}
		if(hasedges) {
			clusters <- generateclasses(edges, min_cluster_size, max_cluster_size);
			ind <- which(!is.na(clusters));
			clusters <- clusters[ind];
			chunk <- chunk[ind, ];
			chunk[, clustername] <- clusters;
		} else {
			chunk <- chunk[c(), ];
		}
		tib[i, chunkname][[1]][[1]] <- chunk;
	}

	df <- tib %>% unnest();
	if(split) df <- df %>% group_by_at(c(groupby, clustername)); #%>% nest(.key='data');

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
		while(sz > 0) {
			children_ <- c();
			for(i in children) {
				e <- edges[[i]];
				if(length(e) == 1) if(is.na(e) || is.null(e)) next;
				children_ <- c(children_, e);
			}
			children <- children_;
			if(length(children) == 0) break;
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

	return(classes);
};