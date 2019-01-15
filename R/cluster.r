#' Data-clustering
#'
#' This package contains methods, which enables clustering in dataframes. Particularly useful for bio-mathematics, cognitive sciences, etc.
#'
#' \code{cluster(df, ...)}
#' @param df Tibble/Dataframe to be clustered. Method also possible with vectors.
#' @param groupby Defaults to \code{c()}. Specificies columns, by which data is to be preliminarily divided into groups, within which the clusters are to be built.
#' @param by Specifies the column(s) for geometric data, according to which the clusters are to be built.
#' @param near A function. This function operates pairs of entries in the columns with geometric data and returns \code{TRUE}/\code{FALSE} if entries are near. Defaults to a Manhattan metric.
#' @param min.dist Defaults to \code{0}. If the default manhattan metric is used for \code{near}, this is the minimum tolerated distance between geometric data.
#' @param max.dist Defaults to \code{Inf}. If the default manhattan metric is used for \code{near}, this is the maximum tolerated distance between geometric data.
#' @param strict Defaults to \code{FALSE}. If the default manhattan metric is used for \code{near}, this sets the proximity to be a strict \code{< dist} or else \code{<= dist}.
#' @param clustername Defaults to \code{'cluster'}. Running \code{df \%>\% clusterby(...)} returns a data frame, which extends \code{df} by 1 column with this name. This column tags the clusters by a unique index.
#' @param min.size Defaults to \code{0}. If a cluster has fewer elements as this, it will not be viewed as a cluster.
#' @param split Defaults to \code{FALSE}. If set to \code{TRUE}, then the output will be group the tibble data by cluster (equivalent to performing \code{\%>\% group_by(...)}).
#' @param is.linear Defaults to \code{FALSE}. If set to \code{TRUE}, then the geometry is assumed to be linear and endowed with a simple difference-metric. This allows for faster computation.
#' @param is.disjoint Defaults to \code{FALSE}. If set to \code{TRUE} in combination with \code{is.linear=TRUE}, then the clusters must occupy disjoint intervals.
#' @export cluster
#' @examples gene %>% cluster(groupby=c('gene','actve'), by='position', dist=400, strict=TRUE, clustername='tag');
#' @examples protein3d %>% cluster(groupby='celltype', by=c('x','y','z'), dist=2.5,);
#' @keywords cluster clustering gene




cluster <- function(data, ...) {
	INPUTVARS <- list(...);
	VARNAMES <- names(INPUTVARS);

	groupby <- INPUTVARS[['groupby']];
	by <- INPUTVARS[['by']];
	clustername <- INPUTVARS[['clustername']];
	min_cluster_size <- INPUTVARS[['min.size']];
	split <- INPUTVARS[['split']];
	near <- INPUTVARS[['near']];
	d_min <- INPUTVARS[['min.dist']];
	d_max <- INPUTVARS[['max.dist']];
	strict <- INPUTVARS[['strict']];
	is_linear <- INPUTVARS[['is.linear']];
	is_disjoint <- INPUTVARS[['is.disjoint']];

	if(!is.vector(groupby)) groupby <- c();
	if(!is.logical(strict)) strict <- FALSE;
	if(!is.logical(split)) split <- FALSE;
	if(!is.logical(is_linear)) is_linear <- FALSE;
	if(!is.logical(is_disjoint)) is_disjoint <- FALSE;
	if(!is.numeric(d_min)) d_min <- 0;
	if(!is.numeric(d_max)) d_max <- Inf;
	if(!is.numeric(min_cluster_size)) min_cluster_size <- 0;
	if(!is.function(near)) {
		if(!is.character(near)) near <- 'Manhattan';
		if(near == 'Euclidean') {
			if(strict) {
				near <- function(x, y) {d <- sqrt(sum((x-y)^2)); return(d_min <= d && d < d_max);};
			} else {
				near <- function(x, y) {d <- sqrt(sum((x-y)^2)); return(d_min <= d && d <= d_max);};
			}
		} else {#if(near == 'Manhattan') {
			bool <- FALSE;
			if(strict) {
				near <- function(x, y) {d <- max(abs(x-y)); return(d_min <= d && d < d_max);};
			} else {
				near <- function(x, y) {d <- max(abs(x-y)); return(d_min <= d && d <= d_max);};
			}
		}
	}


	## Erstellung von Spaltennamen (Klusterspalte + Pufferspalte):
	data <- tibble::as_tibble(data);
	cols <- names(data);
	n <- nrow(data);
	if(!is.character(clustername)) {
		clustername <- 'cluster'; i <- 0; while(clustername %in% cols) {clustername <- paste0('cluster', i); i <- i+1;}
	}
	chunkname <- 'chunk'; i <- 0; while(chunkname %in% c(cols, clustername)) {chunkname <- paste0('chunk', i); i <- i+1;}
	edgesname <- 'edges'; i <- 0; while(edgesname %in% c(cols, clustername)) {edgesname <- paste0('edges', i); i <- i+1;}

	## Prägruppierung der Daten:
	leer <- list(); for(i in c(1:n)) leer[[i]] <- c(NA);
	tib <- add_column(data, !!(clustername):=leer, !!(edgesname):=leer);
	if(is_linear) tib <- tib[order(tib[[by]]), ];

	# Erzeuge Kanten für Klusterbausteine.
	tib <- data %>% group_by_at(groupby) %>% nest(.key=!!(chunkname));
	pos0 <- 0;
	m <- nrow(tib);
	for(i in c(1:m)) {
		chunk <- tib[i, chunkname][[1]][[1]];
		n <- nrow(chunk);
		if(n == 0) next;
		pts <- list(); for(j in c(1:n)) pts[[j]] <- chunk[j, by];
		if(is_linear) {
			edges <- lapply(c(1:n), function(j) {
				e <- c(NA);
				if(j < n) {
					pt <- pts[[j]];
					pt_ <- pts[[j+1]];
					if(near(pt, pt_)) e <- c(pos0 + j + 1);
				}
				return(e);
			});
		} else {
			edges <- lapply(c(1:n), function(j) {
				e <- c();
				if(j < n) {
					pt <- pts[[j]];
					pts_ <- pts[c((j+1):n)];
					bool <- lapply(pts_, function(y) {
						return(near(pt, y));
					});
					e <- pos0 + j + which(unlist(bool));
				}
				if(length(e) == 0) e <- c(NA);
				return(e);
			});
		}
		chunk[[edgesname]] <- edges;
		tib[i, chunkname][[1]][[1]] <- chunk;
		pos0 <- pos0 + n;
	}

	# Erzeuge Kluster aus Kanten.
	tib <- unnest(tib);
	clusters <- generateclasses(tib[[edgesname]], min_cluster_size);
	ind <- which(!is.na(clusters));
	clusters <- clusters[ind];
	tib <- tib[ind, ] %>% select(-c(edgesname));
	data <- add_column(tib, !!(clustername) := clusters);

	if(is_linear && is_disjoint) {
		data <- data[order(data[[by]], data[[clustername]]), ];
		n <- nrow(data);
		print(n);
		if(n > 0) {
			clusters <- data[[clustername]];
			cl_curr <- clusters[1];
			cl_replace <- 1;
			ind <- c();
			i0 <- 1;
			sapply(c(1:n), function(i) {
				cl <- clusters[i];
				if(!(cl == cl_curr)) {
					if(i-i0+1 >= min_cluster_size) ind <<- c(ind, c(i0:i));
					cl_curr <<- cl;
					cl_replace <<- cl_replace + 1;
				}
				data[i, clustername] <<- cl_replace;
				return(TRUE);
			});
			print(ind);
			data <- data[ind, ];
		}
	}

	if(split) data <- group_by_at(data, c(groupby, clustername)); #%>% nest(.key=!!(dataname);

	return(data);
};


generateclasses <- function(edges, min_sz) {
	n <- length(edges);
	if(n > 0) {
		key <- 0;
		ind <- c(1:n);
		classes <- rep(NA,n);

		while(length(ind) > 0) {
			i <- ind[1];
			ind <- ind[-1];
			nodes = c(i);
			children <- nodes;

			while(TRUE) {
				e <- edges[children];
				if(length(e) == 1) {
					grandchildren <- e[[1]];
				} else {
					grandchildren <- apply(cbind(edges[children]), 2, unlist)[,1];
				}
				if(length(grandchildren) == 0) break;
				grandchildren <- unique(grandchildren);
				children <- grandchildren[which(!is.na(grandchildren))];
				filt <- which(children %in% ind);
				if(length(filt) == 0) break;
				children <- children[filt];
				nodes <- c(nodes, children);
				ind <- ind[!(ind %in% children)];
			}

			if(length(nodes) < min_sz) next;

			classes[nodes] <- key;
			key <- key + 1;
		}
	}

	return(classes);
};