#' Data-clustering
#'
#' This package contains methods, which enables clustering in dataframes. Particularly useful for bio-mathematics, cognitive sciences, etc.
#'
#' \code{cluster(df, ...)}
#' @param df Tibble/Dataframe to be clustered. Method also possible with vectors.
#' @param groupby Defaults to \code{c()}. Specificies columns, by which data is to be preliminarily divided into groups, within which the clusters are to be built.
#' @param by Specifies the column(s) for geometric data, according to which the clusters are to be built.
#' @param near A function. This function operates pairs of entries in the columns with geometric data and returns \code{TRUE}/\code{FALSE} if entries are near. Defaults to a Euclidian metric.
#' @param max.dist Defaults to \code{Inf}. If the default euclidean metric is used for \code{near}, this is the maximum tolerated distance between geometric data.
#' @param strict Defaults to \code{TRUE}. If the default euclidean metric is used for \code{near}, this sets the proximity to be a strict \code{< dist} or else \code{<= dist}.
#' @param clustername Defaults to \code{'cluster'}. Running \code{df \%>\% clusterby(...)} returns a data frame, which extends \code{df} by 1 column with this name. This column tags the clusters by a unique index.
#' @param min.size Defaults to \code{0}. If a cluster has fewer elements as this, it will not be viewed as a cluster.
#' @param split Defaults to \code{FALSE}. If set to \code{TRUE}, then the output will be group the tibble data by cluster (equivalent to performing \code{\%>\% group_by(...)}).
#' @export cluster
#' @examples gene %>% cluster(groupby=c('gene','actve'), by='position', dist=400, strict=FALSE, clustername='tag');
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
	d_max <- INPUTVARS[['max.dist']];
	strict <- INPUTVARS[['strict']];

	if(!is.logical(strict)) strict <- FALSE;
	if(!is.logical(split)) split <- FALSE;
	if(!is.numeric(d_max)) d_max <- Inf;
	if(!is.numeric(min_cluster_size)) min_cluster_size <- 0;
	if(!is.function(near)) {
		if(!is.character(near)) near <- 'Euclidean';
		if(near == 'Manhattan') {
			bool <- FALSE;
			if(strict) {
				near <- function(x, y) {return(max(abs(x-y)) < d_max);};
			} else {
				near <- function(x, y) {return(max(abs(x-y)) <= d_max);};
			}
		} else { ## } else if(near == 'Euclidean') {
			if(strict) {
				near <- function(x, y) {return(sqrt(sum((x-y)^2)) < d_max);};
			} else {
				near <- function(x, y) {return(sqrt(sum((x-y)^2)) <= d_max);};
			}
		}
	}
	if(!is.character(clustername)) clustername <- 'cluster';


	## Erstellung von Spaltennamen (Klusterspalte + Pufferspalte):
	cols <- names(data);
	data <- tibble::as_tibble(data);
	chunkname <- 'chunk'; i <- 0; while(chunkname %in% c(cols, clustername)) {chunkname <- paste0('chunk', i); i <- i+1;}
	edgesname <- 'edges'; i <- 0; while(edgesname %in% c(cols, clustername)) {edgesname <- paste0('edges', i); i <- i+1;}

	## Prägruppierung der Daten:
	if(!is.vector(groupby)) groupby <- c();
	leer <- list();
	n <- nrow(data);
	for(i in c(1:n)) leer[[i]] <- c(NA);
	tib <- add_column(data, !!(clustername) := leer, !!(edgesname) := leer);

	# Erzeuge Kanten für Klusterbausteine.
	tib <- nest(group_by_at(data, groupby), .key=!!(chunkname));
	m <- nrow(tib);
	for(i in c(1:m)) {
		chunk <- tib[i, chunkname][[1]][[1]];
		n <- nrow(chunk);
		if(n == 0) next;
		pts <- list(); for(j in c(1:n)) pts[[j]] <- chunk[j, by];
		edges <- lapply(c(1:n), function(j) {
			e <- c();
			if(j < n) {
				pt <- pts[[j]];
				pts_ <- pts[c((j+1):n)];
				bool <- lapply(pts_, function(y) {
					return(near(pt, y));
				});
				e <- j + which(unlist(bool));
			}
			if(length(e) == 0) e <- c(NA);
			return(e);
		});
		chunk[[edgesname]] <- edges;
		tib[i, chunkname][[1]][[1]] <- chunk;
	}

	# Erzeuge Kluster aus Kanten.
	tib <- unnest(tib);
	clusters <- generateclasses(tib[[edgesname]], min_cluster_size);
	ind <- which(!is.na(clusters));
	clusters <- clusters[ind];
	data <- add_column(data[ind,], !!(clustername) := clusters);

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
				ind <- ind[!which(ind %in% children)];
			}

			if(length(nodes) >= min_sz) {
				classes[nodes] <- key;
				key <- key + 1;
			}
		}
	}

	return(classes);
};