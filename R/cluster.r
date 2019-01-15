#' Data-clustering
#'
#' This package contains methods, which enables clustering in dataframes. Particularly useful for bio-mathematics, cognitive sciences, etc.
#'
#' \code{cluster(df, ...)}
#' @param df Tibble/Dataframe to be clustered. Method also possible with vectors.
#' @param by string vector. Specifies the column(s) for geometric data, according to which the clusters are to be built.
#' @param filter.by string vector. Defaults to \code{c()}. Specificies columns, by which data is to be preliminarily divided into groups, within which the clusters are to be built.
#' @param near function. This function operates pairs of entries in the columns with geometric data and returns \code{TRUE}/\code{FALSE} if entries are near. Defaults to a Manhattan metric.
#' @param min.dist a real number. Defaults to \code{0}. If the default manhattan metric is used for \code{near}, this is the minimum tolerated distance between geometric data.
#' @param max.dist a real number. Defaults to \code{Inf}. If the default manhattan metric is used for \code{near}, this is the maximum tolerated distance between geometric data.
#' @param strict boolean. Defaults to \code{FALSE}. If the default manhattan metric is used for \code{near}, this sets the proximity to be a strict \code{< dist} or else \code{<= dist}.
#' @param clustername string. Defaults to \code{'cluster'}. Running \code{df \%>\% clusterby(...)} returns a data frame, which extends \code{df} by 1 column with this name. This column tags the clusters by a unique index.
#' @param min.size a natural number. Defaults to \code{0}. If a cluster has fewer elements as this, it will not be viewed as a cluster.
#' @param split boolean. Defaults to \code{FALSE}. If set to \code{TRUE}, then the output will be group the tibble data by cluster (equivalent to performing \code{\%>\% group_by(...)}).
#' @param is.linear boolean. Defaults to \code{FALSE}. If set to \code{TRUE}, then the geometry is assumed to be linear and endowed with a simple difference-metric. This allows for faster computation.
#' @param is.disjoint boolean. Defaults to \code{FALSE}. If set to \code{TRUE} in combination with \code{is.linear=TRUE}, then the clusters must occupy disjoint intervals.
#' @param summary boolean. Defaults to \code{FALSE}. If set to \code{TRUE} in combination with \code{is.linear=TRUE} and **assuming** the user has presorted the data by the \code{by}-column, then a summary of the clusters as intervalls is provided. This makes most sense, if \code{is.disjoint=TRUE}. This produces the columns \code{filter.by, by, pstart, pend, nstart, nend, n} where \code{pstart}, \code{pend} describes the interval, \code{nstart}, \code{nend} provides the original indices in the input data, and \code{n} is the cluster size (number of points).
#' @export cluster
#' @examples gene %>% cluster(by='position', filter.by=c('gene','actve'), min.size=4, max.dist=400, strict=FALSE, is.linear=TRUE, is.disjoint=TRUE, summary=TRUE);
#' @examples protein3d %>% cluster(by=c('x','y','z'), filter.by='celltype', max.dist=5.8e-7, clustername='segment');
#' @examples soil_data %>% cluster(by=c('x','y'), filter.by=c('density','substance'), max.dist=10e-3, clustername='clump');
#' @keywords cluster clustering gene




cluster <- function(data, ...) {
	INPUTVARS <- list(...);
	VARNAMES <- names(INPUTVARS);

	group_by <- INPUTVARS[['filter.by']];
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
	is_disjoint <- INPUTVARS[['is.disjoint']];
	do_summary <- INPUTVARS[['summary']];

	if(!is.vector(group_by)) group_by <- c();
	if(!is.logical(strict)) strict <- FALSE;
	if(!is.logical(split)) split <- FALSE;
	if(!is.logical(is_linear)) is_linear <- FALSE;
	if(!is.logical(is_disjoint)) is_disjoint <- FALSE;
	if(!is.logical(do_summary)) do_summary <- FALSE;
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
	tib <- tibble::as_tibble(data);
	cols <- names(tib);
	n <- nrow(tib);
	if(!is.character(clustername)) clustername <- uniquecolumnname('cluster', cols);
	chunkname <- uniquecolumnname('chunk', c(cols, clustername));


	if(is_linear && is_disjoint) {
		## Füge Indizes für Gruppennamen und zur Erhaltung der Reihenfolge hinzu:
		indexname <- uniquecolumnname('index', c(cols, clustername));
		groupname <- uniquecolumnname('group', c(cols, clustername));
		tib <- tib %>% add_column(!!(indexname):=c(1:n));

		## Indizes für Präsortierung hinzufügen.
		tib <- tib %>% group_by_at(group_by) %>% nest(.key=!!(chunkname));
		m <- nrow(tib);
		tib <- tib %>% add_column(!!(groupname):=c(1:m));
		tib <- tib %>% unnest();

		## Reihenfolge wiederherstellen.
		tib <- tib[order(tib[[indexname]]), ];
		tib <- tib %>% select(-c(indexname));

		if(do_summary) {
			summcols <- c('pstart','pend','nstart','nend','n');
			strich <- '';
			pref <- function(x) {return(paste0(strich,x));};
			if(max(summcols %in% group_by)) {
				while(max(sapply(summcols, pref) %in% group_by)) strich <- paste0('_',strich);
				summcols <- sapply(summcols, pref)
				warning(paste0("One or more of the columns ",paste(summcols, sep=', ')," is currently occupied. Summary columns will be appropriately renamed. In future you may wish to rename your columns before using the `clusterby` methods."));
			}
			tib_summ <- as_tibble();
			for(col in c(group_by,summcols)) tib_summ <- tib_summ %>% add_column(!!(col):=c(NA));
			tib_summ <- tib_summ[c(), ];
		}

		i0 <- 1;
		cl <- 0;
		while(i0 <= n) {
			i1 <- i0 + 1;
			pt <- tib[i0, by][[1]][[1]];
			g <- tib[i0, groupname][[1]][[1]]
			pt_ <- pt;
			while(i1 <= n) {
				pt__ <- tib[i1, by][[1]][[1]];
				g_ <- tib[i1, groupname][[1]][[1]]
				if(!near(pt, pt__) || !(g == g_)) {
					i1 <- i1 - 1;
					break;
				}
				pt_ <- pt__;
				i1 <- i1 +  1;
			}
			sz <- i1-i0+1;
			if(sz >= min_cluster_size) {
				if(do_summary) {
					cl <- cl + 1;
					for(col in group_by) tib_summ[cl, col] <- tib[i0, col];
					tib_summ[cl, summcols] <- c(pt, pt_, i0, i1, sz);
				} else {
					for(j in c(i0,i1)) tib[j, clustername] <- cl;
					cl <- cl + 1;
				}
			}
			i0 <- i1 + 1;
		}

		if(do_summary) {
			data <- tib_summ;
		} else {
			data <- tib %>% select(-c(groupname));
		}
	} else {
		edgesname <- uniquecolumnname('edges', c(cols, clustername));

		## Prägruppierung der Daten:
		leer <- list(); for(i in c(1:n)) leer[[i]] <- c(NA);
		tib <- tib %>% add_column(!!(clustername):=leer, !!(edgesname):=leer);
		if(is_linear) tib <- tib[order(tib[[by]]), ];

		# Erzeuge Kanten für Klusterbausteine.
		if(length(group_by) > 0) tib <- tib %>% group_by_at(group_by);
		tib <- tib %>% group_by_at(group_by) %>% nest(.key=!!(chunkname));
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
		data <- tib %>% add_column(!!(clustername) := clusters);

		if(split) data <- group_by_at(data, c(group_by, clustername)); #%>% nest(.key=!!(dataname);
	}

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

uniquecolumnname <- function(nom, cols) {
	col <- nom;
	while(col %in% cols) col <- paste0('_',col);
	return(col);
};