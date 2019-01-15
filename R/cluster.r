#' Data-clustering
#'
#' This package contains methods, which enables clustering in dataframes. Particularly useful for bio-mathematics, cognitive sciences, etc.
#'
#' \code{cluster(df, ...)}
#' @param df Tibble/Dataframe to be clustered. Method also possible with vectors.
#' @param by string vector. Specifies the column(s) for geometric data, according to which the clusters are to be built.
#' @param filter.by string vector. Defaults to \code{c()}. Specificies columns, by which data is to be preliminarily divided into groups, within which the clusters are to be built.
#' @param keep string vector. Defaults to \code{c()}. Specificies columns, which should be kept in a summary.
#' @param near symmetric function in two arguments. This function operates pairs of entries in the columns with geometric data and returns \code{TRUE}/\code{FALSE} if entries are near. Defaults to a Manhattan metric.
#' @param min.dist a real number. Defaults to \code{0}. If the default manhattan metric is used for \code{near}, this is the minimum tolerated distance between geometric data.
#' @param max.dist a real number. Defaults to \code{Inf}. If the default manhattan metric is used for \code{near}, this is the maximum tolerated distance between geometric data.
#' @param strict boolean. Defaults to \code{FALSE}. If the default manhattan metric is used for \code{near}, this sets the proximity to be a strict \code{< dist} or else \code{<= dist}.
#' @param cluster.name string. Defaults to \code{'cluster'}. Running \code{df \%>\% clusterby(...)} returns a data frame, which extends \code{df} by 1 column with this name. This column tags the clusters by a unique index.
#' @param min.size a natural number. Defaults to \code{0}. If a cluster has fewer elements as this, it will not be viewed as a cluster.
#' @param split boolean. Defaults to \code{FALSE}. If set to \code{TRUE}, then the output will be group the tibble data by cluster (equivalent to performing \code{\%>\% group_by(...)}).
#' @param is.linear boolean. Defaults to \code{TRUE}, if \code{length(by)=1}, otherwise to \code{FALSE}. If set to \code{TRUE}, then the geometry is assumed to be linear and endowed with a simple difference-metric. This allows for faster computation.
#' @param is.disjoint boolean. Defaults to \code{TRUE}, if \code{length(by)=1}, otherwise to \code{FALSE}. If set to \code{TRUE} in combination with \code{is.linear=TRUE}, then the clusters must occupy disjoint intervals.
#' @param presort boolean. Defaults to \code{TRUE}. If \code{TRUE} and \code{is.linear=TRUE}, then the data will be sorted by the \code{by} column first before the clusters are built. Normally this is desirable, however the option is included to prevent this, by setting \code{presort=FALSE}.
#' @param summary boolean. Defaults to \code{FALSE}. If set to \code{TRUE} in combination with \code{is.linear=TRUE} and **assuming** the user has presorted the data by the \code{by}-column, then a summary of the clusters as intervalls is provided. This makes most sense, if \code{is.disjoint=TRUE}. This produces the columns \code{filter.by, by, pstart, pend, nstart, nend, n} where \code{pstart}, \code{pend} describes the interval, \code{nstart}, \code{nend} provides the original indices in the input data, and \code{n} is the cluster size (number of points).
#' @export cluster
#' @examples gene %>% cluster(by='position', filter.by=c('gene','actve'), min.size=4, max.dist=400, strict=FALSE, is.linear=TRUE, is.disjoint=TRUE, summary=TRUE);
#' @examples protein3d %>% cluster(by=c('x','y','z'), filter.by='celltype', max.dist=5.8e-7, cluster.name='segment');
#' @examples soil_data %>% cluster(by=c('x','y'), filter.by=c('density','substance'), max.dist=10e-3, cluster.name='clump');
#' @keywords cluster clustering gene




cluster <- function(data, ...) {
	INPUTVARS <- list(...);
	VARNAMES <- names(INPUTVARS);

	group_by <- INPUTVARS[['filter.by']];
	keep <- INPUTVARS[['keep']];
	by <- INPUTVARS[['by']];
	clustername <- INPUTVARS[['cluster.name']];
	min_cluster_size <- INPUTVARS[['min.size']];
	split <- INPUTVARS[['split']];
	near <- INPUTVARS[['near']];
	d_min <- INPUTVARS[['min.dist']];
	d_max <- INPUTVARS[['max.dist']];
	strict <- INPUTVARS[['strict']];
	is_linear <- INPUTVARS[['is.linear']];
	is_disjoint <- INPUTVARS[['is.disjoint']];
	presort <- INPUTVARS[['presort']];
	do_summary <- INPUTVARS[['summary']];

	if(!is.vector(group_by)) group_by <- c();
	if(!is.vector(keep)) keep <- c();
	if(length(by) == 1) {
		if(!is.logical(is_linear)) is_linear <- TRUE;
		if(!is.logical(is_disjoint)) is_disjoint <- TRUE;
		if(!is.logical(do_summary) || !is_disjoint) do_summary <- FALSE;
		if(!is.logical(presort)) presort <- TRUE;
	} else {
		is_linear <- FALSE;
		is_disjoint <- FALSE;
		do_summary <- FALSE;
		presort <- FALSE;
	}
	if(!is.numeric(min_cluster_size)) min_cluster_size <- 0;
	if(!is.numeric(d_min)) d_min <- 0;
	if(!is.numeric(d_max)) d_max <- Inf;
	if(!is.logical(strict)) strict <- FALSE;
	if(!is.function(near)) {
		if(!is.character(near)) near <- 'Manhattan';
		if(near == 'Euclidean') {
			if(strict) {
				near <- function(x, y) {d <- sqrt(sum((x-y)^2)); return(d_min <= d && d < d_max);};
			} else {
				near <- function(x, y) {d <- sqrt(sum((x-y)^2)); return(d_min <= d && d <= d_max);};
			}
		} else {#if(near == 'Manhattan') {
			if(strict) {
				near <- function(x, y) {d <- max(abs(x-y)); return(d_min <= d && d < d_max);};
			} else {
				near <- function(x, y) {d <- max(abs(x-y)); return(d_min <= d && d <= d_max);};
			}
		}
	}
	if(!is.logical(split) || !do_summary) split <- FALSE;


	## Erstellung von Spaltennamen (Klusterspalte + Pufferspalte):
	tib <- tibble::as_tibble(data);
	cols <- names(tib);
	n <- nrow(tib);
	if(!is.character(clustername)) clustername <- uniquecolumnname('cluster', cols);
	chunkname <- uniquecolumnname('chunk', c(cols, clustername));


	if(is_linear && is_disjoint) {
		## Füge Indizes für Gruppennamen und zur Erhaltung der Reihenfolge hinzu:
		if(presort) tib <- tib[order(tib[[by]]), ];
		indexname <- uniquecolumnname('index', c(cols, clustername));
		tib <- tib %>% add_column(!!(indexname):=c(1:n));
	}


	## Indizes für Präsortierung hinzufügen.
	groupname <- uniquecolumnname('group', c(cols, clustername));
	tib <- tib %>% group_by_at(group_by) %>% nest(.key=!!(chunkname));
	m <- nrow(tib);
	tib <- tib %>% add_column(!!(groupname):=c(1:m));
	tib <- tib %>% unnest();


	if(is_linear) {
		if(is_disjoint) {
			## Reihenfolge wiederherstellen.
			tib <- tib[order(tib[[indexname]]), ] %>% select(-c(indexname));
		} else if(presort) {
			tib <- tib[order(tib[[groupname]],tib[[by]]), ];
		} else {
			tib <- tib[order(tib[[groupname]]), ];
		}

		if(is_disjoint && do_summary) {
			oldcols <- unique(c(group_by, keep));
			summcols <- c('pstart','pend','nstart','nend','distance','n');
			strich <- '';
			pref <- function(x) {return(paste0(strich,x));};
			if(max(summcols %in% oldcols)) {
				while(max(sapply(summcols, pref) %in% oldcols)) strich <- paste0('_',strich);
				summcols <- sapply(summcols, pref)
				warning(paste0("One or more of the columns ",paste(summcols, sep=', ')," is currently occupied. Summary columns will be appropriately renamed. In future you may wish to rename your columns before using the `clusterby` methods."));
			}
			tib_summ <- as_tibble();
			for(col in c(oldcols,summcols)) tib_summ <- tib_summ %>% add_column(!!(col):=c(NA));
			tib_summ <- tib_summ[c(), ];
		}

		i0 <- 1;
		cl <- 0;
		sel <- c();
		while(i0 <= n) {
			i1 <- i0 + 1;
			pt <- tib[i0, by][[1]][[1]];
			g <- tib[i0, groupname][[1]][[1]]
			pt_ <- pt;
			while(i1 <= n) {
				g_ <- tib[i1, groupname][[1]][[1]]
				pt__ <- tib[i1, by][[1]][[1]];
				if(!(g == g_) || !near(pt, pt__)) break;
				pt_ <- pt__;
				i1 <- i1 +  1;
			}
			i1 <- i1 - 1;

			sz <- i1-i0+1;
			if(sz >= min_cluster_size) {
				if(do_summary) {
					cl <- cl + 1;
					for(col in oldcols) tib_summ[cl, col] <- tib[i0, col];
					tib_summ[cl, summcols] <- c(pt, pt_, i0, i1, abs(pt_-pt) + 1, sz);
				} else {
					dsel <- c(i0:i1);
					sel <- c(sel, dsel);
					tib[dsel, clustername] <- cl;
					cl <- cl + 1;
				}
			}
			i0 <- i1 + 1;
		}

		if(do_summary) {
			data <- tib_summ; # %>% arrange_(group_by);
		} else {
			data <- tib[sel, ] %>% select(-c(groupname));
		}
	} else {
		edgesname <- uniquecolumnname('edges', c(cols, clustername));

		## Prägruppierung der Daten:
		leer <- list(); for(i in c(1:n)) leer[[i]] <- c(NA);
		tib <- tib %>% add_column(!!(clustername):=leer, !!(edgesname):=leer);

		# Erzeuge Kanten für Klusterbausteine.
		tib <- tib %>% group_by_at(groupname) %>% nest(.key=!!(chunkname));
		pos0 <- 0;
		m <- nrow(tib);
		for(i in c(1:m)) {
			chunk <- tib[i, chunkname][[1]][[1]];
			n <- nrow(chunk);
			if(n < min_cluster_size) next;
			pts <- list(); for(j in c(1:n)) pts[[j]] <- chunk[j, by];
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
			chunk[[edgesname]] <- edges;
			tib[i, chunkname][[1]][[1]] <- chunk;
			pos0 <- pos0 + n;
		}

		# Erzeuge Kluster aus Kanten.
		tib <- tib %>% unnest();
		clusters <- generateclasses(tib[[edgesname]], min_cluster_size);
		ind <- which(!is.na(clusters));
		clusters <- clusters[ind];
		data <- tib[ind, ] %>%
			select(-c(groupname, edgesname)) %>%
			add_column(!!(clustername) := clusters);
	}

	if(split) data <- group_by_at(data, c(group_by, clustername)); #%>% nest(.key=!!(dataname);

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