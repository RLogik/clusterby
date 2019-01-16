#' Data-clustering
#'
#' This package contains methods, which enables clustering in dataframes. Particularly useful for bio-mathematics, cognitive sciences, etc.
#'
#' \code{cd <- clusterby::clustereddata(df)}
#' \code{cd$buildclusters(...)}
#' \code{cd$getclusters(...)}
#' \code{cd$summarise(...)}
#'
#' @param df Tibble/Dataframe to be clustered. Method also possible with vectors.
#' @param by string vector. Specifies the column(s) for geometric data, according to which the clusters are to be built.
#' @param filter.by string vector. Defaults to \code{c()}. Specificies columns, by which data is to be preliminarily divided into groups, within which the clusters are to be built.
#' @param keep string vector. Defaults to \code{c()}. Specificies columns, which should be kept in a summary.
#' @param near symmetric function in two arguments. This function operates pairs of entries in the columns with geometric data and returns \code{TRUE}/\code{FALSE} if entries are near. Defaults to a Manhattan metric.
#' @param min.dist a real number. Defaults to \code{0}. If the default manhattan metric is used for \code{near}, this is the minimum tolerated distance between geometric data.
#' @param max.dist a real number. Defaults to \code{Inf}. If the default manhattan metric is used for \code{near}, this is the maximum tolerated distance between geometric data.
#' @param strict boolean. Defaults to \code{FALSE}. If the default manhattan metric is used for \code{near}, this sets the proximity to be a strict \code{< dist} or else \code{<= dist}.
#' @param cluster.name string. Defaults to \code{'cluster'}. Running \code{df \%>\% clusterby(...)} returns a data frame, which extends \code{df} by 1 column with this name. This column tags the clusters by a unique index.
#' @param min.size a natural number. Defaults to \code{1}. If a cluster has fewer elements as this, it will not be viewed as a cluster.
#' @param split boolean. Defaults to \code{FALSE}. If set to \code{TRUE}, then the output will be group the tibble data by cluster (equivalent to performing \code{\%>\% group_by(...)}).
#' @param is.lexical boolean. Defaults to \code{TRUE} if \code{length(by)=1}, otherwise to \code{FALSE}. If set to \code{TRUE}, then the geometry is assumed to be linear and endowed with a simple difference-metric. This allows for faster computation.
#' @param is.disjoint boolean. Defaults to \code{TRUE} if \code{is.lexical=TRUE}, otherwise to \code{FALSE}. If set to \code{TRUE} in combination with \code{is.lexical=TRUE}, then the clusters must occupy disjoint intervals.
#' @param presort boolean. Defaults to \code{TRUE} if \code{is.lexical=TRUE}, otherwise to \code{FALSE}. If \code{TRUE} and \code{is.lexical=TRUE}, then the data will be sorted by the \code{by} column first before the clusters are built. Normally this is desirable, however the option is included to prevent this, by setting \code{presort=FALSE}. Regardless of this setting, in the non-summary mode, all data will be reordered to correspond to the sequence of the input data.
#' @param summary boolean. Defaults to \code{FALSE}. If set to \code{TRUE} in combination with \code{is.lexical=TRUE} and **assuming** the user has presorted the data by the \code{by}-column, then a summary of the clusters as intervalls is provided. This makes most sense, if \code{is.disjoint=TRUE}. This produces the columns \code{filter.by, by, pstart, pend, nstart, nend, n} where \code{pstart}, \code{pend} describes the interval, \code{nstart}, \code{nend} provides the original indices in the input data, and \code{n} is the cluster size (number of points).
#' @param as.interval boolean. Defaults to \code{TRUE} if \code{is.lexical=TRUE} and \code{is.disjoint=TRUE}, otherwise defaults to \code{FALSE}. If \code{TRUE} and, then summaries provide information as interval end points. If \code{FALSE}, then summaries are provided as lists.
#' @export clustereddata
#' @examples cd <- clusterby::clustereddata(gene); cd$buildclusters(by='position', filter.by=c('gene','actve'), min.size=4, max.dist=400, strict=TRUE, is.lexical=TRUE, is.disjoint=TRUE);
#' @examples cd <- clusterby::clustereddata(protein3d); cd$buildclusters(by=c('x','y','z'), filter.by='celltype', max.dist=5.8e-7, cluster.name='segment');
#' @examples cd <- clusterby::clustereddata(soil_data); cd$buildclusters(by=c('x','y'), filter.by=c('density','substance'), max.dist=10e-3, cluster.name='clump');
#' @keywords cluster clustering gene
# #' @param tib tibble data which has been grouped.
# #' @examples tib %>% clusterby::clusterwise(concentration=mean, names=c('set',';'), attributes=c('list',';'));
# #' @export clusterwise
# #' @keywords cluster summarise



clustereddata <- setRefClass('clustereddata',
	fields = list(
		isclustered='logical',
		cluster.name='character',
		cluster.group.by='ANY',
		data='ANY',
		cluster.data='ANY',
		cluster.summary='ANY'
	),
	methods = list(
		initialize = function(...) {
			INPUTVARS = list(...);
			if(length(INPUTVARS) >= 1) {
				.self$data <- INPUTVARS[[1]];
			} else {
				.self$data <- tibble::as_tibble(list());
			}

			.self[['cluster.data']] <- tibble::as_tibble(list());
			.self[['cluster.summary']] <- tibble::as_tibble(list());
			.self[['cluster.name']] <- 'cluster';
			.self[['cluster.group.by']] <- c();
			.self$isclustered <- FALSE;
		},
		buildclusters = function(...) {
			obj <- .self$data %>% buildclusters____(...);
			.self[['cluster.name']] <- obj[['cluster.name']];
			.self[['cluster.group.by']] <- obj[['cluster.group.by']];
			.self[['cluster.data']] <- obj[['cluster.data']];
			.self[['cluster.summary']] <- obj[['cluster.summary']];
			.self$isclustered <- TRUE;
		},
		getclusters = function(...) {
			if(!.self$isclustered) return(tibble::as_tibble(list()));

			INPUTVARS = list(...);
			summary <- INPUTVARS[['summary']];
			if(!is.logical(summary)) summary <- FALSE;
			if(summary) return(.self[['cluster.summary']]);
			return(.self[['cluster.data']])
		},
		summarise = function(...) {
			if(!.self$isclustered) return(tibble::as_tibble(list()));

			INPUTVARS <- list(...);
			summcols <- names(INPUTVARS);

			tib <- .self[['cluster.data']];
			cols <- names(tib);
			summcols <- summcols[which(summcols %in% cols)];

			method <- list();
			for(col in summcols) {
				s <- INPUTVARS[[col]];
				opt <- s[1];
				if(is.function(opt)) {
					f <- opt;
				} else if(is.character(opt)) {
					if(opt == 'set') {
						sep = ','; if(length(s) > 1) sep = s[2];
						f <- function(x) {return(paste0('[',paste(unique(x), collapse=','),']'));};
					} else if(opt == 'list') {
						sep = ','; if(length(s) > 1) sep = s[2];
						f <- function(x) {return(paste0('[',paste(x, collapse=','),']'));};
					} else if(opt == 'length') {
						f <- length;
					} else if(opt == 'min') {
						f <- function(x) {return(min(x, na.rm=TRUE));};
					} else if(opt == 'max') {
						f <- function(x) {return(max(x, na.rm=TRUE));};
					} else if(opt == 'range') {
						f <- function(x) {return(paste0('[',paste(c(min(x, na.rm=TRUE), max(x, na.rm=TRUE)), collapse=','),']'));};
					} else if(opt == 'mean') {
						f <- function(x) {return(mean(x, na.rm=TRUE));};
					} else if(opt == 'var') {
						f <- function(x) {return(var(x, na.rm=TRUE));};
					} else if(opt == 'sd') {
						f <- function(x) {return(sd(x, na.rm=TRUE));};
					} else {
						f <- function(x) {return(NA);};
					}
				}
				method[[col]] <- f;
			}

			clustername <- .self[['cluster.name']];
			group_by <- .self[['cluster.group.by']];
			tib_summ <- tibble::as_tibble();
			for(col in c(group_by, clustername, summcols)) tib_summ <- tib_summ %>% add_column(!!(col):=c(NA));
			tib_summ <- tib_summ[c(), ];

			clusterIds <- tib[[clustername]];
			j <- 1;
			for(cl in unique(clusterIds)) {
				ind <- which(tib[[clustername]] == cl);
				cluster <- tib[ind, ];
				i0 <- ind[1];
				for(col in c(group_by, clustername)) tib_summ[j, col] <- tib[i0, col];
				for(col in summcols) {
					f <- method[[col]];
					tib_summ[j, col] <- cluster[[col]] %>% f;
				}
				j <- j + 1 ;
			}
			return(tib_summ);
		}
	)
);




buildclusters____ <- function(data, ...) {
	INPUTVARS <- list(...);
	VARNAMES <- names(INPUTVARS);

	group_by <- INPUTVARS[['filter.by']];
	keep <- INPUTVARS[['keep']];
	by <- INPUTVARS[['by']];
	clustername <- INPUTVARS[['cluster.name']];
	min_cluster_size <- INPUTVARS[['min.size']];
	near <- INPUTVARS[['near']];
	d_min <- INPUTVARS[['min.dist']];
	d_max <- INPUTVARS[['max.dist']];
	strict <- INPUTVARS[['strict']];
	is_lexical <- INPUTVARS[['is.lexical']];
	is_disjoint <- INPUTVARS[['is.disjoint']];
	presort <- INPUTVARS[['presort']];
	as_interval <- INPUTVARS[['as.interval']];


	dim_by <- length(by);
	if(!is.vector(group_by)) group_by <- c();
	if(!is.vector(keep)) keep <- c();
	if(!is.logical(is_lexical)) is_lexical <- (dim_by == 1);
	if(!is.logical(is_disjoint)) is_disjoint <- is_lexical;
	if(!is.logical(presort)) presort <- is_lexical;
	if(!is.logical(as_interval)) as_interval <- (is_lexical && is_disjoint);
	if(!is.numeric(min_cluster_size)) min_cluster_size <- 1;
	if(!is.numeric(d_min)) d_min <- 0;
	if(!is.numeric(d_max)) d_max <- Inf;
	if(is_lexical) near <- 'lexical';
	if(!is.logical(strict)) strict <- FALSE;
	metric_lex <- function(x, y) {
		d <- abs(y-x);
		i <- min(which(d>0),dim_by);
		return(d[i])
	};
	if(!is.function(near)) {
		if(!is.character(near)) near <- 'Manhattan';
		if(near == 'lexical') {
			metric <- metric_lex;
		} else if(near == 'Euclidean') {
			metric <- function(x, y) {return(sqrt(sum((x-y)^2)))};
		} else {#if(near == 'Manhattan') {
			metric <- function(x, y) {return(max(abs(x-y)))};
		}
		if(strict) {
			near <- function(x, y) {d <- metric(x,y); return(d_min <= d && d < d_max);};
		} else {
			near <- function(x, y) {d <- metric(x,y); return(d_min <= d && d <= d_max);};
		}
	}


	## Definiere Keep-Spalten
	tib <- tibble::as_tibble(data);
	cols <- names(tib);
	if('keep' %in% VARNAMES) {
		keep_summ <- unique(c(group_by, keep));
		keep <- unique(c(group_by, by, keep));
	} else {
		keep_summ <- group_by;
		keep <- cols;
	}


	## Erstellung von Spaltennamen (Klusterspalte + Pufferspalte):
	n <- nrow(tib);
	if(!is.character(clustername)) clustername <- uniquecolumnname____('cluster', cols);
	chunkname <- uniquecolumnname____('chunk', c(cols, clustername));
	indexname <- uniquecolumnname____('index', c(cols, clustername));


	## Indexnamen bevor Verarbeitung speichern.
	tib <- tib %>% add_column(!!(indexname):=c(1:n), !!(clustername):=rep(0, n));
	## Indizes für Präsortierung hinzufügen.
	groupname <- uniquecolumnname____('group', c(cols, clustername));
	tib <- tib %>% group_by_at(group_by) %>% nest(.key=!!(chunkname));
	Ng <- nrow(tib);
	tib <- tib %>% add_column(!!(groupname):=c(1:Ng));
	tib <- tib %>% unnest();
	## ursp. Reihenfolge wiederherstellen.
	tib <- tib %>% arrange_(indexname);


	## Initialisiere Summary-Tibble
	if(as_interval) {
		summcols <- c('pstart','pend','distance','nstart','nend','size');
	} else {
		summcols <- c('indices','size');
		if(is_lexical) summcols <- c(summcols,'positions');
	}
	strich <- '';
	pref <- function(x) {return(paste0(strich, x));};
	if(max(summcols %in% keep)) {
		while(max(sapply(summcols, pref) %in% keep)) strich <- paste0('_',strich);
		summcols <- sapply(summcols, pref)
		warning(paste0("One or more of the columns ",paste(summcols, sep=', ')," is currently occupied. Summary columns will be appropriately renamed to ",paste(summcols, sep=', '),". In future you may wish to rename your columns before using the `clusterby` methods."));
	}
	tib_summ <- tibble::as_tibble();
	for(col in c(keep, summcols)) tib_summ <- tib_summ %>% add_column(!!(col):=c(NA));
	tib_summ <- tib_summ[c(), ];


	## HAUPTMETHODE
	if(is_lexical) {
		if(presort) tib <- tib %>% arrange_(by);
		if(!is_disjoint) tib <- tib %>% arrange_(groupname);
		tib_ <- tib;

		i0 <- 1;
		cl <- 1;
		sel <- c();
		pt_to_json <- function(i) {
			return(paste0('[',paste(tib[i, by], collapse=','),']'));
		};
		while(i0 <= n) {
			## Iteriere bis Gruppe sich ändert, oder nächster Punkt zu weit weg vom Kluster liegt.
			i1 <- i0 + 1;
			pt <- tib[i0, by][[1]];
			g <- tib[i0, groupname][[1]]
			pt_ <- pt;
			while(i1 <= n) {
				g_ <- tib[i1, groupname][[1]]
				pt__ <- tib[i1, by][[1]];
				if(!(g == g_) || !near(pt_, pt__)) break;
				pt_ <- pt__;
				i1 <- i1 +  1;
			}
			i1 <- i1 - 1;

			## Fasse als Kluster zusammen.
			sz <- i1-i0+1;
			if(sz >= min_cluster_size) {
				for(col in keep) tib_summ[cl, col] <- tib[i0, col];
				if(as_interval) {
					dp <- metric_lex(pt, pt_);
					if(dim_by > 1) {
						pt <- pt_to_json(i0);
						pt_ <- pt_to_json(i1);
					}
					tib_summ[cl, summcols] <- c(pt, pt_, dp, i0, i1, sz);
				} else {
					idx <- tib[c(i0:i1), indexname][[1]];
					s <- order(idx);
					if(dim_by == 1) {
						pts <- tib[c(i0:i1), by][[1]];
					} else {
						pts <- sapply(c(i0:i1), pt_to_json);
					}
					idx <- paste0('[',paste(idx[s], collapse=','),']');
					pts <- paste0('[',paste(pts[s], collapse=','),']');
					tib_summ[cl, summcols] <- c(idx, sz, pts);
				}

				dsel <- c(i0:i1);
				sel <- c(sel, dsel);
				tib[dsel, clustername] <- cl;

				cl <- cl + 1;
			}
			i0 <- i1 + 1;
		}
	} else {
		tib <- tib %>% arrange_(groupname);
		cl <- 1;
		sel <- c();
		for(g in c(1:Ng)) {
			ind <- which(tib[[groupname]] == g);
			chunk <- tib[ind, ];
			n <- nrow(chunk);
			if(n < min_cluster_size) next;

			## Graphen mit Kanten erstellen, wenn Punkte „nahe“ liegen.
			pts <- list(); for(j in c(1:n)) pts[[j]] <- chunk[j, by];
			edges <- lapply(c(1:n), function(j) {
				e <- c();
				if(j < n) {
					pt <- pts[[j]];
					bool <- lapply(pts, function(pt_) {
						return(near(pt, pt_));
					});
					e <- which(unlist(bool));
				}
				if(length(e) == 0) e <- c(NA);
				return(e);
			});

			## Erzeuge Kluster aus Kanten.
			clusters <- generateconnectedcomponents____(edges, min_cluster_size, cl);
			dsel <- which(!is.na(clusters));
			if(length(dsel) > 0) cl <- max(clusters[dsel]) + 1;

			## Fasse Kluster zusammen
			for(cl_ in unique(clusters)) {
				if(is.na(cl_)) next;
				dsel <- which(clusters==cl_);
				i0 <- dsel[1];
				for(col in keep) tib_summ[cl_, col] <- chunk[i0, col];
				idx <- sort(chunk[dsel, indexname][[1]]);
				idx <- paste0('[',paste(idx, collapse=','),']');
				sz <- length(dsel);
				tib_summ[cl_, summcols] <- c(idx, sz);
			}

			dsel <- which(!is.na(clusters));
			sel <- c(sel, ind[dsel]);
			tib[ind, clustername] <- clusters;
		}
	}


	## Reihenfolge wiederherstellen und Spalten filtrieren.
	# tib_summ <- tib_summ %>% arrange_(group_by) %>% select(c(keep, summcols));
	tib <- tib[sel, ] %>% arrange_(indexname) %>% select(c(keep, clustername));
	tib_summ <- tib_summ %>% select(c(keep, summcols));
	if(!as_interval) {
		names(tib_summ)[which(names(tib_summ)==summcols[1])] <- indexname;
		if(is_lexical) names(tib_summ)[which(names(tib_summ)==summcols[3])] <- by;
	}

	return(list(
		cluster.group.by=group_by,
		cluster.name=clustername,
		cluster.data=tib,
		cluster.summary=tib_summ
	));
};

generateconnectedcomponents____ <- function(edges, min_sz, key0) {
	n <- length(edges);
	classes <- rep(NA,n);
	if(n == 0) return(classes);

	key <- key0;
	ind <- c(1:n);
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

	return(classes);
};

uniquecolumnname____ <- function(nom, cols) {
	col <- nom;
	while(col %in% cols) col <- paste0('_',col);
	return(col);
};