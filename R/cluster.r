#' Data-clustering
#'
#' This package contains methods, which enables clustering in dataframes. Particularly useful for bio-mathematics, cognitive sciences, etc.
#'
#' \code{cd <- clusterby::clusterdataframe(df)}
#' \code{cd$build(...)}
#' \code{cd$summarise(...)}
#' \code{cd$getoriginal()}
#' \code{cd$getclusters(...)}
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
#' @export clusterdataframe
#' @examples cdf <- clusterby::clusterdataframe(gene); cdf$build(by='position', filter.by=c('gene','actve'), min.size=4, max.dist=400, strict=TRUE, is.lexical=TRUE, is.disjoint=TRUE);
#' @examples cdf <- clusterby::clusterdataframe(protein3d); cdf$build(by=c('x','y','z'), filter.by='celltype', max.dist=5.8e-7, cluster.name='segment');
#' @examples cdf <- clusterby::clusterdataframe(soil_data); cdf$build(by=c('x','y'), filter.by=c('density','substance'), max.dist=10e-3, cluster.name='clump');
#' @examples data <- cdf$getclusters();
#' @examples tib <- cdf$getclusters();
#' @examples tib <- cdf$getclusters(summary=FALSE);
#' @examples tib_summ <- cdf$getclusters(summary=TRUE);
#' @keywords cluster clustering gene
# #' @param tib tibble data which has been grouped.
# #' @examples tib %>% clusterby::clusterwise(concentration=mean, names=c('set',';'), attributes=c('list',';'));
# #' @export clusterwise
# #' @keywords cluster summarise



clusterdataframe <- setRefClass('clusterdataframe',
	fields = list(
		isclustered='logical',
		cluster.name='character',
		cluster.group.by='ANY',
		cluster.by='ANY',
		cluster.keep='ANY',
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

			.self[['cluster.name']] <- 'cluster';
			for(key in c(
				'cluster.by',
				'cluster.group.by',
				'cluster.keep'
			)) .self[[key]] <- c();
			for(key in c(
				'cluster.data',
				'cluster.summary'
			)) .self[[key]] <- tibble::as_tibble(list());
			.self$isclustered <- FALSE;
		},
		get = function(key, ...) {
			if(key %in% c(
				'cluster.name',
				'cluster.by',
				'cluster.group.by',
				'cluster.keep'
			)) return(.self[[key]]);

			if(key == 'original') return(.self[['data']]);

			if(key == 'clusters') {
				if(!.self$isclustered) return(tibble::as_tibble(list()));

				INPUTVARS = list(...);
				summary <- INPUTVARS[['summary']];
				if(!is.logical(summary)) summary <- FALSE;
				if(summary) return(.self[['cluster.summary']]);
				return(.self[['cluster.data']])
			}

			return(NULL);
		},
		groupby = function() {
			if(!.self$isclustered) return(tibble::as_tibble(list()));
			return(.self[['cluster.data']] %>% group_by_at(.self[['cluster.name']]));
		},
		build = function(...) {
			obj <- .self$data %>% buildclusters____(...);
			for(key in c(
				'cluster.name',
				'cluster.by',
				'cluster.group.by',
				'cluster.keep',
				'cluster.data',
				'cluster.summary'
			)) .self[[key]] <- obj[[key]];
			.self$isclustered <- TRUE;
		},
		summarise = function(...) {
			if(!.self$isclustered) return(tibble::as_tibble(list()));

			INPUTVARS <- list(...);
			summcols <- names(INPUTVARS);

			tib <- .self[['cluster.data']];
			clustername <- .self[['cluster.name']];
			group_by <- .self[['cluster.group.by']];
			keep <- .self[['cluster.keep']];
			# cols <- names(tib);
			# if(length(keep) == 0) keep <- cols;
			keep <- unique(c(clustername, group_by, keep));
			keep <- keep[which(!(keep %in% summcols))];

			instructions <- list();
			for(col in keep) instructions[[col]] <- list(
				'col'=col,
				'method'='pick'
			);

			for(col in summcols) {
				s <- INPUTVARS[[col]];
				instructions[[col]] <- s;
				if(is.list(s)) next;
				instructions[[col]] <- list(
					'col'=col,
					'method'=s
				);
			}

			method <- list();
			k <- 1;
			fnames <- c();
			colnames <- c();

			for(col in c(keep,summcols)) {
				if(col == clustername) next;

				s <- instructions[[col]];
				m <- s[['method']];

				if(is.function(m)) {
					f <- opt;
				} else if(is.character(m)) {
					opt <- m[1];
					if(opt == 'pick') {
						f <- function(x) {return(x[1]);};
					} else if(opt %in% c('set','list')) {
						sep <- ',';
						if('sep' %in% names(s)) {
							sep <- s[['sep']];
						} else if(length(m) > 1) {
							sep <- m[2];
						}
						if(opt == 'set') {
							f <- function(x) {return(paste0('[',paste(unique(x), collapse=sep),']'));};
						} else if(opt == 'list') {
							f <- function(x) {return(paste0('[',paste(x, collapse=sep),']'));};
						}
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
						# f <- function(x) {return(NA);};
						next;
					}
				} else {
					next;
				}

				method[[k]] <- f;
				vnames <- paste(s[['col']], collapse=',');
				fnames[k] <- paste0('method[[',k,']](',vnames,')');
				colnames[k] <- col;
				k <- k + 1;
			}

			tib_summ <- .self$groupby() %>% dplyr::summarise_(.dots=setNames(fnames, colnames));

			return(tib_summ);
		}
	)
);




buildclusters____ <- function(tib, ...) {
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

	tib <- tibble::as_tibble(tib);
	cols <- names(tib);

	dim_by <- length(by);
	if(!is.vector(group_by)) group_by <- c();
	if(!is.vector(keep)) keep <- c();
	keep_orig <- keep;
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
	keep_summ <- unique(c(group_by, keep));
	if(!('keep' %in% VARNAMES)) keep <- cols;
	keep <- unique(c(group_by, by, keep));


	## Erstellung von Spaltennamen (Klusterspalte + Pufferspalte):
	n <- nrow(tib);
	if(!is.character(clustername)) clustername <- uniquecolumnname____('cluster', cols);
	chunkname <- uniquecolumnname____('chunk', c(cols, clustername));
	indexname <- uniquecolumnname____('index', c(cols, clustername));


	## Definiere Typen von Spalten
	types <- list(
		data=list(),
		summary=list()
	);
	if(as_interval) {
		types$summary <- list(
			'nstart'='integer',
			'nend'='integer',
			'size'='integer',
			'pstart'='double',
			'pend'='double',
			'distance'='double'
		)
	} else if(is_lexical) {
		types$summary <- list(
			'indices'='character',
			'size'='integer',
			'positions'='character'
		);
	} else {
		types$summary <- list(
			'indices'='character',
			'size'='integer'
		);
	}
	summcols <- names(types$summary);
	types$data[[clustername]] <- 'integer';
	types$summary[[clustername]] <- 'integer';


	## Indexnamen bevor Verarbeitung speichern.
	tib <- tib %>% tibble::add_column(!!(indexname):=c(1:n), !!(clustername):=rep(0, n));
	## Indizes für Präsortierung hinzufügen.
	groupname <- uniquecolumnname____('group', c(cols, clustername));
	tib <- tib %>% group_by_at(group_by) %>% nest(.key=!!(chunkname));
	Ng <- nrow(tib);
	tib <- tib %>% tibble::add_column(!!(groupname):=c(1:Ng));
	tib <- tib %>% unnest();
	## ursp. Reihenfolge wiederherstellen.
	tib <- tib %>% dplyr::arrange_(indexname);


	## Initialisiere Summary-Tibble
	strich <- '';
	pref <- function(x) {return(paste0(strich, x));};
	if(max(summcols %in% keep_summ)) {
		while(max(sapply(summcols, pref) %in% keep_summ)) strich <- paste0('_',strich);
		summcols <- sapply(summcols, pref)
		warning(paste0("One or more of the columns ",paste(summcols, sep=', ')," is currently occupied. Summary columns will be appropriately renamed to ",paste(summcols, sep=', '),". In future you may wish to rename your columns before using the `clusterby` methods."));
	}
	tib_summ <- (tib %>% dplyr::select(keep_summ))[c(1), ];
	for(col in c(clustername, summcols)) tib_summ <- tib_summ %>% tibble::add_column(!!(col):=c(NA));
	tib_summ <- tib_summ[c(), ];

	## HAUPTMETHODE
	if(is_lexical) {
		if(presort) tib <- tib %>% dplyr::arrange_(by);
		if(!is_disjoint) tib <- tib %>% dplyr::arrange_(groupname);
		tib_ <- tib;

		i0 <- 1;
		cl <- 1;
		sel <- c();
		pt_to_json <- function(x) {
			vals <- sapply(x, function(val) {
				if(is.character(val)) val <- paste0('"',val,'"');
				return(val);
			});
			return(paste0('[',paste(vals, collapse=','),']'));
		};
		while(i0 <= n) {
			## Iteriere bis Gruppe sich ändert, oder nächster Punkt zu weit weg vom Kluster liegt.
			i1 <- i0 + 1;
			g <- tib[i0, groupname][[1]];
			pt <- tib[i0, by];
			pt_ <- pt;
			while(i1 <= n) {
				g_ <- tib[i1, groupname][[1]];
				pt__ <- tib[i1, by];
				if(!(g == g_) || !near(pt_, pt__)) break;
				pt_ <- pt__;
				i1 <- i1 +  1;
			}
			i1 <- i1 - 1;

			## Fasse als Kluster zusammen.
			sz <- i1-i0+1;
			if(sz >= min_cluster_size) {
				tib_summ[cl, clustername] <- cl;
				for(col in keep) tib_summ[cl, col] <- tib[i0, col];
				if(as_interval) {
					dp <- metric_lex(pt, pt_);
					if(dim_by > 1) {
						pt <- pt_to_json(tib[i0, by]);
						pt_ <- pt_to_json(tib[i1, by]);
					}
					tib_summ[cl, summcols[1]] <- i0;
					tib_summ[cl, summcols[2]] <- i1;
					tib_summ[cl, summcols[3]] <- sz;
					tib_summ[cl, summcols[4]] <- pt;
					tib_summ[cl, summcols[5]] <- pt_;
					tib_summ[cl, summcols[6]] <- dp;
				} else {
					idx <- tib[c(i0:i1), indexname][[1]];
					s <- order(idx);
					idx <- paste0('[',paste(idx[s], collapse=','),']');
					if(dim_by == 1) {
						pts <- pt_to_json(tib[c(i0:i1)[s], by]);
					} else {
						pts <- sapply(c(i0:i1)[s], function(i) {return(pt_to_json(tib[i, by]));});
						pts <- paste0('[',paste(pts, collapse=','),']');
					}
					tib_summ[cl, summcols[1]] <- idx;
					tib_summ[cl, summcols[2]] <- sz;
					tib_summ[cl, summcols[3]] <- pts;
				}

				dsel <- c(i0:i1);
				sel <- c(sel, dsel);
				tib[dsel, clustername] <- cl;

				cl <- cl + 1;
			}
			i0 <- i1 + 1;
		}
	} else {
		tib <- tib %>% dplyr::arrange_(groupname);
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
				tib_summ[cl_, clustername] <- cl_;
				for(col in keep) tib_summ[cl_, col] <- chunk[i0, col];
				idx <- sort(chunk[dsel, indexname][[1]]);
				idx <- paste0('[',paste(idx, collapse=','),']');
				sz <- length(dsel);
				tib_summ[cl_, summcols[1]] <- idx;
				tib_summ[cl_, summcols[2]] <- sz;
			}

			dsel <- which(!is.na(clusters));
			sel <- c(sel, ind[dsel]);
			tib[ind, clustername] <- clusters;
		}
	}


	## Reinige Typen.
	for(col in names(types$data)) {
		typ <- types$data[[col]];
		if(typ == 'integer') tib[ , col] <- as.integer(tib[[col]]);
	}
	for(col in names(types$summary)) {
		typ <- types$summary[[col]];
		if(typ == 'integer') tib_summ[ , col] <- as.integer(tib_summ[[col]]);
	}


	## Reihenfolge wiederherstellen und Spalten filtrieren.
	# tib_summ <- tib_summ %>% dplyr::arrange_(group_by) %>% dplyr::select(keep, summcols));
	tib <- tib[sel, ] %>% dplyr::arrange_(indexname) %>% dplyr::select(c(clustername, keep));
	tib_summ <- tib_summ %>% dplyr::select(clustername, keep_summ, summcols);
	if(!as_interval) {
		names(tib_summ)[which(names(tib_summ)==summcols[1])] <- indexname;
		if(is_lexical && dim_by == 1) names(tib_summ)[which(names(tib_summ)==summcols[3])] <- by;
	}

	return(list(
		cluster.name=clustername,
		cluster.by=by,
		cluster.group.by=group_by,
		cluster.keep=keep_orig,
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