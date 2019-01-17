#' Data-clustering
#'
#' This package contains methods, which enables clustering in dataframes. Particularly useful for bio-mathematics, cognitive sciences, etc.
#'
#' \code{cd <- clusterby::clusterdataframe(df)}
#' \code{cd$build(...)}
#' \code{cd$summarise(...)}
#' \code{cd$get('original', ...)}
#' \code{cd$get('clusters', summary=<lgl>, ...)}
#'
#' @param df Tibble/Dataframe to be clustered. Method also possible with vectors.
#' @param by string vector. Specifies the column(s) for geometric data, according to which the clusters are to be built.
#' @param filter.by string vector. Defaults to \code{c()}. Specificies columns, by which data is to be preliminarily divided into groups, within which the clusters are to be built.
#' @param keep string vector. Defaults to \code{c()}. Specificies columns, which should be kept when using the $get.
#' @param near symmetric function in two arguments. This function operates pairs of entries in the columns with geometric data and returns \code{TRUE}/\code{FALSE} if entries are near. Defaults to a Manhattan metric.
#' @param min.dist a real number. Defaults to \code{0}. If the default manhattan metric is used for \code{near}, this is the minimum tolerated distance between geometric data.
#' @param max.dist a real number. Defaults to \code{Inf}. If the default manhattan metric is used for \code{near}, this is the maximum tolerated distance between geometric data.
#' @param strict boolean. Defaults to \code{FALSE}. If the default manhattan metric is used for \code{near}, this sets the proximity to be a strict \code{< dist} or else \code{<= dist}.
#' @param cluster.name string. Defaults to \code{'cluster'}. Running \code{df \%>\% clusterby(...)} returns a data frame, which extends \code{df} by 1 column with this name. This column tags the clusters by a unique index.
#' @param min.size a natural number. Defaults to \code{1}. If a cluster has fewer elements as this, it will not be viewed as a cluster.
#' @param max.size a natural number. Defaults to \code{Inf}, determining the maximum allowable size of a cluster.
#' @param split boolean. Defaults to \code{FALSE}. If set to \code{TRUE}, then the output will be group the tibble data by cluster (equivalent to performing \code{\%>\% group_by(...)}).
#' @param is.lexical boolean. Defaults to \code{TRUE} if \code{length(by)=1}, otherwise to \code{FALSE}. If set to \code{TRUE}, then the geometry is assumed to be linear and endowed with a simple difference-metric. This allows for faster computation.
#' @param no.overlaps boolean. Defaults to \code{TRUE} if \code{is.lexical=TRUE}, otherwise to \code{FALSE}. If set to \code{TRUE} in combination with \code{is.lexical=TRUE}, then the clusters must occupy intervals that do not overlap.
#' @param presort boolean. Defaults to \code{TRUE} if \code{is.lexical=TRUE}, otherwise to \code{FALSE}. If \code{TRUE} and \code{is.lexical=TRUE}, then the data will be sorted by the \code{by} column first before the clusters are built. Normally this is desirable, however the option is included to prevent this, by setting \code{presort=FALSE}. Regardless of this setting, in the non-summary mode, all data will be reordered to correspond to the sequence of the input data.
#' @param summary boolean. Defaults to \code{FALSE}. If set to \code{TRUE} in combination with \code{is.lexical=TRUE} and **assuming** the user has presorted the data by the \code{by}-column, then a summary of the clusters as intervalls is provided. This makes most sense, if \code{no.overlaps=TRUE}. This produces the columns \code{filter.by, by, pstart, pend, nstart, nend, n} where \code{pstart}, \code{pend} describes the interval, \code{nstart}, \code{nend} provides the original indices in the input data, and \code{n} is the cluster size (number of points).
#' @param as.interval boolean. Defaults to \code{TRUE} if \code{is.lexical=TRUE} and \code{no.overlaps=TRUE}, otherwise defaults to \code{FALSE}. If \code{TRUE} and, then summaries provide information as interval end points. If \code{FALSE}, then summaries are provided as lists.
#' @export clusterdataframe
#' @examples cdf <- clusterby::clusterdataframe(gene); cdf$build(by='position', filter.by=c('gene','actve'), min.size=4, max.dist=400, strict=TRUE, is.lexical=TRUE, no.overlaps=TRUE);
#' @examples cdf <- clusterby::clusterdataframe(protein3d); cdf$build(by=c('x','y','z'), filter.by='celltype', max.dist=5.8e-7, cluster.name='segment');
#' @examples cdf <- clusterby::clusterdataframe(soil_data); cdf$build(by=c('x','y'), filter.by=c('density','substance'), max.dist=10e-3, cluster.name='clump');
#' @examples data <- cdf$get('clusters');
#' @examples tib <- cdf$get('clusters', keep=c('colour','age'));
#' @examples tib <- cdf$get('clusters', summary=FALSE);
#' @examples tib_summ <- cdf$get('clusters', summary=TRUE, as.interval=TRUE);
#' @examples tib_summ <- cdf$get('clusters', summary=TRUE, as.interval=FALSE);
#' @keywords cluster clustering gene
# #' @param tib tibble data which has been grouped.
# #' @examples tib %>% clusterby::clusterwise(concentration=mean, names=c('set',';'), attributes=c('list',';'));
# #' @export clusterwise
# #' @keywords cluster summarise



clusterdataframe <- setRefClass('clusterdataframe',
	fields = list(
		is.clustered='logical',
		is.lexical='logical',
		index.name='character',
		cluster.name='character',
		cluster.loc='ANY',
		cluster.group.by='ANY',
		cluster.by='ANY',
		cols.keep='ANY',
		cols.keep.set='logical',
		data='ANY',
		data.cols='ANY',
		cluster.data='ANY'
	),
	methods = list(
		initialize = function(...) {
			INPUTVARS = list(...);
			if(length(INPUTVARS) >= 1) {
				.self$data <- INPUTVARS[[1]];
				.self$data.cols <- names(.self$data);
			} else {
				.self$data <- tibble::as_tibble(list());
				.self$data.cols <- c();
			}

			.self$index.name <- 'index';
			.self$cluster.name <- 'cluster';
			for(key in c(
				'cluster.loc',
				'cluster.by',
				'cluster.group.by',
				'cols.keep'
			)) .self[[key]] <- c();
			.self$cols.keep.set <- FALSE;
			.self$cluster.data <- tibble::as_tibble(list());
			.self$is.clustered <- FALSE;
			.self$is.lexical <- FALSE;
		},
		get = function(key, ...) {
			if(key %in% c(
				'is.clustered',
				'is.lexical',
				'index.name',
				'cluster.name',
				'cluster.loc',
				'cluster.by',
				'cluster.group.by',
				'cols.keep',
				'cols.keep.set'
			)) return(.self[[key]]);

			if(key == 'original') return(.self$data);

			if(key == 'clusters') {
				if(!.self$is.clustered) return(tibble::as_tibble(list()));

				INPUTVARS = list(...);
				summary <- INPUTVARS[['summary']];
				if(!is.logical(summary)) summary <- FALSE;

				by <- .self$cluster.by;
				indexname <- .self$index.name;
				keep <- INPUTVARS[['cols.keep']];
				if(is.character(keep)) .self$setkeep(keep);

				if(summary) {
					as_interval <- INPUTVARS[['as.interval']];
					with_pts <- INPUTVARS[['with.points']];
					sep <- INPUTVARS[['sep']];
					with_keys <- INPUTVARS[['with.keys']];
					with_braces <- INPUTVARS[['with.braces']];
					if(!is.logical(as_interval)) as_interval <- .self$is.lexical;
					if(!is.logical(with_pts)) with_pts <- (length(by) == 1);
					if(!is.logical(sep)) sep <- ';';
					if(!is.logical(with_keys)) with_keys <- FALSE;
					if(!is.logical(with_braces)) with_braces <- FALSE;

					if(as_interval) {
						if(.self$is.lexical && with_pts) {
							tib <- .self$summarise(
								'n.start' = list(col=indexname, method='min'),
								'n.end' = list(col=indexname, method='max'),
								'size' = list(col=indexname, method='length'),
								'p.start' = list(col=by, method='lex:min', sep=sep, with.keys=with_keys, with.braces=with_braces),
								'p.end' = list(col=by, method='lex:max', sep=sep, with.keys=with_keys, with.braces=with_braces),
								'distance' = list(col=by, method='lex:range', mode='distance', sep=sep, with.keys=with_keys, with.braces=with_braces)
							);
						} else {
							tib <- .self$summarise(
								'n.start' = list(col=indexname, method='min'),
								'n.end' = list(col=indexname, method='max'),
								'size' = list(col=indexname, method='length')
							);
						}
					} else {
						if(with_pts) {
							tib <- .self$summarise(
								'index' = list(col=indexname, method='list', sep=sep, with.braces=with_braces),
								'size' = list(col=indexname, method='length'),
								'positions' = list(col=by, method='list:points', sep=sep, with.keys=with_keys, with.braces=with_braces)
							);
						} else {
							tib <- .self$summarise(
								'index' = list(col=indexname, method='set', sep=sep, with.braces=with_braces),
								'size' = list(col=indexname, method='length')
							);
						}
					}

					return(tib);
				} else {
					cols <- .self$data.cols;
					clustername <- .self$cluster.name;
					group_by <- .self$cluster.by;
					group_by <- .self$cluster.group.by;
					keep <- .self$cols.keep;
					sel <- .self$cluster.loc;
					print(sel);

					if(!.self$cols.keep.set) keep <- cols;
					keep <- unique(c(clustername, group_by, by, keep));

					return(.self$cluster.data[sel, ] %>% dplyr::select(keep));
				}
			}

			return(NULL);
		},
		setkeep = function(keep) {
			cols <- .self$data.cols;
			keep <- keep[which(keep %in% cols)];
			.self$cols.keep <- keep;
			.self$cols.keep.set <- TRUE;
		},
		groupby = function(only_clusters=TRUE) {
			if(!.self$is.clustered) return(tibble::as_tibble(list()));
			tib <- .self$cluster.data;
			clustername <- .self$cluster.name;
			if(only_clusters) {
				sel <- .self$cluster.loc;
				tib <- tib[sel, ];
			}
			return(tib %>% group_by_at(clustername));
		},
		build = function(...) {
			obj <- .self$data %>% buildclusters____(...);
			for(key in c(
				'is.lexical',
				'index.name',
				'cluster.name',
				'cluster.loc',
				'cluster.by',
				'cluster.group.by',
				'cluster.data'
			)) .self[[key]] <- obj[[key]];
			.self$is.clustered <- TRUE;
		},
		summarise = function(...) {
			if(!.self$is.clustered) return(tibble::as_tibble(list()));

			INPUTVARS <- list(...);
			summcols <- names(INPUTVARS);

			clustername <- .self$cluster.name;
			group_by <- .self$cluster.group.by;
			keep <- .self$cols.keep;
			cols <- .self$data.cols;
			if(!.self$cols.keep.set) keep <- c();
			keep <- unique(c(group_by, keep));
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

			n <- 0;
			for(col in c(keep,summcols)) {
				if(col == clustername) next;

				s <- instructions[[col]];
				f <- (function(s, col) {
					f <- NULL;
					vars <- s[['col']];
					m <- s[['method']];
					d <- length(vars);

					mode <- s[['mode']];
					if(!is.character(mode)) mode <- '';
					sep <- ';';
					if('sep' %in% names(s)) sep <- s[['sep']];
					with_braces <- FALSE;
					if('with.braces' %in% names(s)) with_braces <- s[['with.braces']];
					with_keys <- FALSE;
					if('with.keys' %in% names(s)) with_keys <- s[['with.keys']];

					lbrace <- '';
					rbrace <- '';
					if(with_braces) {
						lbrace <- '[';
						rbrace <- ']';
					}

					llbrace <- '';
					rrbrace <- '';
					if(d > 1) {
						llbrace <- '[';
						rrbrace <- ']';
					}

					if(with_keys) {
						inner_to_json <- function(x) {return(jsonlite::toJSON(x));};
					} else {
						inner_to_json <- function(x) {return(paste0(llbrace,paste(x, collapse=sep),rrbrace));};
					}

					if(is.function(m)) {
						f <- opt;
					} else if(is.character(m)) {
						opt <- m[1];
						if(opt == 'pick') {
							f <- function(x) {return(x[1]);};
						} else if(opt == 'set') {
							f <- function(x) {return(paste0(lbrace,paste(unique(x), collapse=sep),rbrace));};
						} else if(opt == 'list') {
							f <- function(x) {return(paste0(lbrace,paste(x, collapse=sep),rbrace));};
						} else if(opt == 'list:points') {
							if(with_keys) {
								f <- function(...) {
									pts <- list(...);
									names(pts) <- vars;
									return(jsonlite::toJSON(tibble::as_tibble(pts)));
								};
							} else {
								f <- function(...) {
									pts <- list(...);
									names(pts) <- c(1:length(pts));
									pts <- tibble::as_tibble(pts);
									n <- nrow(pts);
									return(paste0(lbrace, paste(sapply(c(1:n), function(i) {
											return(inner_to_json(pts[i, ]));
									}), collapse=sep), rbrace));
								};
							}
						} else if(opt == 'json') {
							f <- jsonlite::toJSON;
						} else if(opt == 'json:set') {
							f <- function(x) {return(jsonlite::toJSON(unique(x)));};
						} else if(opt == 'length') {
							f <- length;
						} else if(opt == 'min') {
							f <- function(x) {return(min(x, na.rm=TRUE));};
						} else if(opt == 'max') {
							f <- function(x) {return(max(x, na.rm=TRUE));};
						} else if(opt == 'range') {
							if(mode == 'distance') {
								f <- function(x) {
									a <- min(x, na.rm=TRUE);
									b <- max(x, na.rm=TRUE);
									return(b - a);
								};
							} else {
								f <- function(x) {
									a <- min(x, na.rm=TRUE);
									b <- max(x, na.rm=TRUE);
									return(paste0('[',a,';',b,']'));
								};
							}
						} else if(opt == 'lex:min') {
							f <- function(...) {
								pts <- list(...);
								nom <- sapply(c(1:length(pts)), function(i) {return(paste0('c',i));});
								names(pts) <- nom;
								pts <- lexsort____(tibble::as_tibble(pts), nom);
								names(pts) <- vars;
								pt <- pts[1, ];
								return(inner_to_json(pt));
							};
						} else if(opt == 'lex:max') {
							f <- function(...) {
								pts <- list(...);
								nom <- sapply(c(1:length(pts)), function(i) {return(paste0('c',i));});
								names(pts) <- nom;
								pts <- lexsort____(tibble::as_tibble(pts), nom);
								names(pts) <- vars;
								pt <- pts[nrow(pts), ];
								return(inner_to_json(pt));
							};
						} else if(opt == 'lex:range') {
							if(mode == 'distance') {
								f <- function(...) {
									pts <- list(...);
									nom <- sapply(c(1:length(pts)), function(i) {return(paste0('c',i));});
									names(pts) <- nom;
									pts <- lexsort____(tibble::as_tibble(pts), nom);
									names(pts) <- vars;
									a <- pts[1,];
									b <- pts[nrow(pts),];
									d <- b-a;
									return(inner_to_json(d));
								};
							} else {
								if(with_keys) {
									f <- function(...) {
										pts <- list(...);
										nom <- sapply(c(1:length(pts)), function(i) {return(paste0('c',i));});
										names(pts) <- nom;
										pts <- lexsort____(tibble::as_tibble(pts), nom);
										names(pts) <- vars;
										a <- pts[1,];
										a <- inner_to_json(a);
										b <- pts[nrow(pts),];
										b <- inner_to_json(b);
										return(paste0('[',a,',',b,']'));
									};
								} else {
									f <- function(...) {
										pts <- list(...);
										nom <- sapply(c(1:length(pts)), function(i) {return(paste0('c',i));});
										names(pts) <- nom;
										pts <- lexsort____(tibble::as_tibble(pts), nom);
										names(pts) <- vars;
										a <- pts[1,];
										b <- inner_to_json(a);
										b <- pts[nrow(pts),];
										b <- inner_to_json(b);
										return(paste0(lbrace,a,sep,b,rbrace));
									};
								}
							}
						} else if(opt == 'mean') {
							f <- function(x) {return(mean(x, na.rm=TRUE));};
						} else if(opt == 'var') {
							f <- function(x) {return(var(x, na.rm=TRUE));};
						} else if(opt == 'sd') {
							f <- function(x) {return(sd(x, na.rm=TRUE));};
						} else {
							# f <- function(x) {return(NA);};
						}
					}

					return(f);
				})(s, col);

				if(is.null(f)) next;

				method[[k]] <- f;
				vnames <- paste(s[['col']], collapse=',');
				fnames[k] <- paste0('method[[',k,']](',vnames,')');
				colnames[k] <- col;
				if(col == 'size') k0 <- k;
				k <- k + 1;
				n <- n + 1;
			}

			tib <- .self$groupby(TRUE);
			tib_summ <- list();
			if(n > 0) {
				for(k in c(1:n)) {
					col <- colnames[k];
					args <- setNames(fnames[k], col);
					tib_summ[[col]] <- (tib %>% dplyr::summarise_(.dots=args))[[col]];
				}
			}

			return(tibble::as_tibble(tib_summ));
		}
	)
);




buildclusters____ <- function(tib, ...) {
	INPUTVARS <- list(...);
	VARNAMES <- names(INPUTVARS);

	group_by <- INPUTVARS[['filter.by']];
	by <- INPUTVARS[['by']];
	clustername <- INPUTVARS[['cluster.name']];
	min_cluster_size <- INPUTVARS[['min.size']];
	max_cluster_size <- INPUTVARS[['max.size']];
	near <- INPUTVARS[['near']];
	d_min <- INPUTVARS[['min.dist']];
	d_max <- INPUTVARS[['max.dist']];
	strict <- INPUTVARS[['strict']];
	is_lexical <- INPUTVARS[['is.lexical']];
	no_overlaps <- INPUTVARS[['no.overlaps']];
	presort <- INPUTVARS[['presort']];

	tib <- tibble::as_tibble(tib);
	cols <- names(tib);

	dim_by <- length(by);
	if(!is.vector(group_by)) group_by <- c();
	if(!is.logical(is_lexical)) is_lexical <- (dim_by == 1);
	if(!is.logical(no_overlaps)) no_overlaps <- is_lexical;
	if(!is.logical(presort)) presort <- is_lexical;
	if(!is.numeric(min_cluster_size)) min_cluster_size <- 1;
	if(!is.numeric(max_cluster_size)) max_cluster_size <- Inf;
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


	## Erstellung von Spaltennamen (Klusterspalte + Pufferspalte):
	n <- nrow(tib);
	if(!is.character(clustername)) clustername <- uniquecolumnname____('cluster', cols);
	chunkname <- uniquecolumnname____('chunk', c(cols, clustername));
	indexname <- uniquecolumnname____('index', c(cols, clustername));


	## Definiere Typen von Spalten
	types <- list()
	types[[indexname]] <- 'integer';
	types[[clustername]] <- 'integer';


	## Indexnamen bevor Verarbeitung speichern.
	tib <- tib %>% tibble::add_column(!!(indexname):=c(1:n), !!(clustername):=rep(NA, n));
	## Indizes für Präsortierung hinzufügen.
	groupname <- uniquecolumnname____('group', c(cols, clustername));
	tib <- tib %>% group_by_at(group_by) %>% nest(.key=!!(chunkname));
	Ng <- nrow(tib);
	tib <- tib %>% tibble::add_column(!!(groupname):=c(1:Ng));
	tib <- tib %>% unnest();
	## ursp. Reihenfolge wiederherstellen.
	tib <- tib %>% dplyr::arrange_(indexname);


	## HAUPTMETHODE
	if(is_lexical) {
		if(presort) tib <- tib %>% lexsort____(by);
		if(!no_overlaps) tib <- tib %>% lexsort____(groupname);
		tib_ <- tib;

		i0 <- 1;
		cl <- 1;
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
			if(min_cluster_size <= sz && sz <= max_cluster_size) {
				dsel <- c(i0:i1);
				tib[dsel, clustername] <- cl;
				cl <- cl + 1;
			}
			i0 <- i1 + 1;
		}
	} else {
		tib <- tib %>% lexsort____(groupname);
		cl <- 1;
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
			clusters <- generateconnectedcomponents____(edges, min_cluster_size, max_cluster_size, cl);
			dsel <- which(!is.na(clusters));
			if(length(dsel) > 0) cl <- max(clusters[dsel]) + 1;
			tib[ind, clustername] <- clusters;
		}
	}


	## Reinige Typen.
	for(col in names(types)) {
		typ <- types[[col]];
		if(typ == 'integer') tib[ , col] <- as.integer(tib[[col]]);
	}

	## Reihenfolge wiederherstellen und Spalten filtrieren.
	tib <- tib %>% dplyr::arrange_(indexname) %>% dplyr::select(-c(groupname));
	sel <- which(!is.na(tib[[clustername]]));

	return(list(
		is.lexical=is_lexical,
		index.name=indexname,
		cluster.name=clustername,
		cluster.loc=sel,
		cluster.by=by,
		cluster.group.by=group_by,
		cluster.data=tib
	));
};


## LOKALE METHODEN
generateconnectedcomponents____ <- function(edges, min_sz, max_sz, key0) {
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
		nodes <- sort(nodes);
		while(length(nodes) >= max_sz) {
			classes[nodes[c(1:max_sz)]] <- key;
			key <- key + 1;
			nodes <- nodes[-c(1:max_sz)];
		}
		if(length(nodes) < min_sz) next;
		classes[nodes] <- key;
		key <- key + 1;
	}

	return(classes);
};

uniquecolumnname____ <- function(nom, cols) {
	col <- nom;
	while(col %in% cols) col <- paste0(col,'_');
	return(col);
};

lexsort____ <- function(df, nom) {
	## data %>% arrange_(c(#,#,...,#)); funktioniert nicht.
	expr <- paste0("df[with(df, order(",paste(nom, collapse=','),")), ]");
	return(eval(parse(text=expr)));
};