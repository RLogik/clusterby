#' Data-clustering
#'
#' This package contains methods, which enables clustering in dataframes. Particularly useful for bio-mathematics, cognitive sciences, etc.
#'
#' \code{cluster(df, ...)}
#' @param df Tibble/Dataframe to be clustered. Method also possible with vectors.
#' @param groupby Defaults to \code{c()}. Specificies columns, by which data is to be preliminarily divided into groups, within which the clusters are to be built.
#' @param by Specifies the column(s) for geometric data, according to which the clusters are to be built.
#' @param isnear A function. This function operates pairs of entries in the columns with geometric data and returns \code{TRUE}/\code{FALSE} if entries are near. Defaults to a Euclidian metric.
#' @param max.dist Defaults to \code{Inf}. If the default euclidean metric is used for \code{isnear}, this is the maximum tolerated distance between geometric data.
#' @param strict Defaults to \code{TRUE}. If the default euclidean metric is used for \code{isnear}, this sets the proximity to be a strict \code{< dist} or else \code{<= dist}.
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
	splqit <- INPUTVARS[['split']];
	isnear <- INPUTVARS[['near']];
	d_max <- INPUTVARS[['max.dist']];
	strict <- INPUTVARS[['strict']];

	if(!is.logical(strict)) strict <- FALSE;
	if(!is.logical(split)) split <- FALSE;
	if(!is.numeric(d_max)) d_max <- Inf;
	if(!is.numeric(min_cluster_size)) min_cluster_size <- 0;
	if(!is.function(isnear)) {
		if(strict) {
			isnear <- function(x, y) {
				dx <- sqrt(sum((x-y)^2));
				# dx <- max(abs(x-y));
				return(dx < d_max);
			};
		} else {
			isnear <- function(x, y) {
				dx <- sqrt(sum((x-y)^2));
				# dx <- max(abs(x-y));
				return(dx <= d_max);
			};
		}
	}
	if(!is.character(clustername)) clustername <- 'cluster';


	## Erstellung von Spaltennamen (Klusterspalte + Pufferspalte):
	cols <- names(data);
	data <- tibble::as_tibble(data);
	chunkname <- 'chunk'; i <- 0; while(chunkname %in% c(cols, clustername)) {chunkname <- paste0('chunk', i); i <- i+1;}
	edgesname <- 'edges'; i <- 0; while(edgesname %in% c(cols, edgesname)) {edgesname <- paste0('edges', i); i <- i+1;}

	## Prägruppierung der Daten:
	if(!is.vector(groupby)) groupby <- c();
	leer <- list();
	n <- nrow(data);
	for(i in c(1:n)) leer[[i]] <- NA;
	data[, clustername] <- leer;
	data[, edgesname] <- leer;
	# data <- add_column(data, !!(clustername) := leer);
	# data <- add_column(data, !!(edgesname) := leer);

	# Erzeuge Kanten für Klusterbausteine.
	tib <- data %>% group_by_at(groupby) %>% nest(.key=!!(chunkname));
	m <- nrow(tib);
	for(i in c(1:m)) {
		chunk <- tib[i, chunkname][[1]][[1]];
		n <- nrow(chunk);
		for(i in c(1:(n-1))) {
			r <- 0;
			e <- c();
			pt_1 <- chunk[i, by];
			for(j in c((i+1):n)) {
				pt_2 <- chunk[j, by];
				if(isnear(pt_1, pt_2)) {
					r <- r + 1;
					e[r] <- j;
				}
			}
			chunk[i, edgesname] <- e;
		}
		tib[i, chunkname][[1]][[1]] <- chunk;
	}

	# Erzeuge Kluster aus Kanten.
	tib <- tib %>% unnest();
	clusters <- generateclasses(tib[, edgesname], min_cluster_size);
	ind <- which(!is.na(clusters));
	clusters <- clusters[ind];
	tib <- chunk[ind, ];
	tib[, clustername] <- clusters;

	if(split) tib <- tib %>% group_by_at(c(groupby, clustername)); #%>% nest(.key=!!(dataname);

	return(tib);
};


generateclasses <- function(edges, min_sz) {
	n <- nrow(edges);
	key <- 0;
	ind <- c(1:n);
	classes <- rep(NA,n);
	while(length(ind) > 0) {
		i <- ind[1];
		ind <- ind[-1];

		nodes = c();
		children <- c(i);
		while(TRUE) {
			grandchildren <- c();
			for(j in children) {
				e <- edges[j,1];
				if(length(e) == 1) if(is.na(e)) next;
				grandchildren <- c(grandchildren, e);
			}
			children <- grandchildren;
			if(length(children) == 0) break;
			filt <- which(children %in% ind);
			children <- children[filt];
			nodes <- c(nodes, children);
			ind <- ind[-c(children)];
		}

		if(length(nodes) >= min_sz) {
			classes[nodes] <- key;
			key <- key + 1;
		}
	}

	return(classes);
};