#' Data-clustering
#'
#' This package contains methods, which enables clustering in dataframes. Particularly useful for bio-mathematics, cognitive sciences, etc.
#' 
#' \code{clusterby(df, ...)}
#' @param df Dataframe to be clustered. Method also possible with vectors.
#' @param by Specifies the column name for geometric data, according to which the clusters are to be built.
#' @param group Defaults to \code{c()}. Specificies columns, by which data is to be preliminarily divided into groups, within which the clusters are to be built.
#' @param dist Defaults to \code{Inf}. For a column \code{x} of positions, a cluster is a subset \code{y} which satisfies \bold{(i)} between every two points there is a path within \code{y} whose steps are 'small'; \bold{(ii)} y is maximal under condition (i). 'Small' here means the the distance between two points is strictly smaller than \code{d}.
#' @param strict Defaults to \code{TRUE}. If set to \code{FALSE} the 'smaller than' is replaced by 'smaller than or equal to'.
#' @param clustername Defaults to 'cluster'. Running \code{df \%>\% clusterby(...)} returns a data frame, which extends \code{df} by 1 column with this name. This column tags the clusters by a unique index.
#' @param dim New feature to be added. Allow points to be arbitrary data, allow input of a distance function/neighbourhoodscheme.
#' @keywords cluster clustering gene
#' @export clusterby
#' @export split
#' @export summarise
#' @examples df %>% clusterby(by='position', group=c('gene','actve'), dist=400, strict=FALSE, clustername='tag');
#' @examples p %>% clusterby(dist=10); # p is a numeric vector




################################################################
#### HAUPTMETHODE ##############################################
clusterby <- function(data, ...) {
	INPUTVARS <- list(...);
	VARNAMES <- names(INPUTVARS);
	d <- INPUTVARS[['dist']];
	tagcols <- INPUTVARS[['group']];
	by <- INPUTVARS[['by']];
	clustername <- INPUTVARS[['clustername']];
	mode <- INPUTVARS[['mode']];
	strict <- INPUTVARS[['strict']];
	summary <- INPUTVARS[['summary']];
	format <- INPUTVARS[['format']];

	strict <- ('strict' %in% VARNAMES) && (strict==TRUE);
	if(!('format' %in% VARNAMES)) {
		if('summary' %in% VARNAMES) {
			format <- 'list';
		} else {
			format <- 'original';
		}
	}
	if(!('dist' %in% VARNAMES) || !is.numeric(d)) d <- Inf;

	if(is.vector(data)) {
		## Geometrische Werkzeuge:
		mengendist <- function(X,Y) {D<-c(); for(x in X) D[length(D)+1] <- min(abs(Y-x)); return(D);};
		if(strict) {
			ball <- function(X,A,r) {return(which(mengendist(X,A) < r));};
		} else {
			ball <- function(X,A,r) {return(which(mengendist(X,A) <= r));};
		}

		## Iterativ die Kluster bilden:
		n <- length(data);
		index <- c(1:n);
		tags <- rep(NA,n);
		tag <- 0;
		firsttime <- TRUE;

		while(n > 0) {
			if(firsttime) {
				## Initialisiere Klustur.
				## Beschränke Daten auf übrig gebliebene Indizes.
				D <- data[index];
				## Ordner dieser Menge Indizes zu.
				I <- c(1:n);
				## Bilde einen Kluster bzgl. Indizes.
				J <- c(1);
				component <- index[J];
				firsttime <- FALSE;
			} else {
				Jneu <- ball(D,D[J],d);
				dJ <- setdiff(Jneu,J);
				dI <- index[dJ];
				component <- c(component,dI);

				## Entferne alte Indizes.
				index <- index[-c(J)];
				## Beschränke Daten auf übrig gebliebene Indizes.
				D <- data[index];
				## Indizes entsprechen dem Wachstum des Klusturs:
				J <- which(index %in% dI);
				n <- length(D);

				if(length(dI) == 0) {
					# Markiere die Indizes der Datenpunkte, die zum Kluster gehören
					tags[component] <- tag;
					# Starte mit einem neuen Kluster
					tag <- tag + 1;
					firsttime <- TRUE;
				}
			}
		}

		out <- tags;
		if(format == 'list') {
			
		}

		return(tags);
	} else if(is.data.frame(data)) {
		## Erstellung von Spaltennamen (Pufferspalte + Klusterspalte):
		cols <- names(data);
		if(!('clustername' %in% VARNAMES)) clustername <- 'cluster';
		tagname <- 'tag';
		i <- 0; while(tagname %in% c(cols,clustername)) {tagname <- paste('tag',i,collapse='');i <- i+1;}

		## Prägruppierung der Daten:
		if(!('group' %in% VARNAMES)) tagcols <- c();
		data <- data %>% tagby(by=tagcols, tagname=tagname);

		## Schleife durch alle Gruppierungen und bilde Kluster:
		tags <- unique(data[, tagname]);
		startingname <- 1;		
		for(tag in tags) {
			## Gruppe
			ind <- which(data[, tagname] == tag);
			## Bilde Kluster mit Vektormethode:
			clusters <- data[ind, by] %>% clusterby(dist=d, mode=mode, strict=strict);
			## Schreibe Klustertags in ursprünglichen Dataframe:
			for(i in c(1:length(ind))) data[ind[i], clustername] <- startingname + clusters[i];
			## Verhindere Überschneidungen in den Klustertagnamen:
			startingname <- startingname + max(0,clusters) + 1;
		}

		# for(cl in clusters) out[cl, ] <- df[[cl]][1, ];
		return(data[, !(names(data) %in% c(tagname))]);
	}

	return(NULL);
};


## Verwandelt geklusterte Daten in eine Liste.
split <- function(data, ...) {
	INPUTVARS <- list(...);
	VARNAMES <- names(INPUTVARS);
	clustername <- INPUTVARS[['clustername']];
	if(!('clustername' %in% VARNAMES)) clustername <- 'cluster';

	clusters <- unique(as.vector(data[, clustername]));
	out <- list();
	i <- 0;
	for(cl in clusters) {
		i <- i+1;
		df <- data[which(data[, clustername] == cl), !(names(data) %in% c(clustername))];
		rownames(df) <- NULL;
		out[[i]] <- df;
	}

	return(out);
};


## Fasst eine Liste von Klustern zusammen.
summarise <- function(llist, ...) {
	INPUTVARS <- list(...);
	VARNAMES <- names(INPUTVARS);

	if(length(llist) == 0) return(c());

	cols <- names(llist[[1]]);
	out <- list();
	blank <- rep(NA,length(llist));
	for(col in cols) out[[col]] <- blank;
	out <- as.data.frame(out);
	rownames(out) <- NULL;

	for(col in cols) {
		f <- function(x) {return(NA);};
		if(col %in% VARNAMES) {
			s <- INPUTVARS[[col]];
			if(mode(s) == 'function') {
				f <- s;
			} else if(mode(s) == 'character') {
				if(s[1] == 'set') {
					sep = ';'; if(length(s) > 1) sep = s[2];
					f <- function(x) {return(paste(unique(x),collapse=sep));};
				} else if(s[1] == 'list') {
					sep = ';'; if(length(s) > 1) sep = s[2];
					f <- function(x) {return(paste(x,collapse=sep));};
				} else if(s[1] == 'length') {
					f <- length;
				} else if(s[1] == 'min') {
					f <- min;
				} else if(s[1] == 'max') {
					f <- max;
				} else if(s[1] == 'range') {
					f <- function(x) {return(paste(min(x),max(x),sep='-'));};
				} else if(s[1] == 'mean') {
					f <- mean;
				} else if(s[1] == 'var') {
					f <- var;
				} else if(s[1] == 'sd') {
					f <- sd;
				}
			}
		} else {
			f <- function(x) {return(x[1]);};
		}

		i <- 0;
		for(data in llist) {
			i <- i + 1;
			out[i, col] <- f(as.vector(data[, col]));
		}
	}

	return(out);
};




################################################################
#### LOKALE FUNKTION: tagby(data, ...) #########################
####    markiert Zeilen in einem Dataframe zwecks Gruppierung ##
tagby <- function(data, ...) {
	INPUTVARS <- list(...);
	VARNAMES <- names(INPUTVARS);
	by <- INPUTVARS[['by']];
	tagname <- INPUTVARS[['tagname']];

	data <- data %>% dplyr::arrange_(.dots=by);
	rownames(data) <- NULL;

	tag <- 0;
	firsttime <- TRUE;
	lastrow <- NULL;
	for(r in rownames(data)) {
		if(!firsttime) {
			change <- FALSE;
			for(col in by) change <- change || !(data[r,col] == data[rlast,col]);
			if(change) tag = tag + 1;
			rlast <- r;
		}
		data[r,tagname] <- tag;
		rlast <- r;
		firsttime <- FALSE;
	}

	return(data);
};