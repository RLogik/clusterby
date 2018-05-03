#' Data-clustering
#'
#' This package contains the method \code{clusterby}, which enables clustering in dataframes. Particularly useful for biomathematics, cognitive sciences, etc.
#' 
#' \code{clusterby(df, ...)}
#' @param df Dataframe to be clustered. Method also possible with vectors.
#' @param by Specifies the column name for geometric data, according to which the clusters are to be built.
#' @param group Defaults to \code{c()}. Specificies columns, by which data is to be preliminarily divided into groups, within which the clusters are to be built.
#' @param dist Defaults to \code{Inf}. For a column \code{x} of positions, a cluster is a subset \code{y} which satisfies \bold{(i)} between every two points there is a path within \code{y} whose steps are 'small'; \bold{(ii)} y is maximal under condition (i). 'Small' here means the the distance between two points is strictly smaller than \code{d}.
#' @param strict Defaults to \code{TRUE}. If set to \code{FALSE} the 'smaller than' is replaced by 'smaller than or equal to'.
#' @param clustername Defaults to 'cluster'. Running \code{df \%>\% clusterby(...)} returns a data frame, which extends \code{df} by 1 column with this name. This column tags the clusters by a unique index.
#' @keywords cluster clustering gene
#' @export clusterby
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
	strict <- INPUTVARS[['strict']];

	strict <- ('strict' %in% VARNAMES) && (strict==TRUE);

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

		return(tags);
	} else if(is.data.frame(data)) {
		## Erstellung von Spaltennamen (Pufferspalte + Klusterspalte):
		cols <- names(data);
		if(!('clustername' %in% VARNAMES)) clustername <- 'cluster';
		tagname <- 'tag';
		i <- 0; while(tagname %in% c(cols,clustername)) {tagname <- paste('tag',i,sep='');i <- i+1;}

		## Prägruppierung der Daten:
		if(!('group' %in% VARNAMES)) tagcols <- c();
		data <- data %>% tagby(by=tagcols, tagname=tagname);

		## Schleife durch alle Gruppierungen und bilde Kluster:
		tags <- unique(data[, tagname]);
		startingname <- 0;		
		for(tag in tags) {
			## Gruppe
			ind <- which(data[, tagname] == tag);
			## Bilde Kluster mit Vektormethode:
			clusters <- data[ind, by] %>% clusterby(dist=d, strict=strict);
			## Schreibe Klustertags in ursprünglichen Dataframe:
			for(i in c(1:length(ind))) data[ind[i], clustername] <- startingname + clusters[i];
			## Verhindere Überschneidungen in den Klustertagnamen:
			startingname <- startingname + max(0,clusters) + 1;
		}

		return(data[, !(names(data) %in% c(tagname))]);
	}
	return(NULL);
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
