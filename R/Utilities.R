#' @title MatchTable
#'
#' @description This function finds the lines in a table that matches strings
#'
#' @param Table Name of the table
#' @param Char1 First character string to matrix
#' @param Char2 Second character string to matrix
#' @param Char3 Third character string to matrix
#' @param Char4 Fourth character string to matrix
#' @param Char5 Fifth character string to matrix
#'
#' @return vector of matching line indices
#' @export
#'
#' @examples
#' \dontrun{
#' }
#'
MatchTable<-function(Table,Char1=NULL,Char2=NULL,Char3=NULL,Char4=NULL,Char5=NULL)
{
  ii <- rep(T,length(Table[,1]))
  if (!is.null(Char1)) ii <- ii & (Table[,1]==Char1)
  if (!is.null(Char2)) ii <- ii & (Table[,2]==Char2)
  if (!is.null(Char3)) ii <- ii & (Table[,3]==Char3)
  if (!is.null(Char4)) ii <- ii & (Table[,4]==Char4)
  if (!is.null(Char5)) ii <- ii & (Table[,5]==Char5)
  ii <- seq(1:length(Table[,1]))[ii]
  return(ii)
}


#' Expand Selectivity parameter rows for time-varying yearlink duplicates
#'
#' Internal helper shared by BuildInputFiles and MakeOutPut so the
#' decimal-yearlink duplication logic for time-varying selectivity blocks
#' (e.g. yearlink 1 splitting into 1, 1.1, 1.2... when a fleet's selectivity
#' pattern changes partway through the time series) can't drift between the
#' two. Returns one row per parameter, in the same order parameters are
#' written to SELEXSPEC.dat, so the result can be indexed directly against
#' SelPars_<n> in the model output.
#'
#' @param wb An openxlsx workbook object (ModelStructure.xlsx)
#' @param startseason Numeric, first season of the model
#' @param endseason Numeric, last season of the model
#' @return A data.frame with one row per Selectivity parameter (including
#'   duplicated decimal-yearlink blocks), containing lwr, upr, par, phase,
#'   Link, useprior, mnprior, sdprior, hash, form, id2, yearlink, uniq --
#'   in write order. `id2` is the composite label (uniq + id + comment)
#'   identical to the ID field written to SELEXSPEC.dat.
#' @keywords internal
ExpandSelectPars <- function(wb, startseason, endseason){

  egap <- readWorkbook(wb, sheet='Selectivity', startRow = 2)
  egappar <- egap %>%
    filter(!is.na(yearlink)) %>%
    mutate(uniq=paste(sex,yearlink)) %>%
    dplyr::select(!starts_with('fleet')) %>%
    mutate(Sex=ifelse(sex=='F',0,1), Sex=Sex-min(Sex))

  ## Decimal yearlink values used anywhere in the fleet/season grid --
  ## these are the time-varying duplicate blocks (e.g. 1.1, 1.2...)
  fleetyr <- egap %>% dplyr::select(starts_with('fleet'))
  if(ncol(fleetyr)==1){
    tfleetyr <- fleetyr
    colnames(tfleetyr) <- 'fleet999'
    fleetyr <- cbind(fleetyr, tfleetyr)
  }
  fleetyr <- fleetyr[egap$season %in% startseason:endseason,]
  pars <- sort(unique(as.vector(as.matrix(fleetyr))))

  ## Per-parameter row table, exactly as written to SELEXSPEC.dat
  tegappar <- egappar %>%
    mutate(order=1:nrow(egappar), hash='#', id2=paste(uniq, id, comment)) %>%
    dplyr::select(order, lwr, upr, par, phase, Link, useprior, mnprior,
                  sdprior, hash, form, id2, yearlink, uniq) %>%
    arrange(order) %>%
    mutate(phase=ifelse(Link<=0, phase, -abs(phase))) %>%
    dplyr::select(-order)

  ## Duplicate blocks for any decimal yearlink not already present
  ## (time-varying selectivity within a yearlink group)
  tegapparog <- tegappar
  for(p in pars){
    if(!p %in% tegapparog$yearlink){
      reps <- which(tegapparog$yearlink == floor(p))
      tmpe <- tegapparog[reps,]
      tmpe$uniq <- paste(substr(tmpe$uniq,1,1), p)
      tmpe$yearlink <- p
      tegappar <- rbind(tegappar, tmpe)
    }
  }

  tegappar
}
