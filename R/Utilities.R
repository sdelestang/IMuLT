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

  fleetyr <- egap %>% dplyr::select(starts_with('fleet'))
  if(ncol(fleetyr)==1){
    tfleetyr <- fleetyr
    colnames(tfleetyr) <- 'fleet999'
    fleetyr <- cbind(fleetyr, tfleetyr)
  }
  fleetyr <- fleetyr[egap$season %in% startseason:endseason,]
  pars <- sort(unique(as.vector(as.matrix(fleetyr))))

  ## Per-parameter row table, exactly as written to SELEXSPEC.dat
  ## (keep raw `id` here -- needed for the order safety check below)
  tegappar <- egappar %>%
    mutate(order=1:nrow(egappar), hash='#', id2=paste(uniq, id, comment)) %>%
    dplyr::select(order, lwr, upr, par, phase, Link, useprior, mnprior,
                  sdprior, hash, form, id, id2, yearlink, uniq) %>%
    arrange(order) %>%
    mutate(phase=ifelse(Link<=0, phase, -abs(phase))) %>%
    dplyr::select(-order)

  ## Duplicate blocks for any decimal yearlink not already present
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

  ## --- Safety check: enforce canonical parameter order within each block ---
  ## Downstream code matches parameters by string content (id2), not row
  ## position, but the *written* order in SELEXSPEC.dat (and therefore the
  ## SelPars_n indexing TMB uses) depends on whatever row order survives
  ## here. Force a fixed, known-good order per selectivity form so a
  ## differently-ordered Excel sheet can't silently scramble the .dat file.
  canonical_order <- list(
    logistic       = c('p1','p2'),
    doublelogistic = c('p1','p2','p3','p4')
  )

  block_ids <- unique(tegappar$uniq)
  ordered_list <- vector("list", length(block_ids))
  for(b in seq_along(block_ids)){
    blk <- tegappar[tegappar$uniq == block_ids[b], ]
    form_type <- unique(blk$form)

    if(length(form_type) != 1){
      warning(sprintf("ExpandSelectPars: block '%s' has mixed/ambiguous form values (%s) -- order not checked",
                      block_ids[b], paste(form_type, collapse=', ')))
      ordered_list[[b]] <- blk
      next
    }

    expected <- canonical_order[[form_type]]
    if(is.null(expected)){
      warning(sprintf("ExpandSelectPars: block '%s' has unrecognised form '%s' -- order not checked",
                      block_ids[b], form_type))
      ordered_list[[b]] <- blk
      next
    }

    if(!setequal(blk$id, expected) || length(blk$id) != length(expected)){
      warning(sprintf("ExpandSelectPars: block '%s' (form=%s) parameter ids don't match expected set.\n  Found:    %s\n  Expected: %s",
                      block_ids[b], form_type,
                      paste(blk$id, collapse=', '), paste(expected, collapse=', ')))
      ordered_list[[b]] <- blk
      next
    }

    ordered_list[[b]] <- blk[match(expected, blk$id), ]
  }
  tegappar <- do.call(rbind, ordered_list)

  tegappar
}

#' Get descriptive names for Selectivity link patterns
#'
#' Internal helper mirroring the pattern-numbering logic in
#' BuildInputFiles (egappar_sum) so MakeOutPut's selectivity plot
#' legends can label each pattern with its Excel comment directly,
#' rather than re-parsing whitespace-flattened text out of
#' SELEXSPEC.dat (which breaks whenever comment word-counts differ
#' between blocks, e.g. "Area 22" vs "Cameras").
#'
#' @param wb workbook object
#' @param startseason,endseason numeric
#' @return data.frame(link, name) -- one row per selectivity pattern.
#'   `link` is 0-indexed, matching the link/pointer values written to
#'   SELEXSPEC.dat and read back via fleet2.
#' @keywords internal
GetSelectPatternNames <- function(wb, startseason, endseason){

  egap <- readWorkbook(wb, sheet='Selectivity', startRow = 2)
  egappar <- egap %>%
    filter(!is.na(yearlink)) %>%
    mutate(uniq=paste(sex,yearlink)) %>%
    dplyr::select(!starts_with('fleet')) %>%
    mutate(Sex=ifelse(sex=='F',0,1), Sex=Sex-min(Sex))

  fleetyr <- egap %>% dplyr::select(starts_with('fleet'))
  if(ncol(fleetyr)==1){
    tfleetyr <- fleetyr
    colnames(tfleetyr) <- 'fleet999'
    fleetyr <- cbind(fleetyr, tfleetyr)
  }
  fleetyr <- fleetyr[egap$season %in% startseason:endseason,]
  pars <- sort(unique(as.vector(as.matrix(fleetyr))))

  ## One row per pattern (yearlink/sex/form block) -- same grouping and
  ## ordering as BuildInputFiles' egappar_sum, plus the descriptive comment
  egappar_sum <- egappar %>%
    group_by(yearlink, Sex, form, uniq) %>%
    summarise(num=length(Sex), comment=dplyr::first(comment), .groups='drop') %>%
    as.data.frame() %>%
    arrange(Sex) %>%
    mutate(pattern=as.numeric(rownames(.))-1) %>%
    dplyr::select(pattern, Sex, uniq, yearlink, comment)

  ## Duplicate blocks for decimal yearlinks (time-varying selectivity),
  ## carrying the source block's comment through to the duplicate
  egappar_sumog <- egappar_sum
  for(p in pars){
    if(!p %in% egappar_sum$yearlink){
      reps <- which(egappar_sumog$yearlink == floor(p))
      tmpe <- egappar_sum[reps,]
      tmpe$uniq <- paste(substr(tmpe$uniq,1,1), p)
      tmpe$pattern <- max(egappar_sum$pattern) + (1:nrow(tmpe))
      tmpe$yearlink <- p
      egappar_sum <- rbind(egappar_sum, tmpe)
    }
  }

  data.frame(link = egappar_sum$pattern, name = egappar_sum$comment)
}
