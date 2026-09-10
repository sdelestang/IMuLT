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

#' Get Selectivity Parameter Names from a Fitted Model's SELEXSPEC.DAT
#'
#' Derives readable parameter names for SelPars_n directly from the comment
#' text in SELEXSPEC.DAT (as loaded into `selx`), rather than reconstructing
#' block order from the Selectivity workbook sheet. This guarantees names
#' line up with the actual TMB parameter order even if the workbook has
#' since been edited/reordered relative to when the model was fitted.
#'
#' @param selx Data frame as read from SELEXSPEC.DAT (fixed-width columns,
#'   e.g. via read.table(..., col.names=1:200))
#' @return Character vector of names, one per SelPars_n, in file order,
#'   prefixed with 'Sel_'
GetSelectParNames <- function(selx){

  sel_comments <- findNclean(c('#','Selectivity','Parameters'), selx, 2, char = TRUE)
  comment_cols <- names(sel_comments)[10:ncol(sel_comments)]

  raw_comment <- apply(sel_comments[, comment_cols], 1, function(r) {
    r <- r[!is.na(r) & trimws(r) != ""]
    r <- r[-1]                            # drop "logistic"
    n <- length(r)
    r <- r[c(1:3, 5:(n-3))]               # keep p1/p2, M, number, description; drop dup p1/p2 and trailing "id M n"
    paste(r, collapse = " ")
  })

  paste0('Sel_', trimws(raw_comment))
}

#### Determine data-driven ceiling for CPUE sigma (SigmaCpueCeiling) ####
## SigmaCpue (and SigmaCpueUse, which this ceiling caps) is a
## dimensionless multiplier on the input CV -- see IMuLT.cpp's
## CpueLikelihood(): the residual it computes is log(obs/pred)/CV_input,
## and SigmaCpue = sqrt(mean(residual^2)) is the RMS of that ratio, so
## SigmaCpue ~= 1 means the input CVs already explain the year-to-year
## noise. A ceiling picked once by hand risks clipping genuine index
## noise in a future dataset where that noise happens to be larger, or
## being needlessly loose where it's smaller -- so this derives it fresh
## from Udat every time Filebuilder runs.
##
## For each CpueInd (the same grouping "Treatment of sigma" uses), fits
## log(Index) to a low-order polynomial trend in Year with no population
## model involved at all, and computes the RMS of residual/CV -- the same
## statistic SigmaCpue computes, but against a smooth trend instead of a
## fitted model. That's a model-independent floor on plausible extra
## variance: even a population model fitting a series perfectly on trend
## would still show at least this much SigmaCpue, purely from index
## noise. Series with too few points (dof < min_dof) are excluded from
## driving the ceiling -- with only a couple of residual degrees of
## freedom the estimate is too noisy to be informative either way.
EstimateCpueSigmaCeiling <- function(Udat, min_dof = 5, safety_factor = 1.2,
                                     min_ceiling = 1.5, verbose = TRUE) {

  selfRMS <- Udat %>%
    group_by(CpueInd) %>%
    group_modify(~{
      d <- .x %>% arrange(Year)
      n <- nrow(d)
      logidx <- log(d$Index)
      deg <- min(3, max(1, floor((n - 2) / 3)))
      deg <- min(deg, n - 2)
      fit <- lm(logidx ~ poly(d$Year, deg, raw = TRUE))
      resid_std <- residuals(fit) / d$CV
      dof <- n - (deg + 1)
      data.frame(n = n, deg = deg, dof = dof,
                 selfRMS = sqrt(mean(resid_std^2)))
    }) %>%
    ungroup()

  #if (verbose) print(selfRMS)

  reliable <- selfRMS %>% filter(dof >= min_dof)
  if (nrow(reliable) == 0) {
    warning("No CPUE series has enough points (dof >= ", min_dof, ") to ",
            "estimate a data-driven SigmaCpueCeiling -- falling back to ",
            min_ceiling, ". Lower min_dof, or set the ceiling by hand for ",
            "this run.")
    return(min_ceiling)
  }

  ceiling <- max(max(reliable$selfRMS) * safety_factor, min_ceiling)
  # if (verbose) {
  #   cat("Data-driven SigmaCpueCeiling:", round(ceiling, 3),
  #       "(max selfRMS", round(max(reliable$selfRMS), 3),
  #       "x safety factor", safety_factor, ")\n")
  # }
  ceiling
}

## Model-independent self-consistency check, same method used to derive
## SigmaCpueCeiling in Filebuilder: fit each series' log(Observed) to a
## low-order polynomial trend in Year (no population model involved), and
## compute the RMS of residual/Relative_CV -- the same statistic Sigma
## computes, but against a smooth trend instead of the fitted model. The
## gap between Sigma_used and this tells you how much of a series' capped
## resistance is genuine index noise vs. still-unresolved model tension.
EstimateCpueSelfRMS <- function(cpue_raw, min_dof = 5) {
  cpue_raw %>%
    group_by(Data_set) %>%
    group_modify(~{
      d <- .x %>% arrange(Year)
      n <- nrow(d)
      if (n < 3) return(data.frame(n = n, dof = NA_real_, selfRMS = NA_real_))
      deg <- min(3, max(1, floor((n - 2) / 3)))
      deg <- min(deg, n - 2)
      fit <- lm(log(d$Observed) ~ poly(d$Year, deg, raw = TRUE))
      dof <- n - (deg + 1)
      data.frame(n = n, dof = dof,
                 selfRMS = sqrt(mean((residuals(fit) / d$Relative_CV)^2)))
    }) %>% ungroup()
}
