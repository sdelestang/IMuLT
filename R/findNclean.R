#' Find and Extract a Data Block from a Control or Data File
#'
#' Searches a data file (read as a data.frame) for a row matching one or more
#' keyword patterns, then extracts and cleans the block of data between that
#' row and the next section header (denoted by \code{#}). Comment rows, empty
#' rows, and trailing empty columns are stripped. By default all columns are
#' converted to numeric; columns to protect from conversion can be specified
#' via \code{convert}, or the entire block can be returned as character with
#' \code{char = TRUE}.
#'
#' @param KeyWord Character vector of length 1--4. Each element is matched
#'   (via \code{\link[base]{grepl}}) against the corresponding column of
#'   \code{DataFile}. All patterns must match in the same row for that row
#'   to be identified as the target header.
#' @param DataFile A data.frame, typically produced by
#'   \code{read.csv(..., header = FALSE)} or equivalent, representing a
#'   CTL / DAT input file with \code{#}-delimited section headers.
#' @param Offset Integer (default \code{1}). Number of rows to skip after
#'   the keyword row before the data block begins. Useful when the keyword
#'   row is followed by a column-name row that should also be skipped.
#' @param convert Integer vector (default \code{0}). Column indices to
#'   \emph{exclude} from numeric conversion (i.e.\ keep as character).
#'   The default of \code{0} means all columns are converted to numeric.
#' @param char Logical (default \code{FALSE}). If \code{TRUE}, return the
#'   raw character data.frame with no numeric conversion and no column names
#'   applied.
#'
#' @return One of:
#' \itemize{
#'   \item A \code{data.frame} with column names taken from the header row
#'         and numeric conversion applied (the typical case).
#'   \item A \code{data.frame} of character values when \code{char = TRUE}.
#'   \item A numeric vector when the extracted block has a single column or
#'         a single row that collapses to a vector.
#'   \item \code{NA} if the keyword is not found or the block is empty.
#' }
#'
#' @details
#' The function locates section boundaries by finding rows whose first column
#' contains \code{#}. The data block is everything between the keyword row
#' (plus \code{Offset}) and the next \code{#} boundary. Within that block,
#' comment rows (containing \code{#}) and blank rows are removed, trailing
#' \code{NA}/empty columns are trimmed, and \code{NaN} strings are replaced
#' with \code{0}. If the number of extracted column-name tokens does not match
#' the number of data columns, synthetic names (\code{a1, a2, ...}) are
#' appended to fill the gap.
#'
#' @examples
#' \dontrun{
#' # Read a CTL file as raw character columns
#' ctl <- read.csv("model.ctl", header = FALSE, fill = TRUE,
#'                 stringsAsFactors = FALSE)
#'
#' # Extract a single-keyword block
#' growth <- findNclean("Growth", ctl)
#'
#' # Two-keyword match (must match columns 1 and 2)
#' sel <- findNclean(c("Selectivity", "Fleet1"), ctl)
#'
#' # Keep column 1 as character (e.g. area labels)
#' catch <- findNclean("Catch", ctl, convert = 1)
#'
#' # Return raw character block for further parsing
#' raw <- findNclean("Maturity", ctl, char = TRUE)
#' }
#'
#' @name findNclean
#' @export
findNclean <- function(KeyWord, DataFile, Offset = 1, convert = 0, char = FALSE) {


  ## --- locate all section boundaries (rows starting with #) ---
  hash <- c(which(grepl("#", DataFile[, 1])), nrow(DataFile))

  ## --- find the keyword row ---
  nk <- length(KeyWord)
  if (nk == 1L) {
    pos1 <- which(grepl(KeyWord, DataFile[, 1]))
  } else if (nk == 2L) {
    pos1 <- which(2 == (grepl(KeyWord[1], DataFile[, 1]) +
                          grepl(KeyWord[2], DataFile[, 2])))
  } else if (nk == 3L) {
    pos1 <- which(3 == (grepl(KeyWord[1], DataFile[, 1]) +
                          grepl(KeyWord[2], DataFile[, 2]) +
                          grepl(KeyWord[3], DataFile[, 3])))
  } else if (nk == 4L) {
    pos1 <- which(4 == (grepl(KeyWord[1], DataFile[, 1]) +
                          grepl(KeyWord[2], DataFile[, 2]) +
                          grepl(KeyWord[3], DataFile[, 3]) +
                          grepl(KeyWord[4], DataFile[, 4])))
  }

  if (length(pos1) == 0L) return(NA)

  ## --- determine the end of the data block ---
  if (!(pos1 + 1) %in% hash) {
    pos2 <- hash[hash > (pos1 + Offset)][1] - 1
    adj  <- 0
  } else {
    pos2 <- hash[hash > (pos1 + 1 + Offset)][1] - 1
    adj  <- 1
  }

  if (is.na(pos2) || pos2 <= pos1) return(NA)

  ## --- extract and clean the block ---
  tmp <- DataFile[(pos1 + Offset):pos2, ]
  tmp <- tmp[!grepl("#", tmp[, 1], fixed = TRUE), ]
  tmp <- tmp[tmp[, 1] != "", ]

  if (nrow(tmp) == 0L) return(NA)

  ## column names from the header row
  rname <- DataFile[(pos1 + adj), ]
  rname <- gsub("#", "", rname)
  rname <- rname[!is.na(rname) & rname != "" & rname != "NA"]

  ## trim trailing empty columns
  maxcol <- max(which(!is.na(tmp) & tmp != "", arr.ind = TRUE)[, 2])
  tmp <- tmp[!is.na(tmp[, 1]), 1:maxcol]
  tmp[tmp == "NaN"] <- 0

  ## pad column names if they don't match data width
  if (!is.null(dim(tmp))) {
    if (length(rname) != ncol(tmp)) {
      rname <- c(rname, paste0("a", seq_len(200)))[1:ncol(tmp)]
    }
  }

  ## single vector case (no dim)
  if (is.null(ncol(tmp))) return(as.numeric(tmp))

  ## --- numeric conversion ---
  convert_cols <- setdiff(seq_len(ncol(tmp)), convert)
  tmp     <- data.frame(tmp)
  chartmp <- tmp

  if (nrow(tmp) == 1L) {
    suppressWarnings(
      tmp[convert_cols] <- lapply(tmp[, convert_cols, drop = FALSE],
                                  function(q) as.numeric(as.character(q)))
    )
  } else {
    suppressWarnings(
      tmp[, convert_cols] <- data.frame(
        apply(as.matrix(tmp[, convert_cols]), 2,
              function(q) as.numeric(as.character(q)))
      )
    )
  }

  tmp <- tmp[!is.na(tmp[, 1]), !is.na(tmp[1, ])]

  if (isTRUE(char)) return(chartmp)

  if (!is.null(dim(tmp))) colnames(tmp) <- rname[seq_len(ncol(tmp))]
  return(tmp)
}
