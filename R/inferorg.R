
#' Summarize for a given organism
#'
#' @keywords internal
#'
#' @return A vector of length 3 giving the proportion of elements
#' in \code{gene_vec} that match the symbols, ensembl IDs and entrez
#' IDs for the organism that \code{gene_vec} corresponds to.
summarize_organism <- function(df, gene_vec) {

  existing_symbols <- df$symbol[!is.na(df$symbol)]
  existing_ens <- df$ensgene[!is.na(df$ensgene)]
  existing_entrez <- df$entrez[!is.na(df$entrez)]

  prop_symbol <- mean(as.character(gene_vec) %in% existing_symbols, na.rm = TRUE)
  prop_ens <- mean(as.character(gene_vec) %in% existing_ens, na.rm = TRUE)
  prop_entrez <- mean(as.character(gene_vec) %in% existing_entrez, na.rm = TRUE)

  c(`symbol` = prop_symbol, `ensgene` = prop_ens, `entrez` = prop_entrez)
}

#' Infer the organism and ID format given a list of genes.
#'
#' @import annotables
#'
#' @export
#'
#' @details
#'
#' This function takes in a vector of gene IDs, and guesses the format, which is one of:
#' \itemize{
#' \item \code{symbol} (e.g. PTPRC)
#' \item \code{ensgene} (e.g. ENSG00000081237)
#' \item \code{entrez} (e.g. 5788)
#' }
#' along with the organism, which is one of:
#' \itemize{
#' \item \code{human}
#' \item \code{mouse}
#' \item \code{fruit_fly}
#' \item \code{macaque}
#' \item \code{worm}
#' \item \code{chicken}
#' \item \code{rat}
#' }
#' In addition, a confidence score for both the organism and format is computed
#' (see vignette for details).
#'
#' @param gene_vec Vector of input gene IDs to guess the organism and format of
#'
#' @return A list with four entries:
#' \itemize{
#' \item \code{organism} The inferred organism.
#' \item \code{format} The inferred ID format.
#' \item \code{confidence_organism} The confidence in the organism inferred, between 0 and 1.
#' \item \code{confidence_format} The confidence in the format inferred, between 0 and 1.
#' }
#' @examples
#' genes <- c("HLA-A", "HLA-B", "HLA-C")
#' inferorg(genes)
inferorg <- function(gene_vec) {

  if(all(is.na(gene_vec))) {
    stop("Input genes are all NA")
  }

  if(! (is.character(gene_vec) || is.numeric(gene_vec))  ) {
    stop("Input genes must be either character or integer vector")
  }

  gene_vec <- as.character(gene_vec)

  genomes <- get_genomes()

  tbl <- sapply(genomes, summarize_organism, gene_vec)

  max_val <- max(tbl)
  guess_index <- which(tbl == max_val, arr.ind=TRUE)[1,]

  guess_format <- rownames(tbl)[guess_index[1]]
  guess_organism <- colnames(tbl)[guess_index[2]]

  format_confidence <- max_val * tbl[guess_index[1], guess_index[2]] / sum(tbl[, guess_index[2]])
  organism_confidence <- max_val * tbl[guess_index[1], guess_index[2]] / sum(tbl[guess_index[1], ])

  if(is.nan(format_confidence)) {
    format_confidence <- 0
  }
  if(is.nan(organism_confidence)) {
    organism_confidence <- 0
  }

  list(
    'organism' = guess_organism,
    'format' = guess_format,
    'confidence_organism' = organism_confidence,
    'confidence_format' = format_confidence
  )
}

#' Automatically convert genes to a given format
#'
#' Infers both the organism and genome ID format before converting
#' to the desired format
#'
#' @param gene_vec Vector of input gene IDs
#' @param to ID format to convert gene IDs to
#'
#' @return A vector of the same size as gene_vec with the mapped IDs. Any
#' IDs that are unable to be mapped are converted to \code{NA}.
#'
#' @importFrom plyr mapvalues
#'
#' @export
#'
#' @examples
#' autoconvert(c("Tapbp", "B2m"), to='ensgene')
#' autoconvert(c(3105, 3106), to='symbol')
autoconvert <- function(gene_vec, to = c('symbol', 'ensgene', 'entrez')) {
  to <- match.arg(to)

  io <- inferorg(gene_vec)

  if(io$confidence_format == 0 || io$confidence_organism == 0) {
    return(rep(NA, length(gene_vec)))
  }

  if(io$format == to) {
    return(gene_vec)
  }

  df <- get_genomes()[[ io$organism ]]
  df <- df[df[[io$format]] %in% gene_vec,]

  new_ids <- suppressMessages({
    plyr::mapvalues(gene_vec,
                    from = df[[io$format]],
                    to = df[[to]])
  })
  new_ids[new_ids == gene_vec] <- NA
  new_ids
}

#' @keywords internal
get_genomes <- function() {
  list(
    human = annotables::grch38[,c('symbol', 'ensgene', 'entrez')],
    mouse = annotables::grcm38[,c('symbol', 'ensgene', 'entrez')],
    worm = annotables::wbcel235[,c('symbol', 'ensgene', 'entrez')],
    `fruit_fly` = annotables::bdgp6[,c('symbol', 'ensgene', 'entrez')],
    `macaque` = annotables::mmul801[,c('symbol', 'ensgene', 'entrez')],
    chicken = annotables::galgal5[,c('symbol', 'ensgene', 'entrez')],
    rat = annotables::rnor6[,c('symbol', 'ensgene', 'entrez')]
  )
}

