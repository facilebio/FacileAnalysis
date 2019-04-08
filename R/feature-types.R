#' Guesses the type of identifiers provided.
#'
#' This was extracted from the `guess_id_type` in the DenaliSigDb package.
#' We are using it as temporary bandaid to extract the "feature space" from
#' results. These should be straightforward to extract when all analyses are
#' going through some type of FacileDataStore.
#'
#' @description
#' A two-column data.frame is returned for id_type and organism. Organism
#' is "unknown" for identifiers where there this can't be inferred (like Refseq).
#'
#' If an identifier matches more than one id_type, the id_type is set to
#' `"ambiguous"`. If the identifier doesn't match any guesses, then `"unknown"`.
#'
#' @param x a character vector of ids
#' @return data.frame with `id` (`x`) and `id_type`. If `with_organism = TRUE`,
#'   a third `organism` column is added with a guess for the organism.
guess_feature_type <- function(x, with_organism = TRUE, summarize = TRUE) {
  regex <- list(
    ensgid  = "^ENS[A-Z]*G\\d+(\\.\\d+)?$",
    enstid  = "^ENS[A-Z]*?T\\d+(\\.\\d+)?$",
    refseq  = "^[NXW][CGMRP]_\\d+(\\.\\d+)?$",
    entrez  = "^\\d+$")

  bool <- sapply(regex, grepl, x)
  if (length(x) == 1L) bool <- t(bool) # convert to matrix
  nmatch <- rowSums(bool)
  type <- names(regex)[apply(bool, 1, function(vals) which(vals)[1])]
  type <- ifelse(nmatch == 1L, type, "ambiguous")
  type <- ifelse(nmatch == 0L, "unknown", type)

  is.bad <- type %in% c("ambiguous", "unknown")
  if (any(is.bad)) {
    warning(sum(is.bad), " identifiers were either ambiguous or unknown",
            immediate. = TRUE)
  }

  out <- tibble(
    id = x,
    feature_type = type)

  if (with_organism) {
    is.ens <- grepl("^ens", out[["feature_type"]])
    ens <- sub("^ENS", "", out[["id"]])
    is.human <- is.ens & grepl("^[TG]\\d+", ens)
    is.mouse <- is.ens & grepl("MUS[TG]\\d+", ens)
    out[["source_organism"]] <- ifelse(is.human, "Homo sapiens", "unknown")
    out[["source_organism"]] <- ifelse(is.mouse, "Mus musculus", out[["source_organism"]])
    unk <- out[["source_organism"]] == "unknown"
    if (any(unk)) {
      warning(sum(unk), " identifiers could not be matched to an organism",
              immediate. = TRUE)
    }
  }

  if (summarize) {
    idtypes <- unique(out[["feature_type"]])
    if (length(idtypes) != 1L) {
      stop("Ambiguous guess of feature type.")
    }
    out <- out[1, -1, drop=FALSE]
  }

  out
}
