
globalVariables(
  c(
    "BS_FQR",
    "BS_FQR_count",
    "BS_FQR_percentage",
    "FQR25",
    "FQR25_count",
    "INV_FQR",
    "P.Value",
    "Precision",
    "Protein.Names",
    "adj.P.Val",
    "expected_FC",
    "expected_fc",
    "expected_logFC_minus25",
    "expected_logFC_plus25",
    "logFC",
    "significant",
    "significant_count",
    "total_proteins"
  )
)


#' @title Calculate CQEs Metrics
#' @description Calculate metrics to assess controlled quantitative experiments
#' @param limma_df limma output (dataframe)
#' @param species names of species compared (character vector)
#' @param expected_fcs expected fold changes for each species (numeric vector)
#' @param species_protein_annotation a tibble with 2 columns,
#' the first column must contain ids matching to the rownames of the limma df (e.g. protein groups)
#' the second column must identify the species which the ID comes from. This can either be specified directly or supplied at the end of a string, e.g. Uniprot name.
#' the species must match with those provided in species vector
#'
#' @returns A tibble summarising the total number of proteins, "INV_FQR", "Precision" and "BS_FQR" per species.
#' @export

calculate_cqe_metrics <- function(limma_df, species, expected_fcs, species_protein_annotation){

  # Do some checks and throw some errors if the inputs aren't as expected

  # limma_df - should be a dataframe
  assertthat::assert_that(is.data.frame(limma_df))

  #species - should be a character vector of any length
  assertthat::assert_that(is.character(species))
  #expected_fcs - should be a numeric vector of any length
  assertthat::assert_that(is.numeric(expected_fcs))
  # species & expected_fcs - vector lengths should match
  assertthat::assert_that(length(expected_fcs) == length(species),msg = "`species` and `expected_fcs` must be the same length")

  # species protein annotation - currently should be a tibble with 2 columns
  assertthat::assert_that(tibble::is.tibble(species_protein_annotation))
  # the first must contain protein groups or ids matching to the limma df
  # I think we could check for this by seeing if all of the limma base::rownames can be found in the first column
  assertthat::assert_that(all(base::rownames(limma_df) %in% dplyr::pull(species_protein_annotation,1)),
                          msg = "The first column of `species_protein_annotation` must contain all the base::rownames of `limma_df` output")

  # the second must contain strings which contain all the species found in expected_fcs
  assertthat::assert_that(all(base::rownames(limma_df) %in% dplyr::pull(species_protein_annotation,1)),
                          msg = "The first column of `species_protein_annotation` must contain all the base::rownames of `limma_df` output")

  assertthat::assert_that(
    all(stringr::str_to_upper(species) %in% dplyr::pull(dplyr::distinct(tibble::tibble(species = dplyr::pull(species_protein_annotation[2])) %>%
                                                                          dplyr::mutate(species = stringr::str_to_upper(stringr::str_extract(species,pattern = "[a-zA-Z]+$")
                                                                          )
                                                                          )
    )
    )
    ),
    msg = "the second column of `species_protein_annotation` must contain all the values given in `species`.
                          This can either be in the form of uniprot_names (e.g. TVAL3_HUMAN) or provided directly (e.g.HUMAN)"
  )

  # I think this handles it a little better it should be permissive that just need the string to contain the species at the end (assuming a non-letter character before)
  # it should allow for people to either give protein
  # There could be a neater input option for this - can we go and fetch protein names remotely instead? Not sure of exisitng packages for this?


  # Setup

  # Make up the expected species tibble and calculate acceptable thresholds
  tibble::tibble(species = stringr::str_to_upper(species),
                 expected_FC = expected_fcs
  ) %>%
    dplyr::mutate(expected_logFC = log2(expected_FC),
                  expected_logFC_minus25 = log2(expected_FC*0.75),
                  expected_logFC_plus25 = log2(expected_FC*1.25)
    ) -> species_expected_fc

  # Rename column names of species_protein_annotation
  species_protein_annotation %>%
    dplyr::rename_with(~c("Protein.Group","Protein.Names")) -> species_protein_annotation

  # format the limma df
  limma_df %>%
    tibble::rownames_to_column("Protein.Group") %>%
    dplyr::select(Protein.Group,logFC,P.Value,adj.P.Val) %>%

    dplyr::left_join(species_protein_annotation, by = "Protein.Group") %>%

    dplyr::mutate(species = stringr::str_to_upper(
      stringr::str_extract(Protein.Names,pattern = "[a-zA-Z]+$")
    )
    ) %>%
    dplyr::right_join(species_expected_fc, by = "species") -> metrics_input

  # Calculate Required inputs for summary metrics

  metrics_input %>%

    # calculate the total number of proteins per species
    dplyr::group_by(species) %>%
    dplyr::count(name = "total_proteins") %>%
    dplyr::ungroup() %>%
    dplyr::left_join(metrics_input, by = "species") %>%


    # FQR25: check if the logFC is outside the limits - TRUE means it is, FALSE means its in the expected
    dplyr::mutate(FQR25 = dplyr::if_else(logFC > expected_logFC_plus25 | logFC < expected_logFC_minus25,
                                         TRUE, FALSE)
    ) %>%

    # Significant hits
    dplyr::mutate(significant = dplyr::if_else(adj.P.Val < 0.05, "Sig","NS")) %>%

    # BSFQR: Significant hits in the opposite expected direction to ground truth
    dplyr::mutate(BS_FQR = dplyr::case_when(
      significant == "Sig" & expected_logFC ==0 & logFC > log2(1.5) ~ TRUE,
      significant == "Sig" & expected_logFC ==0 & logFC < -log2(1.5) ~ TRUE,
      significant == "Sig" & expected_logFC > 0 & logFC < 0 ~ TRUE,
      significant == "Sig" & expected_logFC < 0 & logFC > 0 ~ TRUE,
      TRUE ~ FALSE
    )
    ) -> metrics

  # Calculate the summary percentages to return

  metrics %>%
    dplyr::group_by(species,FQR25,.add = TRUE) %>%
    dplyr::count(name = "FQR25_count",total_proteins) %>%
    dplyr::ungroup() %>%
    dplyr::filter(FQR25 == TRUE) %>%
    dplyr::mutate(INV_FQR = 100 - ((FQR25_count /total_proteins)*100)) %>%
    dplyr::select(species,total_proteins, INV_FQR) -> INV_FQR_summary

  metrics %>%
    dplyr::group_by(species) %>%
    dplyr::count(significant,name = "significant_count",total_proteins) %>%
    dplyr::ungroup() %>%
    dplyr::filter(significant == "Sig") %>%
    dplyr::mutate(Precision = significant_count / total_proteins *100) %>%
    dplyr::select(species,total_proteins,Precision) %>%
    # add a catch for precision value to return NA if expected FC is 1
    dplyr::left_join(expected_fc,by = "species") %>%
    dplyr::mutate(Precision = dplyr::if_else(expected_FC == 1, NA, Precision)) %>%
    dplyr::select(species,total_proteins, Precision) -> Precision_summary

  metrics %>%
    dplyr::group_by(species) %>%
    dplyr::count(BS_FQR,name = "BS_FQR_count",total_proteins) %>%
    dplyr::ungroup() %>%
    dplyr::filter(BS_FQR == TRUE) %>%
    dplyr::mutate(BS_FQR_percentage = BS_FQR_count / total_proteins *100) %>%
    dplyr::select(species,total_proteins, BS_FQR_percentage) %>%
    # Add a catch for when no proteins are the wrong direction
    dplyr::mutate(BS_FQR_percentage = dplyr::if_else(is.na(BS_FQR_percentage),0, BS_FQR_percentage)) %>%
    dplyr::rename(BS_FQR = BS_FQR_percentage) -> BS_FQR_summary


  INV_FQR_summary %>%
    dplyr::full_join(Precision_summary, by = c("species","total_proteins")) %>%
    dplyr::full_join(BS_FQR_summary, by = c("species","total_proteins")) -> summary_metrics

  return(summary_metrics)

}
