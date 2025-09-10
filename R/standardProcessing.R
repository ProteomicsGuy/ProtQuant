
#' @import methods
#' @import RColorBrewer
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom dplyr matches
#' @importFrom magrittr %>%
#' @importFrom limma normalizeMedianValues
#' @importFrom arrow read_parquet
#' @importFrom diann diann_matrix
#' @importFrom diann diann_maxlfq
#' @importFrom data.table setDT
#' @importFrom data.table fread
#' @importFrom data.table setDT
#' @importFrom utils read.table



globalVariables(c("Protein.Group", "Stripped.Sequence"))

#' @title Preprocessing for CQEs
#' @description Apply a standard pre-processing pipeline.
#' @details This function applies a standard pre-processing pipeline with some
#' useful customization options
#' @docType methods
#' @export
#' @rdname standardProcessing-methods
#' @param report_path Path to initial report (report.parquet)
#' @param pg.q diann_matrix protein group q-value threshold
#' @param parq_qval_cutoff Q-value threshold for filtering input report passed to diann_maxlfq
#' @param metadata_path Path to metadata table. Metadata table should contain a column 'Run' containing
#' the unique values from the 'Run' column in the initial report, and columns 'Group' and 'Replicate'
#' from which unique sample labels can be formed
#' @param pg_path Path to pg file
#' @param pr_path Path to pr file
#' @param min_peptides Minimum number of peptides per protein (used when filtering protein groups)
#' @param filter_contaminants Logical to indicate whether protein groups starting 'Cont_' should be filtered out
#' @param apply_normalization Logical to indicate whether limma::normalizeMedianValues should be applied
#' @param lfq_col Column from diann report file to use in quantitation
#' @returns A dataframe containing annotation and values per protein group for each sample.
setGeneric("standardProcessing", function(report_path, pg.q = 0.01, parq_qval_cutoff = 0.01, metadata_path,
                                          lfq_col = "Precursor.Quantity", pg_path, pr_path, min_peptides = 2,
                                          filter_contaminants = TRUE, apply_normalization = TRUE)
  standardGeneric("standardProcessing"))

#' @rdname standardProcessing-methods
setMethod("standardProcessing", signature(report_path = "character", metadata_path = "character"),
          function(report_path, pg.q = 0.01, parq_qval_cutoff = 0.01, metadata_path,
                   lfq_col = "Precursor.Quantity", pg_path, pr_path, min_peptides = 2,
                   filter_contaminants = TRUE, apply_normalization = TRUE){

            # Read in initial report
            message("Reading report ... ", appendLF = FALSE)
            df_parq <- read_parquet(report_path)
            df_parq$File.Name <- df_parq$Run
            message("Done.")

            # Use the report to produce the precursor and protein group matrices:
            # Create DIANN R precursor matrix
            message("Creating precursor matrix with diann_matrix ... ", appendLF = FALSE)
            precursors <- diann_matrix(df_parq, pg.q = pg.q, quantity.header = "Precursor.Quantity")
            precursors <- as.data.frame(precursors)
            message("Done.")

            # Create DIANN R protein group matrix
            message("Creating protein group matrix with diann_maxlfq ... ", appendLF = FALSE)
            protein.groups <- diann_maxlfq(df_parq[df_parq$Q.Value <= parq_qval_cutoff & df_parq$PG.Q.Value <= parq_qval_cutoff,],
                                           group.header="Protein.Group",
                                           id.header = "Precursor.Id",
                                           quantity.header = lfq_col)
            protein.groups <- as.data.frame(protein.groups)
            message("Done.")

            # Re-order and rename based on metadata
            message("Processing metadata and applying filters ... ", appendLF = FALSE)
            metadata <- read.table(metadata_path, stringsAsFactors = FALSE, header = TRUE, sep = "\t")
            # Check that Run column in metadata contains column names of the precursor and protein matrices
            stopifnot("Run entries in metadata do not match column names of precursor matrix" = all(sort(metadata$Run) == sort(colnames(precursors))))
            stopifnot("Run entries in metadata do not match column names of protein group matrix" = all(sort(metadata$Run) == sort(colnames(protein.groups))))

            # Make unique sample labels
            metadata$sample_label = paste0(metadata$Group, "_", metadata$Replicate)

            # Now we know entries are the same, reorder columns of precursor and protein group matrices to match order in metadata
            protein.groups <- protein.groups[, metadata$Run]
            precursors <- precursors[, metadata$Run]

            # Rename columns of precursor and protein group matrices with sample labels
            colnames(protein.groups) = metadata$sample_label
            colnames(precursors) = metadata$sample_label

            setDT(protein.groups, keep.rownames ="Protein.Group")
            setDT(precursors, keep.rownames ="Precursor.Id")

            # Load in whole pg and pr file to obtain extra information to append to filtered precursor and protein group arrays.
            pg <- fread(pg_path)
            pr <- fread(pr_path)

            pginfo <- pg[, 1:4]
            prinfo <- pr[, 1:10]

            # Create precursor and protein tables with extra information, and merge
            pg_merge <- merge(pginfo, protein.groups, by="Protein.Group")
            pr_merge <- merge(prinfo, precursors, by="Precursor.Id")

            # Filter protein groups by number of peptides per protein
            pep_per_prot = pr_merge %>% dplyr::select(Protein.Group, Stripped.Sequence) %>% dplyr::distinct() %>%
              dplyr::pull(Protein.Group) %>% table()
            pep_to_keep = names(pep_per_prot[pep_per_prot >= min_peptides])
            pg_merge_filtered = dplyr::filter(pg_merge, Protein.Group %in% pep_to_keep)

            # Filter contaminants
            if(filter_contaminants){
              pg_merge_filtered = dplyr::filter(pg_merge_filtered, !grepl("Cont_", Protein.Group))
            }

            message("Done.")

            # Normalization
            if(apply_normalization){

              message("Applying normalisation ... ", appendLF = FALSE)

              # Log transform
              protein.matrix = dplyr::select(pg_merge_filtered, matches(metadata$sample_label))
              protein.matrix.annot = dplyr::select(pg_merge_filtered, -matches(metadata$sample_label))

              protein.matrix.log <- log2(protein.matrix + 1)
              #rownames(protein.matrix.log) = pg_merge_filtered$Protein.Group

              # Apply normalize median values
              protein.matrix.norm <- limma::normalizeMedianValues(protein.matrix.log)
              # Append annotation
              pg_merge_filtered = cbind(protein.matrix.annot, protein.matrix.norm)

              message("Done.")
            }

            message("Preprocessing finished.")

            return(pg_merge_filtered)
          }
)


