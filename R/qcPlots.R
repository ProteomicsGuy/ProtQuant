
#' @title QC plots for CQE analysis
#' @description Produce useful QC plots following pre-processing.
#' @details This function produces boxplots, density plots, correlation plots, PCA plots.
#' @docType methods
#' @export
#' @rdname qcPlots-methods
#' @param metadata_path Path to metadata table. Metadata table should contain a column 'Run' containing
#' the unique values from the 'Run' column in the initial report, and columns 'Group' and 'Replicate'
#' from which unique sample labels can be formed. Should be the same metadata file used in
#' the call to \code{standardProcessing}.
#' @param expr.df expression data frame as produced by \code{standardProcessing}. Annotation columns may need to be removed.
#' @param output_dir Path to directory where figures will be put.
#' @param filename_ext Spliced into into output file names e.g. if filename_ext = "MEDNORM", output files will
#' be \code{log2protboxplotsMEDNORM_all.jpeg} etc. Defaults to filename_ext = "" to leave the base forms unchanged.
#' @param columns_per_plot Used for boxplots
setGeneric("qcPlots", function(metadata_path, expr.df, output_dir, filename_ext = "", columns_per_plot)
  standardGeneric("qcPlots"))

#' @rdname qcPlots-methods
setMethod("qcPlots", signature(metadata_path = "character", expr.df = "data.frame"),
          function(metadata_path, expr.df, output_dir, filename_ext = "", columns_per_plot){

            # Read metadata
            metadata <- read.table(metadata_path, stringsAsFactors = FALSE, header = TRUE, sep = "\t")

            # 1. Define the conditions and colors
            conds <- as.factor(metadata$Group)  # Assuming metadata$Group corresponds to conditions
            if (length(unique(conds)) > 9) {
              cond_colours <- colorRampPalette(brewer.pal(9, "Set1"))(length(unique(conds)))
            } else {
              cond_colours <- brewer.pal(length(unique(conds)), "Set1")
            }
            names(cond_colours) <- levels(conds)

            # 2. Determine the number of plots needed, using columns_per_plot argument (maximum columns per plot)
            num_columns <- ncol(expr.df)
            num_plots <- ceiling(num_columns / columns_per_plot)

            # 3. Save all plots in one JPEG file
            output_file_boxplots <- file.path(output_dir, paste0("log2protboxplots", filename_ext, "_all.jpg"))
            jpeg(filename = output_file_boxplots, width = 1600, height = 1200 * num_plots, res = 150)

            # Adjust layout to stack multiple plots vertically
            par(mfrow = c(num_plots, 1))  # One column, multiple rows
            par(mar = c(12, 6, 2, 2) + 0.1)  # Increased bottom and left margins for better readability

            # 4. Set consistent y-axis limits across all plots (to keep them the same)
            y_min <- min(expr.df, na.rm = TRUE)  # Get the minimum value for all columns
            y_max <- max(expr.df, na.rm = TRUE)  # Get the maximum value for all columns


            # 5. Loop through each chunk of columns and create plots
            for (i in seq_len(num_plots)) {
              # Get the columns for this plot
              start_col <- (i - 1) * columns_per_plot + 1
              end_col <- min(i * columns_per_plot, num_columns)
              cols_to_plot <- expr.df[, start_col:end_col, drop = FALSE]

              # Map the condition colors for the current subset of columns
              col_subset <- cond_colours[conds[start_col:end_col]]

              # Create the boxplot with consistent y-axis limits and tick marks
              boxplot(cols_to_plot,
                      col = col_subset,                       # Assign condition-based colors
                      las = 2,                               # Rotate x-axis labels
                      cex.axis = 1.2,                        # Increase tick label size
                      ylab = "log2 Intensity",              # Y-axis label
                      cex.lab = 1.5,                         # Increase axis label font size
                      ylim = c(y_min, y_max),               # Same y-axis limits for all plots
                      xaxt = 's')                            # Automatically plot x-axis labels

              # If column names are long, make sure the labels don't get cut off
              # this overwrites existing axis,so comment
              # axis(1, at = seq(1, ncol(cols_to_plot), by = 5), labels = colnames(cols_to_plot)[seq(1, ncol(cols_to_plot), by = 5)], las = 2, cex.axis = 0.8)
            }

            dev.off()  # Save the file

            cat("Boxplots saved to ", output_file_boxplots, "\n")
            expr.df

            #2. Density plot
            oldmar <- par("mar")
            par(mar = c(5, 4, 4, 6))

            plotDensities(expr.df, legend=FALSE)
            legend("topright", legend = colnames(expr.df),
                   lty = 1 , cex=0.5, ncol = 2, col=1:ncol(expr.df), xpd = TRUE, inset = c(-0.1, 0))

            output_file_density <- file.path(output_dir, paste0("log2protdensity", filename_ext, "_plot.jpg"))
            jpeg(filename = output_file_density, width = 1600, height = 1200, res = 150)
            plotDensities(expr.df, legend=FALSE)
            legend("topright", legend = colnames(expr.df),
                   lty = 1 , cex=0.75, ncol = 2, col=1:ncol(expr.df), xpd = TRUE)

            dev.off()
            par(mar = oldmar)

            cat("Density plots saved to ", output_file_density, "\n")
            #3. Correlation plots
            # Calculate correlation matrix
            m_all = cor(expr.df, use="pairwise.complete.obs")

            # Save the correlation plot to a JPEG file with high resolution
            output_file_correlation <- file.path(output_dir, paste0("log2protcorrelation", filename_ext, "_plot.jpg"))
            jpeg(filename = output_file_correlation , width = 8000, height = 8000, res = 300)

            # Create the correlation plot
            corrplot(m_all,
                     type = "upper",                    # Only upper triangle of the matrix
                     method = "color",                   # Color-coded correlation values
                     cl.lim = c(0.65, 1),                # Color limit for correlation range
                     is.corr = FALSE,                    # Treat as a general matrix, not a correlation matrix
                     tl.col = "black",                   # Text label color for column/row names
                     tl.cex = 0.8,                       # Increase text size for column/row labels
                     addCoef.col = NULL,                 # Do not add correlation coefficient text
                     number.cex = 0.6, # Adjust the size of correlation values on the plot
                     cl.cex = 3)

            dev.off()  # Close the graphics device and save the file

            cat("Correlation plot saved to ", output_file_correlation, "\n")

            #4. PCA Plot
            pca_data = prcomp(t(na.omit(expr.df)))
            pca_data_perc = round(100 * pca_data$sdev^2 / sum(pca_data$sdev^2), 1)
            df_pca_data = data.frame(PC1 = pca_data$x[, 1], PC2 = pca_data$x[, 2], PC3 = pca_data$x[, 3], sample = colnames(expr.df), condition = metadata$Group)

            # Function to find convex hull for each condition (to highlight areas)
            find_hull = function(df_pca_data) df_pca_data[chull(df_pca_data$PC1, df_pca_data$PC2),]
            hulls = ddply(df_pca_data, "condition", find_hull)

            # Define the number of unique conditions for better color mapping
            num_conditions <- length(unique(df_pca_data$condition))
            condition_palette <- distinctColorPalette(num_conditions)  # Generates as many unique colors as needed

            # Create the PCA plot with the improvements
            PCA <- ggplot(df_pca_data, aes(PC1, PC2, color = condition, fill = condition)) +
              geom_point(linewidth = 3, shape = 21, stroke = 1, alpha = 0.7) +  # Make points cleaner with transparency and outline
              labs(x = paste0("PC1(", pca_data_perc[1], "%)"),
                   y = paste0("PC2(", pca_data_perc[2], "%)"),
                   title = "Protein Group-based") +
              stat_ellipse() +  # Add ellipse for each condition
              theme_minimal(base_size = 16) +  # Use minimal theme with larger text
              scale_color_manual(values = condition_palette) +  # Use custom color palette for the conditions
              scale_fill_manual(values = condition_palette) +  # Ensure fill matches the color
              theme(
                panel.grid = element_blank(),  # Remove gridlines
                axis.ticks = element_line(size = 1, color = "black"),  # Add tick marks to axis borders
                axis.line = element_line(size = 1, color = "black"),  # Make axis lines clearer
                axis.text = element_text(size = 14, color = "black"),  # Make axis labels larger
                axis.title = element_text(size = 16, color = "black"),  # Increase axis title size
                legend.position = "right",  # Move legend to the right of the plot
                legend.title = element_text(size = 14),  # Adjust legend title size
                legend.text = element_text(size = 12)  # Adjust legend text size
              )
            # Print the PCA plot
            PCA
            # Save the plot as a high-resolution jpeg
            output_file_PCA <- file.path(output_dir, paste0("log2protPCA", filename_ext, "_QC.jpg"))
            jpeg(filename = output_file_PCA, width = 2200, height = 1800, res = 300)
            print(PCA)
            dev.off()  # Close the graphics device and save the file

            cat("PCA plot saved as ", output_file_PCA, "\n")

            # 5. Protein Count plot
            # Step 1: Count non-NA values for each column in expr.df
            counts_per_column <- colSums(!is.na(expr.df))

            # # Step 2: Create a dataframe to map counts to metadata
            # counts_df <- data.frame(
            #   Rename = names(counts_per_column),
            #   Counts = counts_per_column
            # )
            #
            # # Step 3: Merge counts with metadata
            # merged_df <- merge(counts_df, metadata, by.x = "Rename", by.y = "Rename")

            # Step 2: Create a dataframe to map counts to metadata
            counts_df <- data.frame(
              sample_label = names(counts_per_column),
              Counts = counts_per_column
            )

            # Step 3: Merge counts with metadata
            merged_df <- merge(counts_df, metadata, by.x = "sample_label", by.y = "sample_label")

            # Step 4: Calculate the mean counts for each Group
            group_means <- merged_df %>%
              group_by(Group) %>%
              summarise(
                mean_count = mean(Counts, na.rm = TRUE)  # Calculate mean, handle NAs
              ) %>%
              ungroup()  # Ensure it is no longer grouped for later use

            # Verify the structure of group_means
            print(group_means)
            # write.table(group_means,"Z:/Trost-group/Andy_F/2_Connect_FTMA_AF_2025/Multiomics/Proteomics/2025_UbiSCAPE/lib_free_search/Analysis/Filtered_protein_count_group_means.df.txt", sep = "\t")
            group_means_output_file <- file.path(output_dir, paste0("Filtered_protein_count_group_means", filename_ext, ".df.txt"))
            write.table(group_means, group_means_output_file, sep = "\t")
            cat("Mean counts per group saved as ", group_means_output_file, "\n")
            histoGRAM <- ggplot() +
              # Bars for the counts of non-NA values per group (mean counts)
              geom_bar(data = group_means, aes(x = Group, y = mean_count),
                       stat = "identity", fill = "skyblue", color = "black",
                       width = 0.6, alpha = 0.7) +
              # Dots for each replicate, aligned with the bar on the histogram
              geom_point(data = merged_df, aes(x = Group, y = Counts, color = Group),
                         size = 3, alpha = 0.8, position = position_jitter(width = 0.2)) +
              # Add labels for axes and title
              labs(x = NULL, y = "Protein Groups",
                   title = paste(filename_ext, "Proteins")) +
              # Custom theme to remove gridlines, thicken axes, and adjust other aesthetics
              theme_minimal(base_size = 14) +
              theme(
                panel.grid = element_blank(),  # Remove gridlines
                axis.ticks = element_line(size = 1.5, color = "black"),  # Thicker axis ticks
                axis.line = element_line(size = 1.5, color = "black"),  # Thicker axis lines
                axis.text = element_text(size = 14, color = "black"),  # Make axis labels larger
                axis.title = element_text(size = 16, color = "black"),  # Larger axis titles
                axis.text.x = element_text(angle = 45, hjust = 1),  # Diagonal x-axis labels
                axis.text.y = element_text(angle = 0, hjust = 1),  # Horizontal y-axis labels
                legend.position = "none",  # Remove legend if not necessary
                plot.title = element_text(hjust = 0.5)  # Center the plot title
              ) +
              # Adjust y-axis to ensure enough space for bars, and change y-ticks to every 1000
              scale_y_continuous(expand = c(0, 0),
                                 limits = c(0, max(group_means$mean_count) * 1.5),
                                 breaks = seq(0, max(group_means$mean_count), by = 1000))

            output_file_histogram <- file.path(output_dir, paste0("Protein_Count", filename_ext, "_post_filtering.jpeg"))
            ggsave(output_file_histogram, plot = histoGRAM,
                   width = 10, height = 8, dpi = 300, units = "in",
                   device = "jpeg")
            cat("Histogram saved as ", output_file_histogram, "\n")
          }
)

