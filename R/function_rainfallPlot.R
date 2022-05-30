#' @title Rainfall plot.
#' @description Visualize the IMD, segments and putative kataegis foci using a rainfall plot.
#'
#' Y-axis represents the 5' intermutation distance (IMD) of each genomic variant in a log10-scale.
#' X-axis represent the IMD for each variant. Putative kataegis foci are highlighted by blue backgrounds.
#'
#' Variants within kataegis foci are bold. Color represent the mutational context whilst horizontal lines represent
#' the detected segments (using the mean IMD). Vertical lines depict the detected changepoints.
#'
#' @param kd (\link[katdetectr]{KatDetect}): KatDetect object.
#' @param showSequence (character): Which sequence(s) should be visualized? Choice of: 'All', 'Kataegis', c('seqname1', 'seqname2').
#' @param showKataegis (logical): Highlight putative kataegis foci and underlying genomic variants?
#' @param showSegmentation (logical): Show changepoints and mean IMD of each segment?
#'
#' @examples
#' syntheticData <- generateSyntheticData(nBackgroundVariants = 200)
#' kd <- detectKataegis(syntheticData)
#'
#' # Visualize the IMD of the genomic variants by constructing a rainfall plot
#' katdetectr::rainfallPlot(kd)
#'
#' # Show the chromosomes which contain one or more kataegis foci
#' katdetectr::rainfallPlot(kd, showSequence = "Kataegis")
#'
#' # Show only chromosome 1 and 2
#' katdetectr::rainfallPlot(kd = kd, showSequence = c("chr1", "chr2"))
#'
#' # Display changepoints and mean IMD per segment
#' katdetectr::rainfallPlot(kd = kd, showSequence = c("chr1", "chr2"), showSegmentation = TRUE)
#'
#' @return (\link[ggplot2]{ggplot}): Returns rainfall plot.
#'
#' @author Daan Hazelaar
#' @author Job van Riet
#' @export
rainfallPlot <- function(kd, showSequence = 'All', showKataegis = TRUE, showSegmentation = FALSE){

    # Input validation ---------------------------------------------------------
    .validateInputRainfallPlot(kd, showSequence, showKataegis, showSegmentation)

    # select the sequences requested by the user -------------------------------
    selectedSequences <- .selectSequences(kd, showSequence)

    # Convert to ggplot2-friendly format ---------------------------------------
    plotDataVariants <- kd |>
        katdetectr::getGenomicVariants() |>
        .convertVariantsToGGplotFormat() |>
        dplyr::filter(.data$seqnames %in% selectedSequences)

    plotDataKataegis <- kd |>
        katdetectr::getKataegisFoci() |>
        tibble::as_tibble() |>
        dplyr::filter(.data$seqnames %in% selectedSequences)

    plotDataSegments <- kd |>
        katdetectr::getSegments() |>
        tibble::as_tibble() |>
        dplyr::filter(.data$seqnames %in% selectedSequences)

    # generate plot ------------------------------------------------------------
    p <- .generateRainfallPlot(plotDataVariants, plotDataKataegis, plotDataSegments, showKataegis, showSegmentation)

    return(p)
}

.validateInputRainfallPlot <- function(kd, showSequence, showKataegis, showSegmentation){

    checkmate::assertClass(kd, classes = 'KatDetect')
    checkmate::assertCharacter(showSequence)
    checkmate::assertLogical(showKataegis)
    checkmate::assertLogical(showSegmentation)

    # Check if requested sequence is in the given data.
    sequencesInData <- kd |> katdetectr::getGenomicVariants() |> GenomeInfoDb::seqlevelsInUse()
    if(base::any(showSequence != 'All') & base::any(showSequence != 'Kataegis')){
        if(!base::all(showSequence %in% sequencesInData)){
            base::stop('Requested sequences are not present in the data.')
        }
    }
}

.selectSequences <- function(kd, showSequence){
    # Check if any chromosome harbors putative kataegis, else show all chromosome by default.
    if(base::length(showSequence) == 1){
        if(showSequence == 'All'){
            selectedSequences <- GenomeInfoDb::seqlevelsInUse(kd@genomicVariants)
        } else if(base::length(kd@kataegisFoci) == 0 & showSequence == 'Kataegis'){
            selectedSequences <- GenomeInfoDb::seqlevelsInUse(kd@genomicVariants)
        } else if(showSequence == 'Kataegis'){
            selectedSequences <- kd |>
                katdetectr::getKataegisFoci() |>
                tibble::as_tibble() |>
                dplyr::distinct(.data$seqnames) |>
                dplyr::pull() |>
                base::as.character()
        } else {selectedSequences <- showSequence}
    } else {selectedSequences <- showSequence}

    return(selectedSequences)
}

.convertVariantsToGGplotFormat <- function(genomicVariants){

    plotDataVariants <- genomicVariants |>
        tibble::as_tibble() |>
        dplyr::mutate(
            IMD = tidyr::replace_na(.data$IMD, -1),
            variantType = dplyr::case_when(
                ref == "C" & alt == "T" ~ "C>T",
                ref == "C" & alt == "G" ~ "C>G",
                ref == "C" & alt == "A" ~ "C>A",
                ref == "T" & alt == "C" ~ "T>C",
                ref == "T" & alt == "A" ~ "T>A",
                ref == "T" & alt == "G" ~ "T>G",
            ),
            variantType = base::ifelse(base::is.na(.data$variantType), "Other", .data$variantType)
        )

    return(plotDataVariants)
}

.generateRainfallPlot <- function(plotDataVariants, plotDataKataegis, plotDataSegments, showKataegis, showSegmentation){

    p <- plotDataVariants |>
        ggplot2::ggplot(ggplot2::aes(x = .data$variantID, y = .data$IMD, group = .data$seqnames)) +

        # Plot 5' IMDs.
        ggplot2::geom_point(
            mapping = ggplot2::aes(fill = .data$variantType, group = .data$seqnames, alpha = .data$putativeKataegis, size = .data$putativeKataegis),
            color ='black', shape = 21) +

        # Color point on mutational type.
        ggplot2::scale_fill_manual(
            name = 'Mutational Type',
            values = c('C>T' = '#E82A25','C>G' ='#010101' , 'C>A' = '#25BDEE', 'T>C' = '#e6e600', 'T>A' = '#CAC9C9', 'T>G' = '#ECC9C8', 'Other' = 'grey50'),
            guide = ggplot2::guide_legend(title.position = 'top', title.hjust = 0.5, nrow = 1, keywidth = 0.5, keyheight = 0.5)) +

        # Labs.
        ggplot2::labs(x = base::sprintf('Sample: %s\nTotal variants: %s, putative kataegis foci: %s', base::unique(plotDataVariants$sampleNames), base::nrow(plotDataVariants), base::nrow(plotDataKataegis)), y = 'Intermutation Distance') +

        # Set alpha and scale of IMD's based on flagged within kataegis foci.
        ggplot2::scale_alpha_manual(values = c('TRUE' = base::ifelse(showKataegis, 1, 0.5), 'FALSE' = .3), guide = 'none') +
        ggplot2::scale_size_manual(values = c('TRUE' = base::ifelse(showKataegis, 2, 1), 'FALSE' = 1), guide = 'none') +

        # Options - Axis.
        ggplot2::scale_x_continuous(expand = c(0, 0)) +
        ggplot2::scale_y_continuous(
            trans = scales::pseudo_log_trans(),
            breaks = c(-1, 0, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000),
            labels = c('NA','0bp', '10bp', '100bp', '1000bp', '0.01Mbp', '0.1Mbp', '1Mbp', '10Mbp', '100Mbp'),
            expand = c(0, 0)) +

        ggplot2::labs(y = 'IMD') +

        # Split on chromosomes.
        ggplot2::facet_grid(.~.data$seqnames, space = 'free_x', scales = 'free_x', drop = TRUE) +

        # Theme.
        ggplot2::theme(
            panel.spacing.x = ggplot2::unit(0.05, 'lines'),
            legend.position = 'bottom',
            legend.direction = 'horizontal',
            text = ggplot2::element_text(size = 9, family = 'Helvetica', face = 'bold'),
            axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            axis.title.y = ggtext::element_textbox_simple(size = 8, orientation = 'left-rotated', width = NULL, halign = .5),
            panel.grid.major.x = ggplot2::element_blank(),
            panel.grid.minor.x = ggplot2::element_blank(),
            panel.grid.major.y = ggplot2::element_line(colour = '#E5E5E5', linetype = 'dotted'),
            panel.grid.minor.y = ggplot2::element_blank(),
            panel.background = ggplot2::element_rect(fill = NA, colour = 'black'),
            panel.border = ggplot2::element_rect(fill = NA, colour = NA),
            strip.background = ggplot2::element_blank(),
            legend.key = ggplot2::element_blank(),
            strip.text = ggplot2::element_text(angle = 60)
        )

    # Highlight kataegis foci.
    if(showKataegis & base::nrow(plotDataKataegis) != 0){
        p <- p + ggplot2::geom_rect(
            data = plotDataKataegis,
            inherit.aes = FALSE,
            ggplot2::aes(xmin = .data$firstVariantID, xmax = .data$lastVariantID, ymin = -1, ymax = Inf, group = .data$seqnames),
            fill = '#0080FF30')
    }

    # show means and segments
    if(showSegmentation){
        p <- p +
            ggplot2::geom_segment(
                data = plotDataSegments,
                mapping = ggplot2::aes(x = .data$firstVariantID, xend = .data$lastVariantID + 1, y = .data$meanIMD, yend = .data$meanIMD),
                col = 'black', lty = 'solid', na.rm = TRUE) +
            ggplot2::geom_segment(
                data = plotDataSegments,
                mapping = ggplot2::aes(x = .data$firstVariantID, xend = .data$firstVariantID, y = -1,yend = Inf),
                col = 'black', lty = 'dotted')
    }

    return(p)
}
