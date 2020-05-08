# Command line usage:
# Rscript readings-density.R --(min|mean|max) plotfile.png data1.csv [data2.csv] [...]
#
# Where:
#   - --min, --mean or --max is the function to be applied to each patient
#   - plotfile.png is the path to write a density plot
#   - data#.csv files contain patient data on each line (without headers)

library(tidyr)
library(dplyr, warn.conflicts = interactive()) # keep silent when run on command line
library(ggplot2)
library(stringr)

main <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    action <- args[1]
    plotfile <- args[2]
    filenames <- args[c(-1, -2)]
    stopifnot(action %in% c("--min", "--mean", "--max"))
    action <- str_remove(action, "--") # remove hyphens
    stopifnot(str_detect(plotfile, "\\.(bmp|jpg|eps|ps|tex|pdf|jpeg|tiff|png|bmp|svg)$"))

    if (length(filenames) == 0) {
        process_many(file("stdin"), action, plotfile)
    } else {
        process_many(filenames, action, plotfile)
    }
}

process <- function(filename, action) {
    dat <- read.csv(file = filename, header = FALSE)

    if (action == "min") {
        values <- apply(dat, 1, min)
    } else if (action == "mean") {
        values <- apply(dat, 1, mean)
    } else if (action == "max") {
        values <- apply(dat, 1, max)
    }
    cat(values, sep = "\n")
    values
}


process_many <- function(files, action, outplot = NULL) {
    # Process many files and plot
    all_dat_list <- list()
    for (f in files) {
        values <- process(f, action)
        all_dat_list[[basename(f)]] <- tibble(!!action := values)
    }
    all_dat <- bind_rows(all_dat_list, .id = "filename")

    ggplot(all_dat, aes_string(action, colour = "filename")) +
        geom_density()
    if(!is.null(outplot)) {
        ggsave(outplot)
    }

    last_plot() # show the plot outside this function
}


if(interactive()) {
    # Interactive example
    process_many(fs::dir_ls("data", glob = "data/inflammation*"), "mean", "interactive.png")
} else {
    # Run when used from command line
    main()
}
