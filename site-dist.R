#!/usr/bin/env Rscript

library(argparse)
library(tidyverse)
library(readr)
library(ggplot2)
library(dplyr)
library(ggthemes)
library(cowplot)


desc <- ""

parser <-
    ArgumentParser(description = desc,
                   formatter_class = "argparse.RawTextHelpFormatter")

parser$add_argument("file1", type = "character", help = "mutants info file 1")
parser$add_argument("file2", type = "character", help = "mutants info file 2")

parser$add_argument("-k", type = "integer", default=21, help = "k-mer size")
parser$add_argument("-m", type = "integer", default=5, help = "m-mer size")
parser$add_argument("-s", type = "integer", default=5, help = "scale")

args <- parser$parse_args()

# --------------------------------------------------------------------

df1 <- read_tsv(
    args$file1,
    comment = "#",
    col_types = list(hash = "i", mpos = "i", mpos2 = "i")
) 

df2 <- read_tsv(
    args$file2,
    comment = "#",
    col_types = list(hash = "i", mpos = "i", mpos2 = "i")
) 

# --------------------------------------------------------------------

theme <- theme_bw() +
    theme(
        panel.border =  element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        
        axis.text = element_text(colour = "grey10", size = 10),
        axis.title = element_text(colour = "grey10", size = 12),
        
        axis.line = element_line(color = "grey10", size = 0.5),
        axis.ticks = element_line(size = 0.6),
        
        strip.background = element_rect(
            colour = "white",
            fill = "grey90",
            size = 0.2
        ),
        strip.text = element_text(size = 10),
        
        legend.position = "right",
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(size = 12),
        legend.margin = margin(0, 0, 0, 0, unit = "cm"),
        
        text = element_text(size = 10),
        plot.title = element_text(size = 15),
        plot.subtitle = element_text(size = 12),
        plot.caption = element_text(size = 11),
    )

df1p <- df1 %>% filter(status > 0)
pos <- df1p %>% distinct(pos) %>% nrow()
total <- df1 %>% distinct(pos) %>% nrow()

bins <- max(df1p$mpos)
p1 <- ggplot(df1p, aes(mpos)) +
    geom_histogram(bins = bins, fill = "#56b4e9", color = "grey30") +
    
    xlab("Position") +
    ylab("Count") +
    labs(title = "Single substitution", 
         subtitle = paste0("K-mers with mutants of same hashes: ", round(pos/total*100,2), "%"),
         caption =  paste0("k=",args$k, ", m=", args$m, ", scale=", args$s, "\n"),
         ) +
    
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme

df2p <- df2 %>% filter(status > 0)
pos2 <- df2p %>% distinct(pos) %>% nrow()
total2 <- df2 %>% distinct(pos) %>% nrow()

df22 <- df2p %>%
    group_by(mpos, mpos2) %>%
    summarise(count = n())

p2 <- ggplot(df22, aes(mpos, mpos2)) +
    geom_point(aes(size = count, color = count)) +
    scale_color_gradient_tableau(palette = "Orange", guide = "legend") +
    scale_x_continuous(position = "top") +
    xlab("Position 1") +
    ylab("Position 2") +
    labs(title = "Two substitutions", 
         subtitle = paste0("K-mers with mutants of same hashes: ", round(pos2/total2*100,2), "%"),
         caption =  paste0("K=",args$k, ", m=", args$m, ", scale=", args$s, "\n"),
    ) +
    # guides(shape = guide_legend(reverse = TRUE)) +
    theme +
    theme(
        panel.grid.major = element_line(color = "grey80"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
    )

p <- plot_grid(p1,
               p2,
               ncol = 2,
               rel_widths = c(1, 1))

w <- 10
h <- 4.8
name <- args$file1

ggsave(p,
       file = paste0(name, ".jpg"),
       width = w,
       height = h)

