#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(sleuth)

so <- readRDS('data/sleuth_genes_MCF10Ca1aWT_vs_MCF10Ca1ashZ.rds')
sleuth::sleuth_live(so)