#!/bin/bash
conda activate deepTools
export bigwigDir="/home/sebastian/Data/Tremethick/Breast/PromoterSeqCap/processed_data/GRCh37_hg19_UCSC/deepTools/bamCoverage/MNase/RPKM"
export plotDir="/home/sebastian/Data/Tremethick/Breast/PromoterSeqCap/processed_data/GRCh37_hg19_UCSC/deepTools/plotHeatmap"
export matrixDir="/home/sebastian/Data/Tremethick/Breast/PromoterSeqCap/processed_data/GRCh37_hg19_UCSC/deepTools/computeMatrix"
export bedFile="/home/sebastian/Data/Tremethick/Breast/PromoterSeqCap/2000randompromotersHG19.bed"

cd ~/Data/Tremethick/Breast/PromoterSeqCap

computeMatrix reference-point -S $bigwigDir/MCF10A_WT_Input.Total.bw\
                                 $bigwigDir/MCF10A_TGFb_Input.Total.bw\
                                 $bigwigDir/MCF10A_shH2AZ_Input.Total.bw\
                                 $bigwigDir/MCF10CA1a_WT_Input.Total.bw\
                              -R $bedFile \
                              --outFileName $matrixDir/InputTotalMatrix.gz\
                              --outFileNameMatrix $matrixDir/InputTotalMatrix.tab\
                              --outFileSortedRegions $matrixDir/InputTotalMatrix.bed\
                              --referencePoint center\
                              --beforeRegionStartLength 1000\
                              --afterRegionStartLength 1000\
                              --binSize 10\
                              --smartLabels\
                              --numberOfProcessors 32\
&

plotHeatmap --matrixFile $matrixDir/InputTotalMatrix.gz\
            --outFileName $plotDir/kmeans/InputTotalMatrix.pdf\
            --outFileNameMatrix $plotDir/kmeans/InputTotalMatrix.tab\
            --outFileSortedRegions $plotDir/kmeans/InputTotalMatrix.bed\
            --kmeans 7\
&

plotHeatmap --matrixFile $matrixDir/InputTotalMatrix.gz\
            --outFileName $plotDir/kmeans/InputTotalMatrix_sorton1.pdf\
            --outFileNameMatrix $plotDir/kmeans/InputTotalMatrix_sorton1.tab\
            --outFileSortedRegions $plotDir/kmeans/InputTotalMatrix_sorton1.bed\
            --kmeans 7\
            --sortUsingSamples 1\
&

computeMatrixOperations subset --matrixFile $matrixDir/InputTotalMatrix.gz\
                               --samples MCF10A_WT_Input.Total\
                               --outFileName $matrixDir/InputTotalMatrix_WT.gz

plotHeatmap --matrixFile $matrixDir/InputTotalMatrix_WT.gz\
            --outFileName $plotDir/kmeans/InputTotalMatrix_WT.pdf\
            --outFileSortedRegions $plotDir/kmeans/InputTotalMatrix_WT.bed\
            --kmeans 7\
&

plotHeatmap --matrixFile $matrixDir/InputTotalMatrix.gz\
            --outFileName $plotDir/kmeans/InputTotalMatrix_sortedOnWT.pdf\
            
&                           


