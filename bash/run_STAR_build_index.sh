STAR --runThreadN 64 \
     --runMode genomeGenerate \
     --genomeDir ./STAR_Index_ERCC/ \
     --genomeFastaFiles /home/skurscheid/Data/References/Genomes/Homo_sapiens/GRCh38_ensembl84/primary/Homo_sapiens.GRCh38.dna.primary_assembly.fa ~/Data/References/Transcriptomes/ERCC/ERCC92.fa \
     --sjdbGTFfile /home/skurscheid/Data/References/Annotations/Homo_sapiens/GRCh38_ensembl84/Homo_sapiens.GRCh38.84_ERCC.gtf \
     --sjdbOverhang 76
