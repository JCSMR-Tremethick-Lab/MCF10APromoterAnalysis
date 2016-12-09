STAR --runThreadN 64 \
     --runMode genomeGenerate \
     --genomeDir ./STAR_Index_ERCC/ \
     --genomeFastaFiles /home/sebastian/Data/References/Genomes/Homo_sapiens/GRCh38_ensembl84/primary/Homo_sapiens.GRCh38.dna.primary_assembly.fa ~/Data/References/Transcriptomes/ERCC/ERCC92.fa \
     --sjdbGTFfile /home/sebastian/Data/References/Annotations/Homo_sapiens/GRCh38_ensembl84/Homo_sapiens.GRCh38.84_ERCC.gtf \
     --sjdbOverhang 76


STAR --runThreadN 96 \
     --runMode genomeGenerate \
     --genomeDir ./STAR_Index_ERCC/ \
     --genomeFastaFiles /home/sebastian/Data/References/Genomes/Homo_sapiens/GRCh37_hg19_ensembl75/primary/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa ~/Data/References/Transcriptomes/ERCC/ERCC92.fa \
     --sjdbGTFfile /home/sebastian/Data/References/Annotations/Homo_sapiens/GRCh37_hg19_ensembl75/Homo_sapiens.GRCh37.75_ERCC.gtf \
     --sjdbOverhang 76

STAR --runThreadN 32 \
     --runMode genomeGenerate \
     --genomeDir ./STAR_Index/ \
     --genomeFastaFiles /home/sebastian/Data/References/Genomes/Homo_sapiens/GRCh37_hg19_ensembl75/primary/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
     --sjdbGTFfile /home/sebastian/Data/References/Annotations/Homo_sapiens/GRCh37_hg19_ensembl75/Homo_sapiens.GRCh37.75.gtf \
     --sjdbOverhang 76
