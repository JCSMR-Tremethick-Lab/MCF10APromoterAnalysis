qsub -pe threads 16 \
     -q hugemem.q \
     -l virtual_free=96 \
      -l h_vmem=96.1 \
      -e /home/skurscheid/qsub_output \
      -o /home/skurscheid/qsub_error qsub.sh
