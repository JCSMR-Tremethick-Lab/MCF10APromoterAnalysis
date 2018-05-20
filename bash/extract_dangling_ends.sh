[#!/bin/bash

((ends=0))
for i in $(find . -name "*.RSstat")
  do sample=$(echo $i| cut -f 2 -d "/")
     sample_old=$sample
     ((ends+=$(grep -e "\<Valid_interaction_pairs\>" $i|cut -f 2)))
     if ($sample -ne $sample_old)
     then
       echo -e $sample_old"\t"$ends
     fi
done
]
