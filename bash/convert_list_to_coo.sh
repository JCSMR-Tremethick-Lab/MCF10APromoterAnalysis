#! /bin/bash
sed 's/-/ /g' $1 | sed 's/k/000/g' | awk '{if($5!=0.0) printf("chr%s\t%s\t%s\n",$2-40000,$4-40000,$5)}' > $2
