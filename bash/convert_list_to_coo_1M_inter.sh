#! /bin/bash
sed 's/-/ /g' $1 | sed 's/M/000000/g' | awk '{if($5!=0.0) printf("%s\t%s\t%s\n",$2-1000000,$4-1000000,$5)}' > $2
