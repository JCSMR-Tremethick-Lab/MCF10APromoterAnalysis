#Remove up and down 1Mb from unmappable
#Add one more argument for hg19
#Change resolution variable based on your input resolution


import pybedtools
import sys
import re
import interval as iv

raw_fname=str(sys.argv[1])
unmap=str(sys.argv[2])
hg19_file=str(sys.argv[3])
first_chr=str(sys.argv[4])
sec_chr=str(sys.argv[5])

resolution=1000000

with open(unmap) as unf:
  unmap_file = unf.readlines()

with open(raw_fname) as rf:
  raw_file=rf.readlines()

raw_file=[x.strip() for x in raw_file]
unmap_file=[x.strip() for x in unmap_file]

unmap_dict = {}
for l in unmap_file:
  l = l.split()
  unmap_start = int(l[1]) - resolution
  if unmap_start < 0:
    unmap_start = int(l[1])

  unmap_end = int(l[2]) + resolution

  #unmap_dict.setdefault(l[0],[]).append(iv.interval[int(l[1]),int(l[2])])
  unmap_dict.setdefault(l[0],[]).append(iv.interval[unmap_start,unmap_end])


#hg19=pybedtools.helpers.get_chromsizes_from_ucsc("hg19")

f=open(hg19_file)
hg19={}
for l in f:
  hg19[l.split()[0]] = int(float(l.split()[1]))

first_chr_size = hg19[first_chr] #[1]
sec_chr_size = hg19[sec_chr] #[1]

f.close()

bedpe = []
for line in raw_file:
  line=line.split()
  end1 = int(line[0]) + resolution
  end2 = int(line[1]) + resolution

  chrA_interval = iv.interval[line[0],end1]
  chrB_interval = iv.interval[line[1],end2]

  flag = 1
  for val in unmap_dict[first_chr]:
    if(val & chrA_interval):
      flag = 0

  if(flag != 0):
    for val in unmap_dict[sec_chr]:
      if(val & chrB_interval):
        flag = 0

  if(flag == 1):
    bedpe.append([first_chr,line[0],end1,sec_chr,line[1],end2,int(float(line[2]))])

for n in bedpe:
  if n[2] > first_chr_size:
    n[2] = first_chr_size
  if n[5] > sec_chr_size:
    n[5] = sec_chr_size
  n=map(str,n)
  print '\t'.join(n)
