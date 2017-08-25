import numpy as np

import mirnylib.genome
import mirnylib.h5dict
import sys
import math

#lib = mirnylib.h5dict.h5dict("/projects/imb-pkbphil/jp/hic/ASC_NEXT500_run3_PCollas_UOslo_HiC/Day-0-rep1_S3.hdf5")
lib = mirnylib.h5dict.h5dict(sys.argv[1])


N=lib['ids1'].shape[0]

mapq = np.full(N, 10, dtype=np.int)

chr2id = lib['misc']['genome']['label2idx']
id2chr = {v: k for k, v in chr2id.items()}
id2chr[-1] = '-1'

ab = np.zeros(N, dtype=[('var1','S41'),('var2',int),('var3','S2'),('var4',int),('var5',int),('var6',int),('var7','S2'),('var8',int),('var9',int), ('var10', int), ('var11',int) ])

chrs1= lib['chrms1']
chrs2= lib['chrms2']


c1 = np.array([id2chr[id] for id in chrs1], dtype="S2")
c2 = np.array([id2chr[id] for id in chrs2], dtype="S2")


filter = np.where(np.logical_and(lib['cuts1'] != -1, lib['cuts2'] != -1))

swap = chrs2 < chrs1

ab['var1'] = lib['ids1']
ab['var2'] = lib['strands1']
ab['var3'] = c1
ab['var4'] = lib['cuts1']
ab['var5'] = lib['rfragIdxs1']
ab['var6'] = lib['strands2']
ab['var7'] = c2
ab['var8'] = lib['cuts2']
ab['var9'] = lib['rfragIdxs2']
ab['var10'] = mapq
ab['var11'] = mapq

# Swap so that lowest chr numbers are to the left:

np.copyto(ab['var3'], c2, where=swap)
np.copyto(ab['var7'], c1, where=swap)

np.copyto(ab['var2'], lib['strands2'], where=swap)
np.copyto(ab['var6'], lib['strands1'], where=swap)

np.copyto(ab['var4'], lib['cuts2'], where=swap)
np.copyto(ab['var8'], lib['cuts1'], where=swap)

np.copyto(ab['var5'], lib['rfragIdxs2'], where=swap)
np.copyto(ab['var9'], lib['rfragIdxs1'], where=swap)

# sort on chromosomes:
#np.argsort(ab, order=['var3', 'var7'])

np.savetxt(sys.stdout, ab[filter], fmt="%s %i %s %i %i %i %s %i %i %i %i")
