#!/bin/bash
source activate seqAnalysis
cd ~/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/

# 95% ID
cd-hit-est -i TOTALcombined_shH2AZ_Inp.fa -o TOTALcombined_shH2AZ_Inp.cluster95 -M 0 -n 8 -c 0.95 -d 0 -T 16
cd-hit-est -i TOTALcombined_CA1a_Inp.fa -o TOTALcombined_CA1a_Inp.cluster95 -M 0 -n 8 -c 0.95 -d 0 -T 16
cd-hit-est -i TOTALcombined_A_H2AZ.fa -o TOTALcombined_A_H2AZ.cluster95 -M 0 -n 8 -c 0.95 -d 0 -T 16
cd-hit-est -i TOTALcombined_A_Inp.fa -o TOTALcombined_A_Inp.cluster95 -M 0 -n 8 -c 0.95 -d 0 -T 16
cd-hit-est -i TOTALcombined_A_TGFb.fa -o TOTALcombined_A_TGFb.cluster95 -M 0 -n 8 -c 0.95 -d 0 -T 16
cd-hit-est -i TOTALcombined_CA1a_H2AZ.fa -o TOTALcombined_CA1a_H2AZ.cluster95 -M 0 -n 8 -c 0.95 -d 0 -T 16

# 90% ID
cd-hit-est -i TOTALcombined_shH2AZ_Inp.fa -o TOTALcombined_shH2AZ_Inp.cluster90 -M 0 -n 8 -c 0.90 -d 0 -T 32

cd-hit-est -i peakCallingTest.fa -o peakCallingTest.cluster95 -M 0 -n 8 -c 0.95 -d 0 -T 16
