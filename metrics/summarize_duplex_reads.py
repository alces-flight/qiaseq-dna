from collections import defaultdict
import sys
import pysam
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

IN = pysam.AlignmentFile(sys.argv[1],"rb")

duplex_by_umi = defaultdict(lambda:defaultdict(int))

for read in IN:
    umi =  read.get_tag("Mi")
    duplex_id = read.qname.split(":")[-1]
    if duplex_id == "NN":
        continue

    duplex_by_umi[umi][duplex_id]+=1

duplex_ratio = []
dropped = 0
umi_plus_duplex=0
for umi in duplex_by_umi:    
    if duplex_by_umi[umi]['CC'] < 3:
        duplex_by_umi[umi]['CC'] = 0
    if duplex_by_umi[umi]['TT'] < 3:
        duplex_by_umi[umi]['TT'] = 0
    if duplex_by_umi[umi]['CC'] == 0 and duplex_by_umi[umi]['TT'] == 0:
        dropped+=1
        continue            

    if duplex_by_umi[umi]['CC'] >= 3:
        umi_plus_duplex+=1
    if duplex_by_umi[umi]['TT'] >= 3:
        umi_plus_duplex+=1
        
    try:        
        duplex_ratio.append(float(duplex_by_umi[umi]['CC'])/(duplex_by_umi[umi]['CC'] + duplex_by_umi[umi]['TT']))
    except Exception as e:
        print umi, duplex_by_umi[umi]
        raise(e)


duplex_ratio = sorted(duplex_ratio)
print "UMIs with all CC\tUMIs with no CC\tUMIs with atleast 3 CC or TT reads\tUMIs with 0-2 CC and TT reads, dropped\tNo. of unique UMIs(after including duplex tag)"
print str(duplex_ratio.count(1.0))+"\t"+str(duplex_ratio.count(0.0))+"\t"+str(len(duplex_ratio))+"\t"+str(dropped)+"\t"+str(umi_plus_duplex)

plt.plot(range(len(duplex_ratio)),duplex_ratio)
plt.xlabel('UMIs')
plt.ylabel('Fraction of Reads with CC')
plt.savefig('%s_duplex_ratio.png'%sys.argv[2])



IN.close()    
