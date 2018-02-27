#!/usr/bin/env python

import shutil
import os
import operator
import sys

# OXT LEU A 106
# cluster.pl  -kclust -radius 1 frame.pdb.* -log cluster.log -selmode heavy > cluster.output
# python ../../fetch_clusters_centroids.py
# for i in *.pdb;do python /home/ccorbi/Work/Beagle/SEED/add_chain_homodimers.py ${i} > ./chain/${i};done
# for i in *.pdb;do mkdir /home/ccorbi/Work/Beagle/SEED/AGAG/${i%.pdb};  python /home/ccorbi/Work/Beagle/SEED/add_chain_homodimers.py ${i} > /home/ccorbi/Work/Beagle/SEED/AGAG/${i%.pdb}/${i};done

clusters = dict()
filename = 'cluster.output'
cluster_pop = dict()

# parse cluster
with open(filename,'r') as input_file:
    for line in input_file:
        if line.startswith('@'):
            info = line.split()
            # init cluster
            clust_id = info[1]
            clusters[clust_id] = list()
            # population inside the cluster
            cluster_pop[clust_id] = int(info[3])



        try:
            data = line.split()
            if len(data) == 3:
                clusters[clust_id].append(data)
        except:
            pass

del clusters['t']
#print(clusters.keys())
del cluster_pop['t']



cluster_pop_sorted = sorted(cluster_pop.items(), key=operator.itemgetter(1), reverse=True)


if sys.argv[1]:
    limit = int(sys.argv[1])
else:
    limit = len(cluster_pop_sorted)

for idx, cluster_id in enumerate(cluster_pop_sorted):
    if idx < limit:
        cluster = cluster_id[0]
        data = clusters[cluster]
#for cluster, data in clusters.items():
        distance = 100
        centroid = ''
        for frame in data:

            if float(frame[2])< float(distance):
                distance = frame[2]
                centroid = frame[1]
        print('{} {} {}'.format(cluster,centroid,distance))

        os.makedirs('../clusters_representation', exist_ok=True)
        shutil.copyfile('../allframes/'+centroid, '../clusters_representation/'+cluster+'_centroid.pdb')
    else:
        pass

print('most populated cluster is {}'.format(cluster_pop_sorted[0]))
for i in range(1,len(cluster_pop_sorted)):
    if i > limit:
        print('#',cluster_pop_sorted[i])
    else:
        print(cluster_pop_sorted[i])
