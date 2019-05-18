#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 21 15:14:37 2019

@author: dorgham
"""

###############################################################################
####            k-BIPS simulations for toroidal graphs                    #####
###############################################################################

import matplotlib.pyplot as plt
from pylab import MaxNLocator
from IPython import get_ipython
from collections import OrderedDict
import networkx as nx
import random
import math


ipy = get_ipython()
if ipy is not None:
    ipy.run_line_magic('matplotlib', 'inline')
plt.close()


dim = 4
N=100
k=2
simulations_repetition = 10
max_increase = 0
longest_process_len = 0


k_avgInfecTime_map = OrderedDict()
k_netIncrease_map = OrderedDict()
y_ticks = [math.log2(N)**2, N, dim*N]
y_ticks_labels = ['$\log^2\ N$', 'N', '$d\ N$']


print('\nGenerating [0,N]^d toroidal graph with N=%d and d=%d....' %(N,dim))

G=nx.grid_graph(dim=[dim-1, N, N], periodic=True)

print('Done generating the graph.')  
all_nodes = list(G.nodes)


while k <= 8:
    infection_times_per_k = []
    net_increases_per_k = []
    
    for s in range(simulations_repetition):
        t=0
        infected_set = set()
        #select a random node to be the persistent infected node
        persistent_node = random.choice(all_nodes)
        infected_set.add(persistent_node)
    
    
        #******************* Begining of the k-BIPS process *******************#
        print('%d-BIPS process is running now...'%k)
        net_increase_per_round = []
        while(len(infected_set) != len(all_nodes)):
            prev_infec_set_size = len(infected_set)
            newly_infected = set()
            newly_uninfected = set()
            t += 1
            
            if prev_infec_set_size < len(all_nodes)/3:
                nodes_to_loop_over = set(infected_set)
                for v in infected_set:
                    for u in G.neighbors(v):
                        nodes_to_loop_over.add(u)
            else:
                nodes_to_loop_over = set(all_nodes)
            #each node samples k random neighbors
            for node in nodes_to_loop_over:
                if node == persistent_node:
                    continue
                
                neighbors = list(G.neighbors(node))
                has_infected_neighbor = False
                if node not in infected_set:
                    for v in neighbors:
                        if v in infected_set:
                            has_infected_neighbor = True
                            break
                    if not has_infected_neighbor:
                        continue    #no need to sample neighbors of this node
                sampled_infected_neighbor = False
                for _ in range(k):
                    sampled_neighbor = random.choice(neighbors)
                    if sampled_neighbor in infected_set:
                        newly_infected.add(node)
                        sampled_infected_neighbor = True
                        break
                if not sampled_infected_neighbor and node in infected_set:
                    newly_uninfected.add(node)
                
            infected_set |= newly_infected
            for node in newly_uninfected:
                infected_set.remove(node)
            #print('End of round %d, Number of infected nodes is %d (Net Increase: %d).' %(t, len(infected_set), len(infected_set)-prev_infec_set_size))
            net_increase_per_round.append(len(infected_set)-prev_infec_set_size)
            if len(infected_set)-prev_infec_set_size > max_increase:
                max_increase = len(infected_set)-prev_infec_set_size
        
        print('Number of rounds until %d nodes are infected: %d\n' %(len(all_nodes), t))
        infection_times_per_k.append(t)
        net_increases_per_k.append(net_increase_per_round)
        if t > longest_process_len:
            longest_process_len = t
        
    k_avgInfecTime_map[k] = round((sum(infection_times_per_k) / len(infection_times_per_k)))
    k_netIncrease_map[k] = net_increases_per_k
            
    k += 1
    
    
fig1 = plt.figure(1, figsize=(9,7.2))
plt.plot(k_avgInfecTime_map.keys(), k_avgInfecTime_map.values(), 'b--', label='average infection time')
#plt.plot(k_avgInfecTime_map.keys(), [math.log2(n)**2/ math.log2(i) for i in k_avgInfecTime_map.keys()], 'g-', label='$\log^2\ n\ /\ \log\ k$')
#plt.yscale('log')
plt.yticks(y_ticks, y_ticks_labels)
xa = plt.gca().get_xaxis()
xa.set_major_locator(MaxNLocator(integer=True))
plt.xlabel('Branching factor (k)')
plt.ylabel('Average infection time (in rounds)')
plt.title('k-BIPS on $[0,N]^d$ toroidal graph with N=%d and d=%d'%(N, dim))
plt.grid(True)
#plt.legend(loc='best')
fig1.savefig('bips_toroidal_d'+str(dim)+'_N'+str(N)+'.png', bbox_inches='tight')
plt.close(fig1)

    