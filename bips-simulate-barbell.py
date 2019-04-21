#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 21 15:14:37 2019

@author: dorgham
"""

###############################################################################
####        k-BIPS simulations for barbell-like regular graphs            #####
###############################################################################

import matplotlib.pyplot as plt
from pylab import MaxNLocator
from IPython import get_ipython
from collections import OrderedDict
import random
import math
import gc


ipy = get_ipython()
if ipy is not None:
    ipy.run_line_magic('matplotlib', 'inline')
plt.close()


n=1000
k=2
simulations_repetition = 10
max_increase = 0
longest_process_len = 0
fig_idx = 1


d_avgInfecTime_map = OrderedDict()
d_netIncrease_map = OrderedDict()
y_ticks = [math.log2(math.log2(n)), math.log2(n), math.log2(n)**2]
y_ticks_labels = ['$\log\ \log\ n$', '$\log\ n$', '$\log^2\ n$']

G_neighbors = dict()

while n < 10*n:
    infection_times_per_d = []
    net_increases_per_d = []
    
    print('\nGenerating barbell-like regular graph of %d nodes....' %(2*n))
    G1_nodes = [v for v in range(n)]
    G2_nodes = [v for v in range(n,2*n)]
    
    for node in G1_nodes:
        G_neighbors[node] = [v for v in G1_nodes if v != node]
    for node in G2_nodes:
        G_neighbors[node] = [v for v in G2_nodes if v != node]

    #remove_edge(0,1)
    G_neighbors[0].remove(1)
    G_neighbors[1].remove(0)
    #remove_edge(n,n+1)
    G_neighbors[n].remove(n+1)
    G_neighbors[n+1].remove(n)
    #add_edge(0,n)
    G_neighbors[0].append(n)
    G_neighbors[n].append(0)
    #add_edge(1,n+1)
    G_neighbors[1].append(n+1)
    G_neighbors[n+1].append(1)
    print('Done generating the graph.')  
    all_nodes = list(G_neighbors.keys())
    
    gc.collect()
    
    for s in range(simulations_repetition):
        t=0
        infected_set = set()
        #select a random node to be the persistent infected node
        persistent_node = random.choice(all_nodes)
        infected_set.add(persistent_node)
    
    
        #******************* Begining of the k-BIPS process *******************#
        print('%d-BIPS process is running now...'%k)
        net_increase_per_round = []
        while(len(infected_set) != n):
            prev_infec_set_size = len(infected_set)
            newly_infected = set()
            newly_uninfected = set()
            t += 1
            
            if prev_infec_set_size < n/3:
                nodes_to_loop_over = set(infected_set)
                for v in infected_set:
                    for u in G_neighbors[v]:
                        nodes_to_loop_over.add(u)
            else:
                nodes_to_loop_over = set(all_nodes)
            #each node samples k random neighbors
            for node in nodes_to_loop_over:
                if node == persistent_node:
                    continue
                
                neighbors = list(G_neighbors[node])
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
        
        print('Number of rounds until %d nodes are infected: %d' %(n, t))
        infection_times_per_d.append(t)
        net_increases_per_d.append(net_increase_per_round)
        if t > longest_process_len:
            longest_process_len = t
        
    d_avgInfecTime_map[2*n] = round((sum(infection_times_per_d) / len(infection_times_per_d)))
    d_netIncrease_map[2*n] = net_increases_per_d
            
    n += 1000
    
    
fig1 = plt.figure(fig_idx, figsize=(9,7.2))
plt.plot(d_avgInfecTime_map.keys(), d_avgInfecTime_map.values(), 'b--', label='average infection time')
#plt.plot(d_avgInfecTime_map.keys(), [math.log2(n)/ math.log2(i) for i in d_avgInfecTime_map.keys()], 'r-', label='$\log\ n\ /\ \log\ k$')
plt.yscale('log')
plt.yticks(y_ticks, y_ticks_labels)
xa = plt.gca().get_xaxis()
xa.set_major_locator(MaxNLocator(integer=True))
plt.xlabel('Number of nodes (n)')
plt.ylabel('Average infection time (in rounds)\nlogarithmic scale')
plt.title('%d-BIPS on Barbell-like regular graphs'%k)
plt.grid(True)
plt.legend(loc='best')
fig1.savefig('bips_barbell-like_variable_degree_k'+str(k)+'.png', bbox_inches='tight')
plt.close(fig1)

    