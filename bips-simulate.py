#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 30 14:47:49 2019

@author: dorgham
"""

import matplotlib.pyplot as plt
import networkx as nx
import random
import sys
import math


#constants definition
d=2
n=100000
k=2
simulations_per_degree = 100
degree_expInfecTime_map = dict()
plt.close()

y_ticks = [math.log2(math.log2(n)), math.log2(n), math.log2(n)**2, n]
y_ticks_labels = ['log log n', 'log n', '$log^2 n$', 'n']

#fix n and k. vary d to see its effect
while d <= round(n**(1/3)):
#the condition on the degree being O(n^(1/3)) is to guarantee the asymptotic uniform sampling from the space of random graphs,
# according to the documentation: https://networkx.github.io/documentation/networkx-1.10/reference/generated/networkx.generators.random_graphs.random_regular_graph.html
    infection_times_per_degree = []
    for s in range(simulations_per_degree):
        t=0
        infected_set = set()
        
        #************** Construct a random regular graph **************#
        n_trials=0
        print('\nGenerating regular graph of %d nodes and degree %d....' %(n,d))
        if d == 2:  #using nx.random_regular_graph() with d=2 does not generate connected graph most of the time
            G = nx.cycle_graph(n)
        else:
            while True:
                G = nx.random_regular_graph(d, n)
                if nx.is_connected(G):
                    break
                else:
                    n_trials += 1
                    print('Graph is not connected. Trying to generate another graph...')
                    if n_trials > 20:
                        sys.exit('Could not construct a connected regular graph with the given degree, node size!')
        print('Done generating a connected graph.')       
        
        all_nodes = list(G.nodes())
        #select a random node to be the persistent infected node
        persistent_node = random.choice(all_nodes)
        infected_set.add(persistent_node)
        
        
        #******************* Begining of the k-BIPS process *******************#
        print('%d-BIPS process is running now...'%k)
        while(len(infected_set) != n):
            prev_infec_set_size = len(infected_set)
            newly_infected = set()
            newly_uninfected = set()
            t += 1
            
            #each node samples k random neighbors
            for node in all_nodes:
                if node == persistent_node:
                    continue
                neighbors = list(G.neighbors(node))
                sampled_infected_neighbor = False
                for _ in range(k):
                    sampled_neighbor = random.choice(neighbors)
                    if sampled_neighbor in infected_set:
                        newly_infected.add(node)
                        sampled_infected_neighbor = True
                if not sampled_infected_neighbor and node in infected_set:
                    newly_uninfected.add(node)
                
            infected_set |= newly_infected
            for node in newly_uninfected:
                infected_set.remove(node)
            #print('End of round %d, Number of infected nodes is %d (Net Increase: %d).' %(t, len(infected_set), len(infected_set)-prev_infec_set_size))
            
        
        print('\nNumber of rounds until %d nodes are infected: %d' %(n, t))
        infection_times_per_degree.append(t)
        
    degree_expInfecTime_map[d] = round((sum(infection_times_per_degree) / len(infection_times_per_degree)))
    #increase d
    if d <= 5:
        d += 1
    elif d <= 10:
        d += 2
    elif d <= 20:
        d += 4
    else:
        d += 6
        
plt.plot(degree_expInfecTime_map.keys(), degree_expInfecTime_map.values(), 'b--')
plt.yscale('log')
plt.yticks(y_ticks, y_ticks_labels)
plt.xlabel('Node degree')
plt.ylabel('Expected infection time (in rounds)\nlogarithmic scale')
plt.title('k=%d, n=%d'%(k,n))
plt.grid(True)
plt.show()
plt.savefig('bips_variable_degree_k'+str(k)+'_n'+str(n)+'.png', bbox_inches='tight')
