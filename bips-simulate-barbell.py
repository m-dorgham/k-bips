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
import networkx as nx
import random
import math
import gc


ipy = get_ipython()
if ipy is not None:
    ipy.run_line_magic('matplotlib', 'inline')
plt.close()

vary_d = True

if vary_d:
    n = 10000
    k=2
    d=3
    simulations_repetition = 20
    max_increase = 0
    longest_process_len = 0
    fig_idx = 1
    
    
    d_avgInfecTime_map = OrderedDict()
    d_netIncrease_map = OrderedDict()
    y_ticks = [math.log2(math.log2(n)), math.log2(n), math.log2(n)**2, n]
    y_ticks_labels = ['$\log\ \log\ n$', '$\log\ n$', '$\log^2\ n$', 'n']
    
    G_neighbors = dict()
    nodes_relabels = {}
    block_size = int(n/2)
    
    while d < block_size:
        infection_times_per_d = []
        net_increases_per_d = []
        
        print('\nGenerating barbell-like regular graph of %d nodes and degree %d....' %(n,d))
        
        if d == block_size-1:
            G1_nodes = [v for v in range(block_size)]
            G2_nodes = [v for v in range(block_size, n)]
            for node in G1_nodes:
                G_neighbors[node] = [v for v in G1_nodes if v != node]
            for node in G2_nodes:
                G_neighbors[node] = [v for v in G2_nodes if v != node]
        else:
            G1 = nx.random_regular_graph(d, block_size)
            for node in G1:
                G_neighbors[node] = list(G1.neighbors(node))
            del G1
            gc.collect()
            
            G2 = nx.random_regular_graph(d, block_size)
            for i in range(block_size):
                nodes_relabels[i] = i + block_size
            nx.relabel_nodes(G2, nodes_relabels, copy=False)
            for node in G2:
                G_neighbors[node] = list(G2.neighbors(node))
            del G2
            gc.collect()

    
        #remove_edge(0,1)
        node0_neighbor = G_neighbors[0][0]
        G_neighbors[0].remove(node0_neighbor)
        G_neighbors[node0_neighbor].remove(0)
        #remove_edge(n,n+1)
        node1_neighbor = G_neighbors[int(n/2)][0]
        G_neighbors[int(n/2)].remove(node1_neighbor)
        G_neighbors[node1_neighbor].remove(int(n/2))
        #add_edge(0,n)
        G_neighbors[0].append(int(n/2))
        G_neighbors[int(n/2)].append(0)
        #add_edge(1,n+1)
        G_neighbors[node0_neighbor].append(node1_neighbor)
        G_neighbors[node1_neighbor].append(node0_neighbor)
        print('Done generating the graph.')  
        
        all_nodes = list(G_neighbors.keys())
        
        
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
            
            print('Number of rounds until %d nodes are infected: %d\n' %(n, t))
            infection_times_per_d.append(t)
            net_increases_per_d.append(net_increase_per_round)
            if t > longest_process_len:
                longest_process_len = t
            
        d_avgInfecTime_map[d] = round((sum(infection_times_per_d) / len(infection_times_per_d)))
        d_netIncrease_map[d] = net_increases_per_d
        
        if d < 0.1*block_size:
            d += 4
        elif d < 0.3*block_size:
            d += 10
        elif d < 0.6*block_size:
            d += 20
        elif d < block_size-1:
            d = block_size-1
        else:
            d = block_size
        
    
    fig1 = plt.figure(fig_idx, figsize=(9,7.2))
    plt.plot(d_avgInfecTime_map.keys(), d_avgInfecTime_map.values(), 'b--', label='average infection time')
    #plt.plot(d_avgInfecTime_map.keys(), [math.log2(n)**2/ math.log2(i)**2 for i in d_avgInfecTime_map.keys()], 'r-', label='$\log^2\ n\ /\ \log^2\ k$')
    plt.yscale('log')
    plt.yticks(y_ticks, y_ticks_labels)
    xa = plt.gca().get_xaxis()
    xa.set_major_locator(MaxNLocator(integer=True))
    plt.xlabel('Degree (d)')
    plt.xlim(3, block_size)
    plt.ylabel('Average infection time (in rounds)\nlogarithmic scale')
    plt.title('%d-BIPS on Barbell-like regular graphs'%k)
    plt.grid(True)
    #plt.legend(loc='best')
    fig1.savefig('bips_barbell-like_variable-d_n'+str(n)+'_k'+str(k)+'.png', bbox_inches='tight')
    plt.close(fig1)
    
    
else:
    n = 10000
    k=2
    simulations_repetition = 10
    max_increase = 0
    longest_process_len = 0
    fig_idx = 1
    
    
    k_avgInfecTime_map = OrderedDict()
    k_netIncrease_map = OrderedDict()
    y_ticks = [math.log2(math.log2(n)), math.log2(n), math.log2(n)**2, n]
    y_ticks_labels = ['$\log\ \log\ n$', '$\log\ n$', '$\log^2\ n$', 'n']
    
    G_neighbors = dict()
    
    print('\nGenerating barbell-like regular graph of %d nodes....' %(n))
    G1_nodes = [v for v in range(int(n/2))]
    G2_nodes = [v for v in range(int(n/2), n)]
    
    for node in G1_nodes:
        G_neighbors[node] = [v for v in G1_nodes if v != node]
    for node in G2_nodes:
        G_neighbors[node] = [v for v in G2_nodes if v != node]
    
    #remove_edge(0,1)
    G_neighbors[0].remove(1)
    G_neighbors[1].remove(0)
    #remove_edge(n,n+1)
    G_neighbors[int(n/2)].remove(int(n/2)+1)
    G_neighbors[int(n/2)+1].remove(int(n/2))
    #add_edge(0,n)
    G_neighbors[0].append(int(n/2))
    G_neighbors[int(n/2)].append(0)
    #add_edge(1,n+1)
    G_neighbors[1].append(int(n/2)+1)
    G_neighbors[int(n/2)+1].append(1)
    print('Done generating the graph.')  
    all_nodes = list(G_neighbors.keys())
    
    del G1_nodes, G2_nodes
    gc.collect()
    
    while k <= n:
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
            
            print('Number of rounds until %d nodes are infected: %d\n' %(n, t))
            infection_times_per_k.append(t)
            net_increases_per_k.append(net_increase_per_round)
            if t > longest_process_len:
                longest_process_len = t
            
        k_avgInfecTime_map[k] = round((sum(infection_times_per_k) / len(infection_times_per_k)))
        k_netIncrease_map[k] = net_increases_per_k
                
        k *= 2
        
        
    fig1 = plt.figure(fig_idx, figsize=(9,7.2))
    plt.plot(k_avgInfecTime_map.keys(), k_avgInfecTime_map.values(), 'b--', label='average infection time')
    plt.plot(k_avgInfecTime_map.keys(), [math.log2(n)**2/ math.log2(i)**2 for i in k_avgInfecTime_map.keys()], 'r-', label='$\log^2\ n\ /\ \log^2\ k$')
    plt.yscale('log')
    plt.yticks(y_ticks, y_ticks_labels)
    xa = plt.gca().get_xaxis()
    xa.set_major_locator(MaxNLocator(integer=True))
    plt.xlabel('Branching factor (k)')
    plt.ylabel('Average infection time (in rounds)\nlogarithmic scale')
    plt.title('k-BIPS on Barbell-like regular graphs')
    plt.grid(True)
    plt.legend(loc='best')
    fig1.savefig('bips_barbell-like_fixed-d_n'+str(n)+'.png', bbox_inches='tight')
    plt.close(fig1)

    