#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 30 14:47:49 2019

@author: dorgham
"""

import matplotlib.pyplot as plt
from pylab import MaxNLocator
from IPython import get_ipython
from collections import OrderedDict
import networkx as nx
import random
import sys
import math
import gc


vary_k = True

ipy = get_ipython()
if ipy is not None:
    ipy.run_line_magic('matplotlib', 'inline')
plt.close()

#constants definition
d=2
n=100000
k=2
simulations_repetition = 10
max_increase = 0
longest_process_len = 0
fig_idx = 1

###############################################################################
####                      fixing d, and varying k                         #####
###############################################################################

if vary_k:
    k_avgInfecTime_map = OrderedDict()
    k_netIncrease_map = OrderedDict()
    y_ticks = [math.log2(n), math.log2(n)**2, n]
    y_ticks_labels = ['$\log\ n$', '$\log^2\ n$', 'n']
    
    while d < round(n**(1/3)):     #repeat the simulations for different degrees
        k = 2
        while k <= 10*d:     #fix n and d. vary k to see its effect
            infection_times_per_k = []
            net_increases_per_k = []
            for s in range(simulations_repetition):
                t=0
                infected_set = set()
                
                #************** Construct a random regular graph **************#
                n_trials=0
                print('\nGenerating regular graph of %d nodes and degree %d....' %(n,d))
                if d == 2:  #using nx.random_regular_graph() with d=2 does not generate connected graph most of the time
                    G = nx.cycle_graph(n)
                elif d == n-1:
                    gc.collect()    #large graphs are memory consuming, so free the memory of the previous graph
                    G = nx.complete_graph(n)
                else:
                    while True:
                        G = nx.random_regular_graph(d, n)
                        if nx.is_connected(G):
                            break
                        else:
                            n_trials += 1
                            print('Graph is not connected. Trying to generate another graph...')
                            if n_trials > 10:
                                sys.exit('Could not construct a connected regular graph with the given degree, node size!')
                print('Done generating a connected graph.')       
                
                all_nodes = list(G.nodes())
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
                
                print('Number of rounds until %d nodes are infected: %d' %(n, t))
                infection_times_per_k.append(t)
                net_increases_per_k.append(net_increase_per_round)
                if t > longest_process_len:
                    longest_process_len = t
                
                if d==2:
                    break
                
            k_avgInfecTime_map[k] = round((sum(infection_times_per_k) / len(infection_times_per_k)))
            k_netIncrease_map[k] = net_increases_per_k
            #increase k
            if k < 5:
                k += 1
            elif k < 10:
                k += 2
            elif k < 20:
                k += 3
            elif k < 200:
                k += 5
            else:
                k += 20
                
        if d != 2:
            y_ticks = [math.log2(math.log2(n)), math.log2(n), math.log2(n)**2]
            y_ticks_labels = ['$\log\ \log\ n$', '$\log\ n$', '$\log^2\ n$']
            
        fig1 = plt.figure(fig_idx, figsize=(9,7.2))
        if d == 2:
            plt.plot(k_avgInfecTime_map.keys(), k_avgInfecTime_map.values(), 'b--', label='average infection time')
            plt.plot(k_avgInfecTime_map.keys(), [n/ math.log2(i) for i in k_avgInfecTime_map.keys()], 'r-', label='$n\ /\ \log\ k$')
        else:
           plt.plot(k_avgInfecTime_map.keys(), k_avgInfecTime_map.values(), 'b--', label='average infection time')
           plt.plot(k_avgInfecTime_map.keys(), [math.log2(n)/ math.log2(i) for i in k_avgInfecTime_map.keys()], 'r-', label='$\log\ n\ /\ \log\ k$')
        plt.yscale('log')
        plt.yticks(y_ticks, y_ticks_labels)
        xa = plt.gca().get_xaxis()
        xa.set_major_locator(MaxNLocator(integer=True))
        plt.xlabel('Branching factor (k)')
        plt.ylabel('Average infection time (in rounds)\nlogarithmic scale')
        plt.title('n=%d, d=%d'%(n,d))
        plt.grid(True)
        plt.legend(loc='best')
        fig1.savefig('bips_variable_branching_d'+str(d)+'_n'+str(n)+'.png', bbox_inches='tight')
        plt.close(fig1)
        
        
        fig_idx += 1
        subfig_idx = 1
        for br_factor in k_netIncrease_map.keys():
            fig = plt.figure(fig_idx, figsize=(20,16))
            ax = fig.add_subplot(2, 3, subfig_idx)
            xa = ax.get_xaxis()
            xa.set_major_locator(MaxNLocator(integer=True))
            plt.xlabel('Rounds until graph infection')
            plt.ylabel('Net increase of infected nodes per round')
            num_plots = 0
            for inc_lst in k_netIncrease_map[br_factor]:
                ax.plot(range(1, len(inc_lst)+1), inc_lst)
                if d == 2:
                    ax.set_ylim(-2.25, 2.25)
                    ax.set_xlim(0, longest_process_len)
                else:
                    ax.set_ylim(0, max_increase+1)
                    ax.set_xlim(0, longest_process_len)
                num_plots += 1
                if num_plots == 5:
                    break   #drwing all plots makes the plot not visible
            plt.title('k=%d'%br_factor)
            if subfig_idx==6:
                if d == 2:
                    fig.suptitle('1 simulation of k-BIPS on %d-regular graph with %d nodes' %(d,n))
                else:
                    fig.suptitle('5 simulations of k-BIPS on %d-regular graph with %d nodes' %(d,n))
                fig.savefig('bips_net_increase_per_round_d'+str(d)+'_n'+str(n)+'_figidx'+str(fig_idx-1)+'.png', bbox_inches='tight')
                fig_idx += 1
                subfig_idx = 1
                #break   #to plot the remaining degrees, comment this line 
            else:
                subfig_idx += 1
        #save the last figure
        if subfig_idx != 1:
            if d == 2:
                fig.suptitle('1 simulation of k-BIPS on %d-regular graph with %d nodes' %(d,n))
                fig.savefig('bips_net_increase_per_round_d'+str(d)+'_n'+str(n)+'_figidx'+str(fig_idx-1)+'.png', bbox_inches='tight')
            else:
                fig.suptitle('5 simulations of k-BIPS on %d-regular graph with %d nodes' %(d,n))
                fig.savefig('bips_net_increase_per_round_d'+str(d)+'_n'+str(n)+'_figidx'+str(fig_idx-1)+'.png', bbox_inches='tight')
            fig_idx += 1
        plt.close(fig)
        
                
        if d == n-1:
            break
        
        if d < 5:
            d += 1
        else:
            d += 5
        
        longest_process_len = 0
        max_increase = 0
    
###############################################################################
####                      fixing k, and varying d                         #####
###############################################################################
else:
    k_avgInfecTime_map = OrderedDict()
    k_netIncrease_map = OrderedDict()
    
    y_ticks = [math.log2(math.log2(n)), math.log2(n), math.log2(n)**2, n]
    y_ticks_labels = ['$\log\ \log\ n$', '$\log\ n$', '$\log^2\ n$', 'n']
    
    #fix n and k. vary d to see its effect
    while d <= round(n**(1/3)):
    #the condition on the degree being O(n^(1/3)) is to guarantee the asymptotic uniform sampling from the space of random graphs,
    # according to the documentation: https://networkx.github.io/documentation/networkx-1.10/reference/generated/networkx.generators.random_graphs.random_regular_graph.html
        infection_times_per_k = []
        net_increases_per_k = []
        for s in range(simulations_repetition):
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
            net_increase_per_round = []
            while(len(infected_set) != n):
                prev_infec_set_size = len(infected_set)
                newly_infected = set()
                newly_uninfected = set()
                t += 1
                
                if prev_infec_set_size < n/3:
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
            
            print('\nNumber of rounds until %d nodes are infected: %d' %(n, t))
            infection_times_per_k.append(t)
            net_increases_per_k.append(net_increase_per_round)
            if t > longest_process_len:
                longest_process_len = t
                
            if d==2:
                break
            
        k_avgInfecTime_map[d] = round((sum(infection_times_per_k) / len(infection_times_per_k)))
        k_netIncrease_map[d] = net_increases_per_k
        #increase d
        if d <= 5:
            d += 1
        elif d <= 10:
            d += 2
        elif d <= 20:
            d += 4
        else:
            d += 10
            
            
    fig1 = plt.figure(1, figsize=(9,7.2))
    plt.plot(k_avgInfecTime_map.keys(), k_avgInfecTime_map.values(), 'b--')
    plt.yscale('log')
    plt.yticks(y_ticks, y_ticks_labels)
    xa = plt.gca().get_xaxis()
    xa.set_major_locator(MaxNLocator(integer=True))
    plt.xlabel('Node degree')
    plt.ylabel('Average infection time (in rounds)\nlogarithmic scale')
    plt.title('k=%d, n=%d'%(k,n))
    plt.grid(True)
    fig1.show()
    fig1.savefig('bips_variable_degree_k'+str(k)+'_n'+str(n)+'.png', bbox_inches='tight')
    
    fig_idx = 2
    subfig_idx = 1
    for deg in k_netIncrease_map.keys():
        fig = plt.figure(fig_idx, figsize=(18,14))
        ax = fig.add_subplot(2, 3, subfig_idx)
        xa = ax.get_xaxis()
        xa.set_major_locator(MaxNLocator(integer=True))
        plt.xlabel('Rounds until graph infection')
        plt.ylabel('Net increase of infected nodes per round')
        num_plots = 0
        for inc_lst in k_netIncrease_map[deg]:
            ax.plot(range(1, len(inc_lst)+1), inc_lst)
            ax.set_ylim(0, max_increase)
            ax.set_xlim(0, longest_process_len)
            num_plots += 1
            if num_plots == 5:
                break   #drwing all plots makes the plot not visible
        plt.title('degree=%d'%deg)
        if subfig_idx==6:
            fig_idx += 1
            subfig_idx = 1
            fig.savefig('bips_net_increase_per_round_k'+str(k)+'_n'+str(n)+'_figidx'+str(fig_idx-1)+'.png', bbox_inches='tight')
            #break   #to plot the remaining degrees, comment this line 
        else:
            subfig_idx += 1
    plt.show()
            

