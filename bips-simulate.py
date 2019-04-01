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


#constants definition
d=90
n=20000
k=2

t=0
infected_set = set()
draw_infection_graph = False

#************** Construct a random regular graph **************#
n_trials=0
print('\nGenerating regular graph of %d nodes and degree %d....' %(n,d))
while True:
    G = nx.random_regular_graph(d, n)
    if nx.is_connected(G):
        break
    else:
        n_trials += 1
        print('Graph is not connected. Trying to generate another graph...')
        if n_trials > 20:
            sys.exit('Could not construct a connected regular graph with the given degree, node size!')
print('Done generating the graph.')       

all_nodes = list(G.nodes())
#select a random node to be the persistent infected node
persistent_node = random.choice(all_nodes)
infected_set.add(persistent_node)


#draw the initial infection graph  at round t=0
if draw_infection_graph:
    plt.close()
    fig, axes = plt.subplots(nrows=4, ncols=4)
    ax = axes.flatten()
    for a in range(len(ax)):
        ax[a].set_axis_off()
    pos = nx.spring_layout(G)
    color_map = []
    for node in G:
        if node == persistent_node:
            color_map.append('darkred')
        elif node in infected_set:
            color_map.append('red')
        else:
            color_map.append('lightgreen')
    nx.draw_networkx(G, ax=ax[0], node_color= color_map, font_color='white', pos=pos)
    ax[0].set_title('t=' + str(t))
        
    
#******************* Begining of the k-BIPS process *******************#
while(len(infected_set) != n):
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
    print('End of round %d, Number of infected nodes is %d' %(t, len(infected_set)))
    
    if draw_infection_graph:
        color_map = []
        for node in G:
            if node == persistent_node:
                color_map.append('darkred')
            elif node in infected_set:
                color_map.append('red')
            else:
                color_map.append('lightgreen')
        nx.draw_networkx(G, ax=ax[t], node_color= color_map, font_color='white', pos=pos)
        ax[t].set_title('t=' + str(t))
    

print('\n\nNumber of rounds until %d nodes are infected: %d' %(n, t))

if draw_infection_graph:
    plt.show()
