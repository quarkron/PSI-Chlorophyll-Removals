# Load pertinent modules
from __future__ import division 
import numpy as np
import scipy
from scipy.constants import hbar, Planck, c, elementary_charge
import pandas as pd 
import os 
import heapq 
import copy 

def overlap_integral(Ei, Ej, v=240, S=160):
    '''
    Calculate the spectral overlap integral necessary for the transfer rate calculation. 
    Parameters:
    Ei = absorption peak of donor, in cm^-1.
    Ej = absorption peak of receptor, in cm^-1. 
    v = Gaussian width. 240 cm^-1 in room temperature
    S = Stokes' shift. 160 cm^-1 in room temperature 
    '''
    return (1/(2*v*np.sqrt(np.pi)))*np.exp((-S**2 - Ei**2 - Ej**2 +2*S*Ei + 2*Ei*Ej - 2*S*Ej)/(4*v**2))

def coupling(node1, node2):
    C=116000
    r_ij = node1.coord() - node2.coord()
    R_ij = np.linalg.norm(r_ij)
    di, dj = node1.tdipVec(), node2.tdipVec()
    return C*(np.inner(di, dj))/((R_ij)**3) - (3*np.inner(r_ij, di)*np.inner(r_ij, dj))/((R_ij)**5)
    
def transfer_rate(node1, node2): #transfer rate in ps^-1
    Ei, Ej = node1.siteEnergy(), node2.siteEnergy()
    return coupling(node1, node2)**2*overlap_integral(Ei, Ej)*(2*np.pi)**2*scipy.constants.c*100*1e-12

def weight_func(node1, node2): #test weight function
    return 1/(transfer_rate(node1, node2))

class Node: #For PSI
    def __init__(self, df, i): #df is a DataFrame of chlorophyll molecules, i is the row number
        self._node = df.iloc[i]
        # if self._node['ID'][-5] != '0': self.name = self._node['ID'][-1] + self._node['ID'][-4:-2]
        # if self._node['ID'][-5] == '0': self.name = 'EC-' + self._node['ID'][-1] + self._node['ID'][-4:-2]
        # if self._node['ID'] == 'CLA1237.A': self.name = 'B37'
        # if self._node['ID'] == 'CLA1402.A': self.name = 'K02'
        # if self._node['ID'] == 'CLA1301.F': self.name = 'J01'
        # if self._node['ID'] == 'CLA1801.A': self.name = 'PL01'
        self.name = self._node['ID']
        self.dist = np.inf 
        self.previous = None 
        
    def coord(self):
        x,y,z = self._node['Center0'], self._node['Center1'], self._node['Center2']  
        return np.array([x,y,z])
    def id(self):
        return self._node['ID']
    def chainId(self):
        return self._node['Chain ID'] 
    def siteEnergy(self):
        return self._node['Energy (cm^-1)']
    def tdipVec(self):
        x,y,z = self._node['Tdip_vec0'], self._node['Tdip_vec1'], self._node['Tdip_vec2']  
        return np.array([x,y,z])
    
    def setDist(self, value): #for running Dijkstra 
        self.dist = value 

    # The following Class methods are addtional layers of Node abstraction to allow the use of priority queues in the Dijsktra algorithm.
    def clearPathMemory(self):
        self.previous = None
        self.dist = np.inf
        
    def __hash__(self):
        # Use a hash of the node's unique identifier (e.g., self.name)
        return hash(self.name)
        
    def __str__(self):
        return self.name
    
    def __eq__(self, other):
        return (self.name == other.name) and (self.dist == other.dist)

    def __ne__(self, other):
        return not (self == other)

    def __lt__(self, other):
        return (self.name < other.name) and (self.dist < other.dist)

    def __gt__(self, other):
        return  (self.name > other.name) and (self.dist > other.dist)
    
    def __le__(self, other):
        return (self < other) or (self == other)
    def __ge__(self, other):
        return (self > other) or (self == other)

    
class Edge:
    def __init__(self, src, dest):
        """ Assumes src and dest are Node objects"""
        self.src, self.dest = src, dest
    def setWeight(self, wght):
        self.weight = wght
    def getSource(self):
        return self.src 
    def getDestination(self):
        return self.dest
    def getWeight(self):
        return self.weight
    def __str__(self):
        return self.src.getName() + '->(' + str(self.weight) + ')' + self.dest.getName()
    
class Graph:
    def __init__(self):
        # self.nodes are lists of Node objects
        # self.edges is a dictionary mapping each node to a list of neighbor nodes
        # self.weight is a dictionary with the edges as keys, and vals as weights.
        self.nodes, self.edges, self.weights = [], {}, {}
        self.k_cs = 1.5
        self.k_diss = 0.001
    
    def P700lookup(self):
        for i in self.nodes:
            if str(i) == 'CLA1011.A':
                self.p700a = i
            if str(i) == 'CLA1021.B':
                self.p700b = i

    def init_Weights(self):
        weights_Total = np.sum(list(self.weights.values()))
        self.Weights = {key: value/weights_Total for key,value in self.weights.items()} #Normalized edge weights for standard network efficiency calculation.
        
    def init_Costs(self):
        weights_Total = np.sum(list(self.weights.values()))
        self.Costs = {key: weights_Total/value for key,value in self.weights.items()} #Edge cost

    def addNode(self, node):
        if node in self.nodes:
            raise ValueError('Duplicate Node')
        else:
            self.nodes.append(node)
            self.edges[node] = []

    def addEdge(self, edge):
        src, dest = edge.getSource(), edge.getDestination()
        if not (src in self.nodes and dest in self.nodes):
            raise ValueError('One or more nodes not in graph')
        to_add = (src, dest)
        if to_add in self.weights:
            raise ValueError('Duplicate Edge')
        self.edges[src].append(dest)
        self.weights[to_add] = edge.getWeight()

    def removeNode(self, node):
        if node not in self.nodes:
            raise ValueError('Node not found in the graph')

        # Remove the node from the nodes list
        self.nodes.remove(node)

        # Remove the node from the edges dictionary and its corresponding edges
        if node in self.edges:
            del self.edges[node]

        # Remove edges involving the node from the weights dictionary
        edges_to_remove = [edge for edge in self.weights if node in edge]
        for edge in edges_to_remove:
            del self.weights[edge]

        # Remove references to the removed node from the neighboring nodes
        for neighbor, neighbors_list in self.edges.items():
            self.edges[neighbor] = [n for n in neighbors_list if n != node]

    def neighborsOf(self, node):
        return self.edges[node]
    
    def hasNode(self, node):
        return node in self.nodes 
    
    def allNodes(self):
        return self.nodes 
    
    def allEdges(self):
        return self.weights 

    def dijkstra(self, start_node): #Dijkstra, optimized for network efficiency calculation. Apply only for antenna nodes
        # Initialize distances to all nodes as infinity
        costs = {node: float('inf') for node in self.nodes}

        # The distance to the source node is set to 0
        costs[start_node] = 0

        # Initialize a priority queue (min-heap) with the source node
        pq = [(0, start_node)]

        # Dictionary to store the previous node on the shortest path
        previous = {}

        while pq:
            current_dist, current_node = heapq.heappop(pq)

            # If you want to store the previous node on the path
            if current_node not in previous:
                previous[current_node] = None

            # Visit neighbors
            for neighbor in self.edges[current_node]:
                # Calculate the tentative distance to this neighbor
                tentative_dist = current_dist + self.Costs[(current_node, neighbor)]

                # If this path is shorter than the recorded distance, update it
                if tentative_dist < costs[neighbor]:
                    costs[neighbor] = tentative_dist
                    previous[neighbor] = current_node
                    heapq.heappush(pq, (tentative_dist, neighbor))

        return costs, previous
    
    def getShortestPath(self, previous, destination): #optimized for Neff calculation
        path = []
        while destination is not None:
            path.insert(0, destination)
            destination = previous[destination]
        return path
    
    def shortestPath(self, start_node, dest_node): #Dijkstra algorithm, optimized for BC calculation
        # dest_node is a list of destination Node objects
        # start_node a Node object
        weights_Total = np.sum(list(self.weights.values()))
        for i in self.nodes:
            if i != start_node:
                i.setDist(np.inf)
        pq = [(0, start_node)] 
        visited = set()
        shortest_path = []
        flag = [dest_node]
        i = 0

        while pq:
            # print('pq1:', [(i[0], str(i[1])) for i in pq])
            current_tuple = heapq.heappop(pq)
            # print('Current node:', (current_tuple[0], str(current_tuple[1]))) 
            
            if current_tuple[1] in visited:
                # print('vertex skipped')
                continue

            visited.add(current_tuple[1])

            if current_tuple[1] == dest_node:
                flag.remove(dest_node)
                break
            
            # if current_tuple[1] not in self.edges:
            #     print('Problem, Houston!')
                
            for adj_node in self.edges[current_tuple[1]]:
                
                
                if adj_node in visited:
                    # print('Node already visited')
                    continue

                # print(str(current_tuple[0]), str(current_tuple[1]), str(adj_node)) 
                tentative_dist = current_tuple[0] + self.Costs[(current_tuple[1], adj_node)]  #Reciprocal weight since the algorithm finds the min path, not max path. Old naive method
                
                if  tentative_dist < adj_node.dist:
                    adj_node.setDist(tentative_dist)
                    # print('Adjacent node:', (tentative_dist, adj_node.dist, str(adj_node)))
                    # print('pqa:', [(i[0], str(i[1])) for i in pq])
                    # print('To push: ', (adj_node.dist, str(adj_node)))
                    adj_node.previous = current_tuple[1] 
                    heapq.heappush(pq, (adj_node.dist, adj_node))

                i += 1
            # print('pq2:', [(i[0], str(i[1])) for i in pq])
            # print('visited nodes:', [str(i) for i in visited], '\n')
        
        if flag: #Return NoneType if no shortest paths exist, occurs if the nodes are disconnected.
            return None

        node = dest_node
        shortest_path=[node]
        while node.previous:
            shortest_path.insert(0, node.previous)
            node = node.previous

        for node in self.nodes: #Clear algorithm memory for all nodes
            node.clearPathMemory()
        
        return shortest_path
    
    def transferToP700(self, node): 
        if str(node) not in ['CLA1011.A', 'CLA1021.B']:
            if self.weights[node, self.p700a] < self.weights[node, self.p700b]:
                return self.weights[node, self.p700b] 
            else:
                return self.weights[node, self.p700a]
        else:
            return 0
    
    def euclidean_node_P700distance(self, node1):
        coordp700a, coordp700b = self.p700a.coord(), self.p700b.coord()
        dist1, dist2 = np.linalg.norm(coordp700a - node1.coord()), np.linalg.norm(coordp700b - node1.coord())
        return np.min([dist1, dist2])
        
    def distance(self, node1, node2):
        SP = self.shortestPath(node1, node2) 
        dist = 0
        if SP:
            for i in range(len(SP)-1):
                dist += self.Weights[(SP[i], SP[i+1])] 
        else:
            dist = np.inf
        return dist 
    
    def costdistance(self, node1, node2):
        SP = self.shortestPath(node1, node2) 
        cost = 0
        if SP:
            for i in range(len(SP)-1):
                cost += self.Costs[(SP[i], SP[i+1])] 
        else:
            cost = np.inf
        return cost 
    
    def meanDegrees(self):
        degrees = 2*sum([len(i) for i in self.edges.values()])/len(self.nodes)
        return degrees
    
    def meanInDegrees(self):
        return sum([len(i) for i in self.edges.values()])/len(self.nodes)
    
    def meanOutDegrees(self):
        return sum([len(i) for i in self.edges.values()])/len(self.nodes)
    
    def meanStrength(self):
        return sum(self.weights.values())/len(self.nodes)
    
    def meanOutStrength(self):
        return self.meanStrength()/2 
    
    def meanInStrength(self):
        return self.meanStrength()/2 

    def betweenessCentrality(self): 
        bc_dict = dict.fromkeys(self.nodes, 0)
        for i in range(len(self.nodes)):
            if str(self.nodes[i]) == 'CLA1011.A':
                ecA1 = self.nodes[i]
            if str(self.nodes[i]) == 'CLA1021.B':
                ecB1 = self.nodes[i]
        for i in self.nodes:
            if str(i) not in ['CLA1011.A', 'CLA1021.B']:
                bc_dict[i] = 0
                sp1=self.shortestPath(i, ecA1) #Compare shortest paths to ec-A1 and ec-B1
                sp2= self.shortestPath(i, ecB1)
                cost1, cost2 = 0, 0
                if (sp1 == None) and (sp2 != None):
                    sp = sp2
                if (sp1 != None) and (sp2 == None):
                    sp = sp1
                else:
                    for i in range(len(sp1)-1):
                        cost1 += self.Costs[(sp1[i], sp1[i+1])]
                    for j in range(len(sp2) -1):
                        cost2 += self.Costs[(sp2[j], sp2[j+1])]
                    if cost1 < cost2: # select the path with the lower cost to either.
                        sp = sp1
                    if cost2 < cost1:
                        sp = sp2 
                for node in sp[1:-1]: #do not include the start and the dest nodes.
                    bc_dict[node] += 1
        return bc_dict

    def nodeOutStrength(self, node):
        runningTotal = 0
        for neighborNode in self.edges[node]: runningTotal += self.Weights[(node, neighborNode)] 
        return runningTotal 
    
    def nodeInStrength(self, node):
        runningTotal = 0 
        for (src,dest) in list(self.weights.keys()):
            if dest == node: runningTotal += self.Weights[(src, dest)] 
        return runningTotal

    def K_ij(self): #Quantum master equation, K_ij units in ps^-1
        K = np.eye(len(self.nodes))
        for i in range(len(self.nodes)):
            for j in range(len(self.nodes)): 
                if i != j:
                    K[i,j] = self.weights[(self.nodes[j], self.nodes[i])]
                if i == j:
                    T_ik = 0 
                    for k in range(len(self.nodes)):
                        if i != k:
                            T_ik += self.weights[(self.nodes[i], self.nodes[k])]
                    if str(self.nodes[i]) in ['CLA1011.A', 'CLA1021.B']:
                        K[i,j] = -T_ik - self.k_diss - self.k_cs 
                    else:
                        K[i,j] = -T_ik - self.k_diss
        self.K = K

    def T_ij(self): #Quantum master equation, T_ij units in ps^-1
        T = np.eye(len(self.nodes))
        for i in range(len(self.nodes)):
            for j in range(len(self.nodes)): 
                if i != j:
                    T[i,j] = self.weights[(self.nodes[i], self.nodes[j])]
                if i == j:
                    T[i,j] = 0 
        self.T = T
    
    def Qeff(self):
        ones = np.ones(len(self.nodes))
        ones_t = ones.reshape(-1,1) 
        K_inv = np.linalg.inv(self.K)
        return 1-(-np.matmul(np.matmul(ones, K_inv), ones_t)*self.k_diss/len(self.nodes)) 
    
    def Neff(self): #standard network efficiency metric, measures only the antenna efficiency
        antenna_indices=[]
        N = len(self.nodes)
        distances = np.full((N, N), np.inf)

        for i in range(len(self.nodes)):
            if str(self.nodes[i]) not in ['CLA1011.A', 'CLA1021.B']:
                antenna_indices.append(i)
        
        M = len(antenna_indices)
        for i in antenna_indices:
            dijkstra_run = self.dijkstra(self.nodes[i])
            for j in antenna_indices:
                if i == j:
                    continue
                distances[i][j] = dijkstra_run[0][self.nodes[j]]
        
        # Calculate the sum of inverses of distances
        inverse_distances_sum = np.sum(1 / distances) 

        if np.isnan(inverse_distances_sum/(M*(M-1))):
            return 0
        else:
            return inverse_distances_sum/(M*(M-1))
    
def interferenceBC(g, node, to_remove):
    G = copy.deepcopy(g)
    #calculate BC interference of node "node" with respect to removing node "to_remove". node and to_remove are Node objects
    bc_dict = G.betweenessCentrality()
    unperturbed_BC = bc_dict[node]
    G.removeNode(to_remove)
    bc_dict_perturbed = G.betweenessCentrality()
    perturbed_BC = bc_dict_perturbed[node]
    return unperturbed_BC - perturbed_BC

def robustnessBC(g, node):
    interferences = []
    for i in g.nodes:
        if i != node and str(i) not in ['CLA1011.A', 'CLA1021.B']:
            interferences.append(interferenceBC(g, node, i))
    return 1/np.abs(max(interferences))