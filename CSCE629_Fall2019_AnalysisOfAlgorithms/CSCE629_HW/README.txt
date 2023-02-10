The report is in 'CSCE_629_Analysis_of_ALgorithms'
The code is in 'main.py'
The results from python is in 'Results'


Code in 'main.py' can be run directly in 'python'
######################Library##################
import time
import random
import numpy as np
import math
#######################Function###################
VERT_SIZE = 5000 #size of vertices

class Edge:# input structure for Kruskal's method     

def ConstructNewGraph1(G1,E1): # Construct the spare graph
              
def ConstructNewGraph2(G2,E2):# Construct the dense graph

def DijkstraNoHeap(G,s,t): # no heap Dijkstra's algorithm
#0 unseen;1 fringe;2 intree
 

def heapifyUp(heap, value, fringe2heap, idx):# make array to be heap from low to up according to value         
            
def heapifyDown(heap, value, fringe2heap, idx):# make array to be heap from up to low according to value        
        
def heapInsert(heap, value, fringe2heap, vert):# insert a vertex(vert) into heap

def heapDelete(heap, value, fringe2heap, idx): # delete a vertex (idx) form heap

def DijkstraHeap( G, s, t): # Dijkstra's algorithm with heap
   
def kruskalsearch(T,capacity, s, t):# find bandwidth with construted T in Kruskal

#def MakeSet(vertex,dad,rank):
#    dad[vertex] = -1
#    rank[vertex] = 0
def Union(r1,r2,dad,rank):

def Find(vertex,dad):

def Kruskal(E, s, t): # first part: construct tree T; second part: kruskalseach
