# Hello World program in Python
#reset 
import time
import random
import numpy as np
import math

VERT_SIZE = 5000

class VertDeg:
    def __init__(self,index,degree):
        self.index = index
        self.degree = degree
    def __repr__(self):
        return repr((self.index, self.degree))

class Edge:
    def __init__(self,vertex1,vertex2,weight):
        self.vertex1 = vertex1
        self.vertex2 = vertex2
        self.weight = weight
    def __repr__(self):
        return repr((self.vertex1, self.vertex2, self.weight))        



def ConstructNewGraph1(G1,E1):
    E1.clear()
    for i in range(0,VERT_SIZE-2):
        G1[i,i] = 0
        wt = random.randint(1,100000)
        e = Edge(i,i+1,wt)
        E1.append(e)
        #e = Edge(i+1,i,wt)
        #E1.append(e)
        G1[i,i+1] = wt
        G1[i+1,i] = wt
        
    NewSize = int((VERT_SIZE-2)*(VERT_SIZE-1)/2)
    NewLen = int(VERT_SIZE*6/2-(VERT_SIZE-1))
    RandEdge = random.sample(range(0,NewSize),NewLen)
    RandEdge.sort()
    for i in range(0,VERT_SIZE-2):
        for j in range(i+2,VERT_SIZE):
            RandPoint = i*VERT_SIZE - (i+1)*(i+2)/2 + 1 + j-(i+2)
            #RandPoint.sort()
            if (len(RandEdge)!=0)and(RandPoint==RandEdge[0]):
                wt = random.randint(1,100000)
                e = Edge(i,j,wt)
                E1.append(e)
                #e = Edge(j,i,wt)
                #E1.append(e)
                G1[i,j] = wt
                G1[j,i] = wt
                RandEdge.pop(0)
                
    #E1.sort(key=lambda Edge: Edge.weight, reverse = True)
              
def ConstructNewGraph2(G2,E2):
    E2.clear()
    for i in range(0,VERT_SIZE-1) :
        G2[i,i] = 0
        wt = random.randint(1,100000)
        e = Edge(i,i+1,wt)
        E2.append(e)
        #e = Edge(i+1,i,wt)
        #E2.append(e)
        G2[i,i+1] = wt
        G2[i+1,i] = wt
    
    for i in range(0,VERT_SIZE-1):
        for j in range(i+1,VERT_SIZE):
            perct = random.randint(0,4)%5
            if(perct==0):
                wt = random.randint(1,100000)
                e = Edge(i,j,wt)
                E2.append(e)
                #e = Edge(j,i,wt)
                #E2.append(e)
                G2[i,j] = wt
                G2[j,i] = wt        
                
    #E2.sort(key=lambda Edge: Edge.weight, reverse = True)



def DijkstraNoHeap(G,s,t):#0 unseen;1 fringe;2 intree
    if s==t:
        return -1
    
    status = np.random.randint(0,1,VERT_SIZE)
    dad = np.random.randint(-1,0,VERT_SIZE)
    capacity = np.zeros(VERT_SIZE)
    status[s] = 2
    for j in range(0,VERT_SIZE):
        if G[s,j]>0:
            status[j] = 1
            capacity[j] = G[s,j]
            dad[j] = s
            
    while status[t]!=2 :
        maxCap = 0
        maxCapIdx = -1
        for i in range(0,VERT_SIZE):
            if (status[i] ==1) and (capacity[i]>maxCap):
                maxCap = capacity[i]
                maxCapIdx = i
                
        v = maxCapIdx
        status[v] = 2
        for w in range(0,VERT_SIZE):
            if G[v,w]>0 :
                if status[w] ==0 :
                    status[w] =1
                    dad[w] = v
                    capacity[w] = min(capacity[v],G[v,w])
                elif (status[w] ==1) and (capacity[w]< min(capacity[v],G[v,w])):
                    dad[w] = v
                    capacity[w] = min(capacity[v],G[v,w])
                    
    return capacity[t]                



#def heapMax(heap):
#    return heap[0]

def heapifyUp(heap, value, fringe2heap, idx):
    if idx >= 1 :
        idxHalf = math.floor(idx/2)
        if value[heap[idx]] > value[heap[idxHalf]]:
            temp = heap[idxHalf]
            heap[idxHalf] = heap[idx]
            heap[idx] = temp
            fringe2heap[heap[idxHalf]] = idxHalf
            fringe2heap[heap[idx]] = idx
            heapifyUp(heap, value, fringe2heap, idxHalf)            
            
def heapifyDown(heap, value, fringe2heap, idx):
    if idx * 2 >= len(heap):
        return
    
    childIdx = 2 * idx
    if idx * 2 + 1 < len(heap):
        if value[heap[2*idx+1]]>value[heap[2*idx]]:
            childIdx = 2 * idx + 1
            
    if value[heap[childIdx]] > value[heap[idx]]:
        temp = heap[childIdx]
        heap[childIdx] = heap[idx]
        heap[idx] = temp
        fringe2heap[heap[childIdx]] = childIdx
        fringe2heap[heap[idx]] = idx
        heapifyDown(heap, value, fringe2heap, childIdx)
        
def heapInsert(heap, value, fringe2heap, vert):
    fringe2heap[vert] = len(heap)
    heap.append(vert)
    heapifyUp(heap, value, fringe2heap, len(heap) - 1)

def heapDelete(heap, value, fringe2heap, idx):
    if len(heap) > 1:
        heap[idx] = heap.pop()
        fringe2heap[heap[idx]] = idx
        if idx < len(heap)-1:
            heapifyDown(heap, value, fringe2heap, idx)
    else:
        heap = []

def DijkstraHeap( G, s, t):
    if s == t:
        return -1
    status = np.random.randint(0,1,VERT_SIZE)
    dad = np.random.randint(-1,0,VERT_SIZE)
    capacity = np.zeros(VERT_SIZE) 
    fringe2heap = np.random.randint(-1,0,VERT_SIZE)
    heap = []
    status[s] = 2
    for j in range(0,VERT_SIZE):
        if G[s,j] > 0:
            status[j] = 1
            capacity[j] = G[s,j]
            dad[j] = s
            heapInsert(heap, capacity, fringe2heap, j)
            
    while status[t] != 2:
        v = heap[0]#heapMax(heap)
        status[v] = 2
        heapDelete(heap, capacity, fringe2heap, 0)
        for w in range(0,VERT_SIZE):
            if G[v,w] > 0:
                if status[w] == 0:
                    status[w] = 1
                    dad[w] = v
                    capacity[w] = min(capacity[v], G[v,w])
                    heapInsert(heap, capacity, fringe2heap, w)
                elif (status[w] == 1) and (capacity[w] < min(capacity[v], G[v,w])):
                    dad[w] = v
                    capacity[w] = min(capacity[v], G[v,w])
                    heapifyUp(heap, capacity, fringe2heap, fringe2heap[w])
                    
    return capacity[t]



def kruskalsearch(T,capacity, s, t):
    for i in range(0,len(T)):
        e = T[i]
        if (e.vertex1 == s) and (capacity[e.vertex2] == -1):
            capacity[e.vertex2] = min(capacity[s], e.weight)
            if(e.vertex2 == t):
                return 1           
            elif (kruskalsearch(T, capacity, e.vertex2, t)!=0):
                return 1
            
        if (e.vertex2 == s) and(capacity[e.vertex1] == -1):
            capacity[e.vertex1] = min(capacity[s], e.weight)
            if(e.vertex1 == t):
                return 1            
            elif(kruskalsearch(T, capacity, e.vertex1, t)!=0):
                return 1
            
    return 0

#def MakeSet(vertex,dad,rank):
#    dad[vertex] = -1
#    rank[vertex] = 0

def Union(r1,r2,dad,rank):
    if r1 ==r2:
        return -1
    
    if rank[r1] > rank[r2] :
        dad[r2] = r1
    elif rank[r1] < rank[r2] :
        dad[r1] = r2
    else : #rank[r1] = rank[r2]
        dad[r2] = r1
        rank[r1] += 1

def Find(vertex,dad):
    w = vertex
    while dad[w] != -1:
        w = dad[w]

    return w

def Kruskal(E, s, t):
    heap = []
    cap = np.zeros(len(E)) 
    fringe2heap = np.random.randint(-1,0,len(E))
    for i in range(0,len(E)):
        cap[i] = E[i].weight
        heapInsert(heap, cap, fringe2heap, i)            
            
    #E.sort(key=lambda Edge: Edge.weight, reverse = True)
    T = []
    #############makeset###########
    rank = np.random.randint(0,1,VERT_SIZE)
    dad = np.random.randint(-1,0,VERT_SIZE)   
    
    for i in range(0,len(E)):
        if len(heap)!=0:
            v = heap[0]#heapMax(heap)
            heapDelete(heap, cap, fringe2heap, 0)
            vert1 = E[v].vertex1#E[i].vertex1
            vert2 = E[v].vertex2#E[i].vertex2
            r1 = Find(vert1,dad)
            r2 = Find(vert2,dad)
            if r1 != r2:
                T.append(E[v])#(E[i])
                Union(r1,r2,dad,rank)
            
    capacity = -1*np.ones(VERT_SIZE)
    capacity[s] = 100000
    kruskalsearch(T, capacity, s, t)
    return capacity[t]





time11 = []
time12 = []
time13 = []
time21 = []
time22 = []
time23 = []
s1 = np.random.randint(-1,0,25)
t1 = np.random.randint(-1,0,25)
s2 = np.random.randint(-1,0,25)
t2 = np.random.randint(-1,0,25)
for i in range(0,5):
    print('Data ',i+1,':\n')
################################################
####################Graph1####################
################################################    
    Graph1 = -1*np.ones((VERT_SIZE,VERT_SIZE))
    Edge1 = []
    startime10 = time.process_time()
    ConstructNewGraph1(Graph1,Edge1)
    print('Time for constructing Graph 1:',time.process_time()-startime10)
    for j in range(0,5):
        st1 = random.sample(range(0,VERT_SIZE),2)
        indx = 5*i+j
        s1[indx] = st1[0]# random.randint(0,VERT_SIZE)
        t1[indx] = st1[1]#random.randint(0,VERT_SIZE)
####################NoHeap##########################			
        startime11 = time.process_time()
        print('\n s=',s1[indx],',t=',t1[indx],',bandwidth=',DijkstraNoHeap(Graph1, s1[indx], t1[indx]))
        endtime11 = time.process_time()
        time11.append(endtime11-startime11)
####################Heap##########################			
        startime12 = time.process_time()
        print('\n s=',s1[indx],',t=',t1[indx],',bandwidth=',DijkstraHeap(Graph1, s1[indx], t1[indx]))
        endtime12 = time.process_time()
        time12.append(endtime12-startime12)
####################Kruskal##########################
        startime13 = time.process_time()
        print('\n s=',s1[indx],',t=',t1[indx],',bandwidth=',Kruskal(Edge1, s1[indx], t1[indx]))
        endtime13 = time.process_time()
        time13.append(endtime13-startime13)
    #print ('\n Dijstra no heap time:',endtime-startime)
    #print ('\n Dijstra with heap time:',endtime-startime)
    #print ('\n kruskal time:',endtime-startime)
        print('\n Dijkstra no heap time, Dijkstra with heap time, kruskal time:')
        print(endtime11-startime11,endtime12-startime12,endtime13-startime13)
        print('\n')
    print('\n\n')
################################################
####################Graph2####################
################################################    
    Graph2 = -1*np.ones((VERT_SIZE,VERT_SIZE))
    Edge2 = []
    startime20 = time.process_time()
    #Graph1, Edge1 = 
    ConstructNewGraph2(Graph2,Edge2)
    print('Time for constructing Graph 2:',time.process_time()-startime20)
    for j in range(0,5):
        indx = 5*i+j
        st2 = random.sample(range(0,VERT_SIZE),2)
        s2[indx] = st2[0]# random.randint(0,VERT_SIZE)
        t2[indx] = st2[1]#random.randint(0,VERT_SIZE)
####################NoHeap##########################
        startime21 = time.process_time()
        print('\n s=',s2[indx],',t=',t2[indx],',bandwidth=',DijkstraNoHeap(Graph2, s2[indx], t2[indx]))
        endtime21 = time.process_time()
        time21.append(endtime21-startime21)
####################Heap##########################			
        startime22 = time.process_time()
        print('\n s=',s2[indx],',t=',t2[indx],',bandwidth=',DijkstraHeap(Graph2, s2[indx], t2[indx]))
        endtime22 = time.process_time()
        time22.append(endtime22-startime22)
####################Kruskal##########################
        startime23 = time.process_time()
        print('\n s=',s2[indx],',t=',t2[indx],',bandwidth=',Kruskal(Edge2, s2[indx], t2[indx]))
        endtime23 = time.process_time()
        time23.append(endtime23-startime23)
    #print ('\n Dijstra no heap time:',endtime-startime)
    #print ('\n Dijstra with heap time:',endtime-startime)
    #print ('\n kruskal time:',endtime-startime)
        print('\n Dijkstra no heap time, Dijkstra with heap time, kruskal time:')
        print(endtime21-startime21,endtime22-startime22,endtime23-startime23)
        print('\n')
    print('\n\n\n')
    
print('Graph1: t, s, Dijkstra no heap time, Dijkstra with heap time, kruskal time\n')
for i in range(0,25):
    print('& ',t1[i],'&',s1[i],'&',round(time11[i],4),'&',round(time12[i],4),'&',round(time13[i],4),'\\\\','\n')
    
print('Graph2: t, s, Dijkstra no heap time, Dijkstra with heap time, kruskal time\n')
for i in range(0,25):
    print(' &',t2[i],'&',s2[i],'&',round(time21[i],4),'&',round(time22[i],4),'&',round(time23[i],4),'\\\\','\n')