import numpy as np
from multiprocessing.dummy import Pool as ThreadPool
from operator import itemgetter
from scipy.spatial import cKDTree

pool = ThreadPool(16)

class Branch:
    ''''
    Class that contains a branch of the fractal cKDTree
    ''''
    def __init__(self,Init_Node,Init_Dir,l,w,Nodes,Brother_Nodes,Nsegments):

    self.Child = [0,0]
    self.Dir = np.array([0.0,0.0,0.0])
    self.Nodes =[]
    self.Queue = []
    self.Growing = True
    Shared_Node = -1
    nodes.update_collision_tree(Brother_Nodes)

    self.Nodes.append(Init_Node)
    self.Queue.append(Nodes.Nodes[Init_Node])
    grad = Nodes.Gradient(self.Queue[0])

    #### Need to decide how to introduce randomness in the initial bifurcation phase (probably theta(rand),phi(rand))
    dir = (Init_Dir+w*grad)/np.linalg.norm(Init_Dir+w*grad)

    for i in range(1,Nsegments):
            collision=nodes.collision(self.queue[i])
            if collision[1]<l/5.:
                print "Collision",i, collision
                self.growing=False
                self.queue.pop()
                self.triangles.pop()
                shared_node=collision[0]
                break
            grad=nodes.gradient(self.queue[i])
            #Project the gradient to the surface
            grad=grad-(np.dot(grad,normal))*normal
            dir=(dir+w*grad)/np.linalg.norm(dir+w*grad)
        nodes_id=nodes.add_nodes(self.queue[1:])
        [self.nodes.append(x) for x in nodes_id]
        if not self.growing:
            nodes.end_nodes.append(self.nodes[-1])
        self.dir=dir
       # #print self.triangles
    #Uncomment the following lines for a closed network
     #   if shared_node is not -1:
      #      self.nodes.append(shared_node)
