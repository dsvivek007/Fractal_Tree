import numpy as np
from multiprocessing.dummy import Pool as ThreadPool
from operator import itemgetter
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#######################################################################################
class Parameters():

    def __init__(self):
        self.Boundary = np.array([[0.,0.,0.],[1.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.],[1.,0.,1.],[1.,0.,1.],[0.,1.,1.]])
        self.init_node=np.array([0. ,0.,0.])
        self.second_node=np.array([0.25,  0.25, 0.25 ])
        self.init_length= np.linalg.norm(self.init_node-self.second_node)
        self.filename='test'
        # Number of generations
        self.N_it =1;
        # Mean length of branches
        self.length = 0.3
        # Standard deviation of branch length
        self.std_length = np.sqrt(0.2)*self.length
        # Minmum length: prevents negative lengths
        self.min_length = self.length/10
        self.branch_angle = 45.0
        self.w = 0.1  # Weight for gradient of distance
        self.c = 0.1  # Randomness parameter for directipo to grow
        # Approximate length of each segment
        self.l_segment = 0.01


        self.Fascicles=True
        ###########################################
        # Fascicles data
        ###########################################
        self.fascicles_angles=[-1.5,.2] #rad
        self.fascicles_length=[.5,.5]
        # Save data?
        self.save=True
        self.save_paraview=True

###################################################################################

class Branch:

    def __init__(self,init_node,init_dir,l,w,c,n_segments):
        self.dir = np.array([0.0,0.0,0.0])
        self.nodes = []
        self.queue = []
        self.growing = True
        shared_node  = -1

        self.nodes.append(init_node)
        self.queue.append(nodes.nodes[init_node])

        #Generate a growth direction from init_dir by perturbing it with a random vector
        dir= c*np.random.randn(3)+init_dir
        dir=dir/np.linalg.norm(dir)


        grad=nodes.gradient([self.queue[0]])

        dir=dir+w*np.transpose(grad)

        dir=dir/np.linalg.norm(dir)
        for i in range(1,n_segments+1):
            self.add_node_to_queue(self.queue[i-1],dir*l/n_segments)
            grad=nodes.gradient(self.queue[i])
            dir=dir+w*np.transpose(grad)
            dir=dir/np.linalg.norm(dir)


        nodes_id=nodes.add_nodes(self.queue[1:])
        [self.nodes.append(x) for x in nodes_id]
        if not self.growing:
            nodes.end_nodes.append(self.nodes[-1])
        self.dir=dir
    def add_node_to_queue(self,init_node,dir):
        point = init_node+dir

        self.queue.append(point)
        return 0






#######################################################################################
class Nodes:

    def __init__(self,init_node):
        self.nodes=[]
        self.nodes.append(init_node)
        self.last_node=0
        self.end_nodes=[]
        self.tree=cKDTree(self.nodes)

    def add_nodes(self,queue):
        """This function stores a list of nodes of a branch and returns the node indices. It also updates the tree to compute distances.

        Args:
            queue (list): a list of arrays containing the coordinates of the nodes of one branch.

        Returns:
            nodes_id (list): the indices of the added nodes.
        """
        nodes_id=[]
        #print(queue[0])
        #for point in queue:
        for i in range(0,len(queue)):
            point = np.array(queue[i])
            #print(point)
            self.nodes.append(point)
            self.last_node+=1
            nodes_id.append(self.last_node)
        for i in range(0,len(self.nodes)):
            self.tree=cKDTree(self.nodes[1])
        return nodes_id
    def distance_from_point(self,point):
        d,node=self.tree.query(point)
        return d

    def gradient(self,point):
        """This function returns the gradient of the distance from the existing points of the tree from any point. It uses a central finite difference approximation.
            Args:
            point (array): the coordinates of the point to calculate the gradient of the distance from.
            Returns:
            grad (array): (x,y,z) components of gradient of the distance.
            """
        delta=0.01
        dx=np.array([delta,0,0])
        dy=np.array([0.0,delta,0.0])
        dz=np.array([0.0,0.0,delta])
        distx_m=self.distance_from_point(point-dx)
        distx_p=self.distance_from_point(point+dx)
        disty_m=self.distance_from_point(point-dy)
        disty_p=self.distance_from_point(point+dy)
        distz_m=self.distance_from_point(point-dz)
        distz_p=self.distance_from_point(point+dz)
        grad= np.array([(distx_p-distx_m)/(2*delta),(disty_p-disty_m)/(2*delta),(distz_p-distz_m)/(2*delta)])
        return grad
#######################################################################################################
param = Parameters()
nodes=Nodes(param.init_node)
#test = Branch(0,np.array([1.,1.,1.]),1.0,0.1,0.1,100)
#print(nodes.nodes[test.nodes[100]])
#####################################################################################################
def Fractal_Tree_3D(param):
    #Define the initial direction
    init_dir=(param.second_node-param.init_node)/np.linalg.norm(param.second_node-param.init_node)

    #Initialize the nodes object, contains the nodes and all the distance functions
    #nodes=Nodes(param.init_node)

    #Initialize the dictionary that stores the branches objects
    branches={}
    last_branch=0
    branches[last_branch]=Branch(0,init_dir,param.init_length,0.0,0.0,int(param.init_length/param.l_segment))
    branches_to_grow=[]
    branches_to_grow.append(last_branch)
    ien=[]
    for i_n in range(len(branches[last_branch].nodes)-1):
        ien.append([branches[last_branch].nodes[i_n],branches[last_branch].nodes[i_n+1]])






    if param.save:
        if param.save_paraview:
            from ParaviewWriter import write_line_VTU
            print ('Finished growing, writing paraview file')
            xyz=np.zeros((len(nodes.nodes),3))
            for i in range(len(nodes.nodes)):
                print(nodes.nodes[i])
                xyz[i,:]=nodes.nodes[i]

            write_line_VTU(xyz, ien, param.filename + '.vtu')

        np.savetxt(param.filename+'_ien.txt',ien,fmt='%d')
        np.savetxt(param.filename+'_xyz.txt',xyz)
        np.savetxt(param.filename+'_endnodes.txt',nodes.end_nodes,fmt='%d')

###################################################################################################
branches = Fractal_Tree_3D(param)
#print(nodes.nodes)
