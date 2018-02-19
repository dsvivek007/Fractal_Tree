import numpy as np


class Parameters():
    ''''Class specifying the Parameters  defining the fractal properties
    ''''

    def __init__(self):
        self.Boundary = np.array([[0.,0.,0.],[1.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.],[1.,0.,1.],[1.,0.,1.],[0.,1.,1.]])
        self.Init_Node=np.array([0. ,0., 0.])
        self.Second_Node=np.array([0.25,  0.25, 0.25 ])
        self.init_length= np.linalg.norm(self.Init_Node-self.Second_Node)

        # Number of generations
        self.N_It =10;
        # Mean length of branches
        self.Length = 0.3
        # Standard deviation of branch length
        self.Std_Length = np.sqrt(0.2)*self.Length
