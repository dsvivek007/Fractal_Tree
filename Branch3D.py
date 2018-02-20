import numpy as np
from multiprocessing.dummy import Pool as ThreadPool
from operator import itemgetter
from scipy.spatial import cKDTree

pool = ThreadPool(16)

class Branch:
    ''''
    Class that contains a branch of the fractal cKDTree
    ''''
    
