import numpy as np

class Target():
    def __init__(self, N_col):
        self.z = (np.zeros(N_col, dtype=np.float64) 
                   + np.zeros(N_col, dtype=np.float64)*1j)
                   
class Sources():
    def __init__(self, Nsour):
        self.z = (np.zeros(Nsour, dtype=np.float64) 
                   + np.zeros(Nsour, dtype=np.float64)*1j)
        self.w = np.zeros(Nsour, dtype=np.float64)
