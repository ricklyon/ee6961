import numpy as np


def dftmtx(N):
    dftmtx = np.fft.fft(np.eye(N))
    return dftmtx

def non_square_dftmtx(data, L, K):
    V = np.empty([data.shape[0], L], dtype=np.complex64)
    for i in range(0, data.shape[0]):
        for j in range(0, L):
            test = data[i]
            V[i][j] = np.exp(-1j*2*np.pi*(data[i]*(j/K)))
    return V
        
