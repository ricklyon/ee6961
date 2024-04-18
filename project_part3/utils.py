import numpy as np


def dftmtx(N):
    dftmtx = np.fft.fft(np.eye(N))
    return dftmtx