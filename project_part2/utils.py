import numpy as np
import matplotlib.pyplot as plt

def qpsk_encoder(data):
    """Generate QPSK symbols from data sequence. Symbols are complex valued at baseband."""

    # 4 symbols per byte of data
    symbols = np.zeros(int(len(data) * 4), dtype='complex128')
    # look up table for qpsk complex sybmols. Counterclockwise from first quadrant: (00, 01, 11, 10)
    qpsk_lut = (1 / np.sqrt(2)) * np.array([1+1j, -1+1j, 1-1j, -1-1j])
    # QPSK encoder. loop over each byte in the data sequence
    for i, d in enumerate(data):
        # 4 symbols per byte, loop over each group of 2 bits
        for sym_i in range(4):
            # get the binary data for the current symbol
            symb_data = (d >> (sym_i*2)) & 0x3
            # map the symbol data into a complex value, and put into the symbol array
            symbols[i*4 + (3-sym_i)] = qpsk_lut[symb_data]

    return symbols

def qpsk_decoder(symbols):
    """ Return binary data from QPSK symbol sequence"""

    # reshape into rows of 4 symbols (corresponding to one byte each)
    symbols = np.reshape(symbols, (len(symbols) // 4, 4))
    data_bytes = np.zeros(len(symbols), dtype=np.uint8)
    
    # loop over each group of 4 symbols
    for i, sym_byte in enumerate(symbols):
        # loop over each symbol in the byte
        for j, d in enumerate(sym_byte):
            # convert qpsk symbol back to data
            d_real, d_imag = np.real(d), np.imag(d)

            if d_real > 0 and d_imag > 0:
                data_j = 0x00 
            elif d_real < 0 and d_imag > 0:
                data_j = 0x01 
            elif d_real < 0 and d_imag < 0:
                data_j = 0x03
            else:# d_real > 0 and d_imag < 0:
                data_j = 0x02
            # put the symbol in the appropriate spot in the data byte
            data_bytes[i] |= (data_j  << 2*(3-j))

    return bytes(data_bytes)

def dtft(xn: np.ndarray, f: np.ndarray, fs: float):
    """
    Compute the DTFT of the discrete time signal (x[n]) over a continuous time frequency range (f)
    [cyles per second]. This requires the sampling frequency (fs) that the x[n] is sampled with.

    Note that if the requested frequencies are greater than fs/2, the spectrum will alias.
    """
    # convert the continuous time frequency into a discrete frequency range. The discrete frequencies
    # are bounded by -0.5 to 0.5.
    # To convert the continuous time frequency (cycles / sec) into the discrete frequency (cycles / sample). 
    # divide it by the sampling rate fs (samples / sec):
    # (cycles / sec) / (samples / sec) = (cycles / sample)
    fn = f / fs
    omega_n = 2 * np.pi * fn

    # create meshgrid with n along the rows and omega in the columns
    n = np.arange(len(xn))
    omega_mesh, n_mesh = np.meshgrid(omega_n, n)

    # broadcast input sequence across all omega
    x_b = np.broadcast_to(xn[..., None], (len(n), len(omega_n)))
    # sum across all n
    Xw = np.sum(x_b * np.exp(-1j* omega_mesh * n_mesh), axis=0)

    return Xw

def srrc_pulse(n: np.ndarray, sps: int, beta: float):
    """
    Square root raised cosine pulse. 
	
    This pulse is given the time domain here: https://en.wikipedia.org/wiki/Root-raised-cosine_filter
	
	The T in this equation is the symbol period.
    To convert to discrete time, replace t with nTs (Ts is symbol period, T is sampling period, sps is samples per symbol),
    t/Ts -> nT/Ts = nT/(T*sps) = n/sps

    Parameters:
    -----------
    n: np.ndarray
        vector of discrete sample indices
    sps: int
        samples per symbol
    beta: float
        roll off factor
    """
    with np.errstate(divide='ignore', invalid='ignore'):
        srrc_num = (np.sin(np.pi*n*(1-beta)/sps)) + (4*beta*n/sps)*(np.cos(np.pi*n*(1+beta)/sps)) 
        srrc_denom = (np.pi*n/sps) * (1 - (4*beta*n/sps)**2)
        srrc = srrc_num / srrc_denom

    # the above expression evaluate to +/- inf when t = +/- Ts / (4 * beta), or
    # n = sps / (4 * beta)
    value_1 = (beta/np.sqrt(2))*( (1+(2/np.pi))*np.sin(np.pi/(4*beta)) + (1-(2/np.pi))*np.cos(np.pi/(4*beta)) )
    srrc = np.where(np.abs(n) == sps / (4 * beta), value_1, srrc)

    # set the value at n=0 (also NaN in the equation above)
    value_0 = (1 + beta * ((4 / np.pi) - 1))
    srrc = np.where(n == 0, value_0, srrc)

    return (1 / np.sqrt(sps)) * srrc

def plot_fft(x, fs, ax = None, figsize=None):
    
    if ax is None:
        fig, ax = plt.subplots(1,1, figsize=figsize)

    Xf = np.fft.fft(x)
    f = (1 / len(Xf))
    fvec = np.linspace(-0.5, 0.5 - f, len(Xf)) * fs
    
    ax.plot(fvec / 1e3, np.fft.fftshift(np.abs(Xf)))
    ax.set_xlabel('Frequency [kHz]')
    ax.set_ylabel('')
    return ax