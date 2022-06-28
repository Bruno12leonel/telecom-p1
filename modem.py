from random import sample
import numpy as np
import matplotlib.pyplot as plt

class Modem:
    def __init__(self, fs, bufsz, ans=False):        
        self.fs = fs  # taxa de amostragem
        self.bufsz = bufsz  # quantidade de amostas que devem ser moduladas por vez
        # frequências de modulação (upload)
        self.tx_omega0 = 2*np.pi*(1080 + 100)
        self.tx_omega1 = 2*np.pi*(1080 - 100)
        # frequências de demodulação (download)
        self.rx_omega0 = 2*np.pi*(1750 + 100)
        self.rx_omega1 = 2*np.pi*(1750 - 100)
        # se o modem estiver atendendo uma ligação
        
        self.amostras = []
        self.phi = 0
        if ans:
            # inverte as frequências
            self.tx_omega0, self.rx_omega0 = self.rx_omega0, self.tx_omega0
            self.tx_omega1, self.rx_omega1 = self.rx_omega1, self.tx_omega1

    # Modulação
    def put_bits(self, bits):       
        for bit in bits:            
            y = np.zeros(self.bufsz)            
            t = 0.0
            omega = (self.tx_omega0 if bit==0 else self.tx_omega1)
            self.phi = self.phi - omega*t            
            for i in range(self.bufsz):
                y[i] = np.sin(omega*t + self.phi)
                t = t + 1/self.fs
            self.phi = omega*t + self.phi
            self.amostras.append(y)
    
    def get_samples(self):
        if self.amostras == []:
            self.put_bits([1])
        res = self.amostras[0]
        self.amostras = self.amostras[1:]
        return res

    # Demodulação
    def put_samples(self, data):
        return 0

    def get_bits(self):
        return []
