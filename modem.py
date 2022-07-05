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
        self.data = []
        ##
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
        self.data.extend(data)
        

    def get_bits(self):
        s = self.data
        v0i=np.zeros(len(s))
        v0r=np.zeros(len(s))
        v1i=np.zeros(len(s))
        v1r=np.zeros(len(s))
        L = self.bufsz
        T = 1/self.fs

        r=0.99

        for n in range(1,len(s)):
            v0r[n] = s[n] - r**L * np.cos(self.rx_omega0*L*T)*s[n-L] + r*np.cos(self.rx_omega0*T)*v0r[n-1] - r*np.sin(self.rx_omega0*T)*v0i[n-1]
            v0i[n] = -r**L*np.sin(self.rx_omega0*L*T)*s[n-L] + r*np.cos(self.rx_omega0*T)*v0i[n-1] + r*np.sin(self.rx_omega0*T)*v0r[n-1]
            v1r[n] = s[n] - r**L * np.cos(self.rx_omega1*L*T)*s[n-L] + r*np.cos(self.rx_omega1*T)*v1r[n-1] - r*np.sin(self.rx_omega1*T)*v1i[n-1]
            v1i[n] = -r**L*np.sin(self.rx_omega1*L*T)*s[n-L] + r*np.cos(self.rx_omega1*T)*v1i[n-1] + r*np.sin(self.rx_omega1*T)*v1r[n-1]
        
            rho = v1r**2+v1i**2+v0r**2+v0i**2   # carrier detection

            c = abs(v1r**2+v1i**2-v0r**2-v0i**2)
            v = np.zeros(len(c))
            y = np.zeros(len(c))
            r=0.9999
            for n in range(1,len(s)):
                v[n] = (1-r)*c[n] + 2*r*np.cos(2*np.pi*300/self.fs)*v[n-1] - r**2*v[n-2]
                y[n] = v[n] - v[n-2]

        return y
