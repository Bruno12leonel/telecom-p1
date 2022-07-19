from random import sample
from this import d
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

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
        self.data = self.bufsz*[0]
        self.bits = []
        self.bits_filt = []
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
            #sin(wt + phi)
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
        count = 46*self.fs//48000
        L = self.bufsz
        T = 1/self.fs
        r = 0.9999

        v0i = np.zeros(len(self.data))
        v0r = np.zeros(len(self.data))
        v1i = np.zeros(len(self.data))
        v1r = np.zeros(len(self.data))
        
        for n in range(1,len(self.data)):
            v0r[n] = self.data[n] - r**L * np.cos(self.rx_omega0*L*T)*self.data[n-L] + r*np.cos(self.rx_omega0*T)*v0r[n-1] - r*np.sin(self.rx_omega0*T)*v0i[n-1]
            v0i[n] = -r**L*np.sin(self.rx_omega0*L*T)*self.data[n-L] + r*np.cos(self.rx_omega0*T)*v0i[n-1] + r*np.sin(self.rx_omega0*T)*v0r[n-1]
            v1r[n] = self.data[n] - r**L * np.cos(self.rx_omega1*L*T)*self.data[n-L] + r*np.cos(self.rx_omega1*T)*v1r[n-1] - r*np.sin(self.rx_omega1*T)*v1i[n-1]
            v1i[n] = -r**L*np.sin(self.rx_omega1*L*T)*self.data[n-L] + r*np.cos(self.rx_omega1*T)*v1i[n-1] + r*np.sin(self.rx_omega1*T)*v1r[n-1]            

        rho = v1r**2+v1i**2+v0r**2+v0i**2
        c = abs(v1r**2+v1i**2+v0r**2+v0i**2)    
        v = np.zeros(len(c))
        y = np.zeros(len(c))

        for n in range(1, len(self.data)):
            v[n]=(1-r)*c[n]+2*r*np.cos(2*np.pi*300/self.fs)*v[n-1]-r**2*v[n-2]
            y[n]=v[n]-v[n-2]

        filt= signal.firwin(40, 300, pass_zero='lowpass', fs=self.fs)
        amostra=1*((y[1:]>0)&(y[:-1]<0))
        for i in range(len(amostra)):
            if amostra[i]!=0:
                self.bits_filt.append(1 if d[i]>0 else 0)
        # return self.bits_filt
        

    def get_bits(self):
        return self.bits_filt