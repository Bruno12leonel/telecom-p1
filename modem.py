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
        self.buffer = self.bufsz*[0]
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
        self.data = data
        

    def get_bits(self):
        s = np.concatenate((self.buffer,self.data))
        v0i = np.zeros(len(s))
        v0r = np.zeros(len(s))
        v1i = np.zeros(len(s))
        v1r = np.zeros(len(s))
        count = 46*self.fs//48000
        L = self.bufsz
        T = 1/self.fs
        r = 0.999

        for n in range(L,len(s)):
            v0r[n] = s[n] - r**L * np.cos(self.rx_omega0*L*T)*s[n-L] + r*np.cos(self.rx_omega0*T)*v0r[n-1] - r*np.sin(self.rx_omega0*T)*v0i[n-1]
            v0i[n] = -r**L*np.sin(self.rx_omega0*L*T)*s[n-L] + r*np.cos(self.rx_omega0*T)*v0i[n-1] + r*np.sin(self.rx_omega0*T)*v0r[n-1]
            v1r[n] = s[n] - r**L * np.cos(self.rx_omega1*L*T)*s[n-L] + r*np.cos(self.rx_omega1*T)*v1r[n-1] - r*np.sin(self.rx_omega1*T)*v1i[n-1]
            v1i[n] = -r**L*np.sin(self.rx_omega1*L*T)*s[n-L] + r*np.cos(self.rx_omega1*T)*v1i[n-1] + r*np.sin(self.rx_omega1*T)*v1r[n-1]            

        c = abs(v1r**2+v1i**2-v0r**2-v0i**2)                   
        v = np.zeros(len(c))
        y = np.zeros(len(c))
        d = v1r**2+v1i**2-v0r**2-v0i**2    
        for n in range(L,len(s)):
            v[n] = (1-r)*c[n] + 2*r*np.cos(2*np.pi*300/self.fs)*v[n-1] - r**2*v[n-2]
            y[n] = v[n] - v[n-2]
            

        filt = signal.firwin(40, 300, pass_zero='lowpass', fs= self.fs)
        amostra = 1*((y[1:]<0)&(y[:-1]>=0))
        for i in range(len(amostra)):
            if amostra[i] != 0:
                self.bits.append(1 if d[i] > 0 else 0)
        

        '''
        

        convolvedArr = np.convolve(amostra, filt)
        for i in range(1,len(y),1):
            if y[i-1]<0 and y[i]>=0:
                if len(convolvedArr) > i + count:                                        
                    self.bits.append(1 if convolvedArr[i + count] > 0 else 0)
                else:    
                    self.bits.append(1 if convolvedArr[-1] > 0 else 0)
        self.data = self.data[L:]            
        '''
        self.buffer = s[:L]
        return self.bits