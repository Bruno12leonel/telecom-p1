from random import sample
from this import d
import numpy as np
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

        self.filt = signal.firwin(40, 300, pass_zero='lowpass', fs= self.fs)
        
        self.amostras = []
        self.phi = 0
        self.data = []
        self.buffer = self.bufsz*[0]
        
        self.vorAnterior = 0
        self.voiAnterior = 0
        self.v1rAnterior = 0
        self.v1iAnterior = 0
        self.vanterior = 0
        self.vantanterior = 0
        self.yant = 0
        self.listD = []
        self.counter_state = 'idle'
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
        self.data.extend(data)
        

    def get_bits(self):
        s = np.concatenate((self.buffer,self.data))        
        L = self.bufsz
        T = 1/self.fs
        r = 0.99999
        bits = []

        #depois do quarto teste ele passa a ser de 320 o len(s)
        for n in range(L,len(s)):
            v0r = s[n] - r**L * np.cos(self.rx_omega0*L*T)*s[n-L] + r*np.cos(self.rx_omega0*T)*self.vorAnterior - r*np.sin(self.rx_omega0*T)*self.voiAnterior
            v0i = -r**L*np.sin(self.rx_omega0*L*T)*s[n-L] + r*np.cos(self.rx_omega0*T)*self.voiAnterior + r*np.sin(self.rx_omega0*T)*self.vorAnterior
            v1r = s[n] - r**L * np.cos(self.rx_omega1*L*T)*s[n-L] + r*np.cos(self.rx_omega1*T)*self.v1rAnterior - r*np.sin(self.rx_omega1*T)*self.v1iAnterior
            v1i = -r**L*np.sin(self.rx_omega1*L*T)*s[n-L] + r*np.cos(self.rx_omega1*T)*self.v1iAnterior + r*np.sin(self.rx_omega1*T)*self.v1rAnterior            

            c = abs(v1r**2+v1i**2-v0r**2-v0i**2)                   
            

            v = (1-r)*c + 2*r*np.cos(2*np.pi*300/self.fs)*self.vanterior - r**2*self.vantanterior

            y = v - self.vantanterior

            if self.counter_state == 'running':
                self.listD.append(v1r**2+v1i**2-v0r**2-v0i**2)
                self.counter -= 1
            if self.counter_state == 'idle' and (y >= 0) and (self.yant < 0):
                self.counter = 46*48000//self.fs
                self.listD = []
                self.counter_state = 'running'               
            elif self.counter_state == 'running' and self.counter == 0:
                bits.append(1 if np.dot(self.filt,self.listD[-40:]) > 0 else 0)                
                self.counter_state = 'idle'              

            self.vorAnterior = v0r
            self.voiAnterior = v0i
            self.v1rAnterior = v1r
            self.v1iAnterior = v1i            
            self.vantanterior = self.vanterior
            self.vanterior = v
            self.yant = y       
        
        
        self.buffer = s[-L:]
        self.data = []
        return bits
               
        