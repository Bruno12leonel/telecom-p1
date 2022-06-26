import numpy as np

class Modem:
    def __init__(self, fs, bufsz, ans=False):
        self.fs = fs  # taxa de amostragem
        self.amostras = []
        self.bufsz = bufsz  # quantidade de amostas que devem ser moduladas por vez
        # frequências de modulação (upload)
        self.tx_omega0 = 2*np.pi*(1080 + 100)
        self.tx_omega1 = 2*np.pi*(1080 - 100)
        # frequências de demodulação (download)
        self.rx_omega0 = 2*np.pi*(1750 + 100)
        self.rx_omega1 = 2*np.pi*(1750 - 100)
        self.y = []
        # se o modem estiver atendendo uma ligação
        if ans:
            # inverte as frequências
            self.tx_omega0, self.rx_omega0 = self.rx_omega0, self.tx_omega0
            self.tx_omega1, self.rx_omega1 = self.rx_omega1, self.tx_omega1

    # Modulação

    def put_bits(self, bits):
        for bit in bits:
            omega = 2*np.pi*(1180 if bit==0 else 980)
            for i in range(self.fs//300):
                phi += omega/self.fs
                self.y.append(np.sin(phi))

    def get_samples(self):
        if not self.y: self.y = [1] * self.bufsz
        res = self.y[0:self.bufsz]
        self.y = self.y[self.bufsz:]
        
        return res
        #return np.zeros(self.bufsz)

    # Demodulação

    def put_samples(self, data):
        pass

    def get_bits(self):
        return []
