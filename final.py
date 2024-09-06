import numpy as np
import matplotlib.pyplot as plt
from scipy import signal


h = 0.001 #krok calkowania
czas_symulacji = 10.0
pi = np.pi


a3, a2, a1, a0 = 0, 0, 0, 0  
b3, b2, b1, b0 = 0, 0, 0, 0  

print("Transmitancja naszego układu= (a1*s + a0)/(b2*s^2+b1*s+b0)")
print("Transmitancja regulatora PI: kp + ki/s")
a11 = float(input("Wartosc parametru a1: "))
a01 = float(input("Wartosc parametru a0: "))
b21 = float(input("Wartosc parametru b2: "))
b11 = float(input("Wartosc parametru b1: "))
b01 = float(input("Wartosc parametru b0: "))
kp = float(input("Wartosc parametru kp: "))
ki = float(input("Wartosc parametru ki: "))
ampl = float(input("Wartosc parametru amplitudy:"))
f = float(input("Liczba okresow sygnalu w czasie symulacji:"))

a3 = 1
a2 = (b11 + (a11 * kp))/b21
a1 = (b01 + (a11 * ki) + (a01 * kp))/b21
a0 = (a01 * ki)/b21
b3 = 0
b2 = (a11 * kp)/b21
b1 = ((a11 * ki) + (a01 * kp))/b21
b0 = (a01 * ki)/b21


A = np.array([[0, 1, 0],
              [0, 0, 1],
              [-a0, -a1, -a2]])

B = np.array([0, 0, 1])

C = np.array([b0, b1, b2])

D = 0



total = int(czas_symulacji / h)

w = 2 * pi * f / czas_symulacji

t = np.arange(0, czas_symulacji, h)


# zerowe warunki poczatkowe
xi1_sin = np.array([0, 0, 0])
xi1_sq = np.array([0, 0, 0])
xi1_tri = np.array([0, 0, 0])
xi2_sin = np.array([0, 0, 0])
xi2_sq = np.array([0, 0, 0])
xi2_tri = np.array([0, 0, 0])

y_sin = np.zeros(total)
y_sq = np.zeros(total)
y_tri = np.zeros(total)

#sygnaly wejsciowe
us = []
uf = []
ut = []

for i in range(total):
    us.append(ampl * np.sin(w * h * i))
    uf.append(np.where(us[i] > 0, ampl, -ampl))
    ut.append(ampl * signal.sawtooth( (w*h*i)+pi/2, 0.5))



for i in range(total):
    Ax_sin = np.dot(A, xi1_sin)
    Bu_sin = B * us[i]
    Cx_sin = np.dot(C, xi1_sin)
    Du_sin = D * us[i]
    xi_sin = (h/2) * (xi2_sin + (Ax_sin + Bu_sin))
    xi_sin = xi_sin + xi1_sin
    xi1_sin = xi_sin
    xi2_sin = Ax_sin + Bu_sin
    y_sin[i] = Cx_sin + Du_sin

    Ax_sq = np.dot(A, xi1_sq)
    Bu_sq = B * uf[i]
    Cx_sq = np.dot(C, xi1_sq)
    Du_sq = D * uf[i]
    xi_sq = (h/2) * (xi2_sq + (Ax_sq + Bu_sq))
    xi_sq = xi_sq + xi1_sq
    xi1_sq = xi_sq
    xi2_sq = Ax_sq + Bu_sq
    y_sq[i] = Cx_sq + Du_sq

    Ax_tri = np.dot(A, xi1_tri)
    Bu_tri = B * ut[i]
    Cx_tri = np.dot(C, xi1_tri)
    Du_tri = D * ut[i]
    xi_tri = (h/2) * (xi2_tri + (Ax_tri + Bu_tri))
    xi_tri = xi_tri + xi1_tri
    xi1_tri = xi_tri
    xi2_tri = Ax_tri + Bu_tri
    y_tri[i] = Cx_tri + Du_tri


plt.figure(figsize=(10, 6))

plt.subplot(3, 1, 1)
plt.plot(t, us, label='Sygnal harmoniczny')
plt.plot(t, uf, label='Sygnal prostokatny')
plt.plot(t, ut, label='Sygnal trojkatny')
plt.legend()
plt.title('Sygnaly wejsciowe')
plt.xlabel('Czas (s)')
plt.ylabel('Amplituda')
plt.grid(True)

plt.subplot(3, 1, 2)
plt.plot(t, y_sin, label='Wyjsciowy harmoniczny')
plt.plot(t, y_sq, label='Wyjsciowy prostokatny')
plt.plot(t, y_tri, label='Wyjsciowy trojkatny')
plt.legend()
plt.title('Sygnaly wyjsciowe')
plt.xlabel('Czas (s)')
plt.ylabel('Amplituda')
plt.grid(True)

plt.tight_layout()
plt.show()
