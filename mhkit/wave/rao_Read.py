import numpy as np
import mhkit.wave.performance as p
import scipy.io as sio
import matplotlib.pyplot as plt

# WAMIT DATA
filename = 'C:/Users/akeeste/Desktop/temp/MHKiT rao test/cyl.4'
data = np.loadtxt(filename, skiprows=1)
w = data[:,0]
nw = len(np.unique(w))

w            = np.reshape(data[:,0],[nw,6])
w = w[:,0]
direction    = np.reshape(data[:,1],[nw,6])
dof          = np.reshape(data[:,2],[nw,6])
raoExp_mag   = np.reshape(data[:,3],[nw,6])
raoExp_phase = np.reshape(data[:,4],[nw,6])
raoExp_real  = np.reshape(data[:,5],[nw,6])
raoExp_imag  = np.reshape(data[:,6],[nw,6])


# BEM DATA
filename2 = 'C:/Users/akeeste/Desktop/temp/MHKiT rao test/bem.mat'
bem_data = sio.loadmat(filename2)
bem_data = bem_data ['data']
vol = bem_data['vol'][0][0][0][0]
rho = bem_data['rho'][0][0][0][0]

mass = vol*rho
M = np.zeros([6,6])
M[0][0] = mass
M[1][1] = mass
M[2][2] = mass

A = bem_data['A'][0][0]
B = bem_data['B'][0][0]
C = bem_data['Khs'][0][0]

Fre = bem_data['ex_re'][0][0]
Fim = bem_data['ex_im'][0][0]
# Fre = bem_data['sc_re'][0][0]
# Fim = bem_data['sc_im'][0][0]
# Fre = bem_data['fk_re'][0][0]
# Fim = bem_data['fk_im'][0][0]
F = Fre + Fim*1j

raoCalc = np.zeros([len(w),6])
for j in range(1,len(w)):
    denom = C - (M + A[:,:,j])*w[j]**2 + 1j*B[:,:,j]*w[j]
    raoCalc[j,:] = np.matmul(np.transpose(F[:,:,j]), np.linalg.inv(denom))


i = 2
plt.plot(w, raoExp_mag[:,i], w, abs(raoCalc[:,i]), '--')



# D = C - (M+A)*w**2 + 1j*B*w

# i = 2
# raoCalc = p.response_amplitude_operator(F[i,0,:], M[i,i], A[i,i,:], B[i,i,:], C[i,i], w[:,i])
# raoCalc = p.response_amplitude_operator(F, M, A, B, C, w)

# overhaul RAO calc to use DOFs correctly. Reference summation in wanan sheng book
# plt.plot(w, raoExp_mag[:,i], w, raoCalc['magnitude'], '--')