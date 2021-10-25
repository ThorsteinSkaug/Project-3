import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Importing from txt file made in the cpp code

l = []
df = pd.read_csv('coordinates_an.txt', sep=' ', header=None)
l.append(df[0].to_numpy())
print(l)

def read(filename):
    df = pd.read_csv(filename, sep=' ', header=None)
    var = []
    for i in range(len(df.columns)):
        var.append(df[i].to_numpy())
    return var

# 1.

pos_an =  read('coordinates_an.txt')
t001 = pos_an[0]
pos_eu = read('coordinates_euler.txt')
pos_rk = read('coordinates_rk.txt')


fig,ax = plt.subplots(1,2, figsize=(20,10))

fig.suptitle('Particle motion in Penning Trap')

ax[0].plot(pos_eu[1],pos_eu[2], '--', color='r', label='Euler')
ax[0].plot(pos_rk[1],pos_rk[2], '--', color='b', label='RK4')
ax[0].plot(pos_an[1],pos_an[2], '-', color='k', label='Analytical')
ax[0].plot(pos_an[1][0], pos_an[2][0], 'o', markersize=8, label='Starting position')
ax[0].plot(pos_an[1][-1], pos_an[2][-1], 'o', markersize=8, label='Final position')
ax[0].plot(pos_rk[1][-1], pos_rk[2][-1], 'o', markersize=8, label='Final position_rk')
ax[0].plot(pos_eu[1][-1], pos_eu[2][-1], 'o', markersize=8, label='Final position_eu')
ax[0].grid()
ax[0].set_ylabel(r'y [$\mu m$]')
ax[0].set_xlabel(r'x [$\mu m$]')
ax[0].legend()
ax[0].axis('equal')

ax[1].plot(t001,pos_eu[3], '--', color='r', label='Euler')
ax[1].plot(t001,pos_rk[3], '--', color='b', label='RK4')
ax[1].plot(t001,pos_an[3], '--', color='k', label='Analytical')
ax[1].grid()
ax[1].set_ylabel(r'z [$\mu m$]')
df = pd.read_csv('resonance.txt', sep=' ', header=None)
w_V = df[0].to_numpy()
f_1 = df[1].to_numpy()
f_2 = df[2].to_numpy()
f_3 = df[3].to_numpy()
ax[1].set_xlabel(r't [$\mu s$]')
ax[1].legend()

fig.tight_layout()
plt.savefig('1particle_motion.pdf', dpi=900)
plt.show()

# 2.

rk2_no =  read('rk2particles_without.txt')
rk2 =  read('rk2particles_with.txt')


fig, ax = plt.subplots(1,2, figsize=(20,10))

fig.suptitle('Particle motions in the xy-plane')

ax[0].plot(rk2_no[1],rk2_no[2], '--', color='r', label='First particle')
ax[0].plot(rk2_no[7],rk2_no[8], '--', color='b', label='Second Particle')
ax[0].grid()
ax[0].set_ylabel(r'x [$\mu m$]')
ax[0].set_xlabel(r'y [$\mu m$]')
ax[0].legend()
ax[0].axis('equal')
ax[0].set_title('Without Particle interactions')

ax[1].plot(rk2[1],rk2[2], '--', color='r', label='First particle')
ax[1].plot(rk2[7],rk2[8], '--', color='b', label='Second Particle')
ax[1].grid()
ax[1].set_xlabel(r'x [$\mu m$]')
ax[1].set_ylabel(r'y [$\mu m$]')
ax[1].legend()
ax[1].axis('equal')
ax[1].set_title('With Particle interactions')

fig.tight_layout()
plt.savefig('2particle_motion.pdf', dpi=900)
plt.show()

# 3.

fig,ax = plt.subplots(1,3, figsize=(30,10))

fig.suptitle('Phase space in Penning Trap without particle interactions')

for i in range(3):
    ax[i].plot(rk2_no[i+1],rk2_no[i+4], '--', color='r', label='First particle')
    ax[i].plot(rk2_no[i+7],rk2_no[i+10], '--', color='b', label='Second Particle')
    ax[i].grid()
    ax[i].set_ylabel(r'x [$\mu m$]')
    ax[i].set_xlabel(r'$v_x$ [$\mu m / s$]')
    ax[i].legend()

fig.tight_layout()
plt.savefig('phase_spaces_without.pdf', dpi=900)
plt.show()


fig,ax = plt.subplots(1,3, figsize=(30,10))

fig.suptitle('Phase space in Penning Trap with particle interactions')

for i in range(3):
    ax[i].plot(rk2[i+1],rk2[i+4], '--', color='r', label='First particle')
    ax[i].plot(rk2[i+7],rk2[i+10], '--', color='b', label='Second Particle')
    ax[i].grid()
    ax[i].set_ylabel(r'x [$\mu m$]')
    ax[i].set_xlabel(r'$v_x$ [$\mu m / s$]')
    ax[i].legend()

fig.tight_layout()
plt.savefig('phase_spaces_with.pdf', dpi=900)
plt.show()

#Plotting

from mpl_toolkits import mplot3d

fig = plt.figure(figsize=(20,10))
plt.suptitle('Particle motions')

ax = fig.add_subplot(1,2,1,projection='3d')

ax.set_title('No particle interactions')
ax.plot(rk2_no[1], rk2_no[2], rk2_no[3], label='Particle 1')
ax.plot(rk2_no[7], rk2_no[8], rk2_no[9], label='Particle 2')
ax.set_xlabel(r'x [$\mu m$]')
ax.set_ylabel(r'y [$\mu m$]')
ax.set_zlabel(r'z [$\mu m$]')
ax.legend()

ax = fig.add_subplot(1,2,2, projection='3d')

ax.set_title('With particle interactions')
ax.plot(rk2[1], rk2[2], rk2[3], label='Particle 1')
ax.plot(rk2[7], rk2[8], rk2[9], label='Particle 2')
ax.set_xlabel(r'x [$\mu m$]')
ax.set_ylabel(r'y [$\mu m$]')
ax.set_zlabel(r'z [$\mu m$]')
ax.legend()

plt.savefig('3d.pdf', dpi=900)
plt.show()


rk_h1 = read('coordinates_rk0.100000.txt')
rk_h05 = read('coordinates_rk0.050000.txt')
rk_h01 = read('coordinates_rk0.010000.txt')
rk_h005 = read('coordinates_rk0.005000.txt')
rk_h001 = read('coordinates_rk0.001000.txt')
eu_h1 = read('coordinates_eu0.100000.txt')
eu_h05 = read('coordinates_eu0.050000.txt')
eu_h01 = read('coordinates_eu0.010000.txt')
eu_h005 = read('coordinates_eu0.005000.txt')
eu_h001 = read('coordinates_eu0.001000.txt')

t = [rk_h1[0],rk_h05[0],rk_h01[0],rk_h005[0],rk_h001[0]]
r_rk = [rk_h1[1],rk_h05[1],rk_h01[1],rk_h005[1],rk_h001[1]]
r_eu = [eu_h1[1],eu_h05[1],eu_h01[1],eu_h005[1],eu_h001[1]]

h = [0.1,0.05,0.01,0.005,0.001]


fig,ax = plt.subplots(1,2, figsize=(20,10))
for i in range(5):
    ax[0].plot(t[i], r_rk[i],label=('h = %.3f' %(h[i])))
ax[0].legend()
ax[0].set_title('Relative error from Runge-Kutta method')
ax[0].set_xlabel(r'Time [$\mu$ s]')
ax[0].set_ylabel(r'Relative error')

for i in range(5):
    ax[1].plot(t[i], r_eu[i],label=('h = %.3f' %(h[i])))
ax[1].legend()
ax[1].set_title('Relative error from Forward Euler method')
ax[1].set_xlabel(r'Time [$\mu$ s]')
ax[1].set_ylabel(r'Relative error')

plt.savefig('error.pdf', dpi=900)

print('RK4')
for i in range(5):
    print(h[i], np.mean(r_rk[i]))

print('Euler')
for i in range(5):
    print(h[i], np.mean(r_eu[i]))


#Calculate r_err
delta_rk = [rk_h1[2],rk_h05[2],rk_h01[2],rk_h005[2],rk_h001[2]]
delta_eu = [eu_h1[2],eu_h05[2],eu_h01[2],eu_h005[2],eu_h001[2]]
r_err_rk = 0
r_err_eu = 0
for i in range(1,5):
    delta_max_i_rk = max(delta_rk[i-1])
    delta_max_i1_rk = max(delta_rk[i])
    r_err_rk += np.log(delta_max_i1_rk/ delta_max_i_rk) / np.log(h[i]/h[i-1])
    delta_max_i_eu = max(delta_eu[i-1])
    delta_max_i1_eu = max(delta_eu[i])
    r_err_eu += np.log(delta_max_i1_eu/ delta_max_i_eu) / np.log(h[i]/h[i-1])
r_err_rk = r_err_rk*1/4
r_err_eu = r_err_eu*1/4

print(r_err_rk)
print(r_err_eu)


res = read('resonance.txt')

fig, ax = plt.subplots()
ax.plot(res[0], res[1], label='f=0.1')
ax.plot(res[0], res[2], label='f=0.4')
ax.plot(res[0], res[3], label='f=0.7')
ax.set_xlabel(r'$\omega_V$ [MHz]')
ax.set_ylabel('Num. particles', rotation=90)
ax.set_title('Numper of particles left in box for a range of frequencies')
ax.legend()
plt.savefig('resonance.pdf', dpi=900)
plt.show()


res = read('resonance_zoomed.txt')

fig, ax = plt.subplots()
ax.plot(res[0], res[1], label='Without particle interactions')
ax.plot(res[0], res[2], label='With particle interactions')
ax.set_xlabel(r'$\omega_V$ [MHz]')
ax.set_ylabel('Num. particles', rotation=90)
ax.set_title('Numper of particles left in box for a range of frequencies')
ax.legend()
plt.savefig('resonance_zoomed.pdf', dpi=900)
plt.show()
