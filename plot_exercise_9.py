import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Importing from txt file made in the cpp code

# 1.
df = pd.read_csv('coordinates_an.txt', sep=' ', header=None)
x_an = df[1].to_numpy()
y_an = df[2].to_numpy()
z_an = df[3].to_numpy()

df = pd.read_csv('coordinates_euler.txt', sep=' ', header=None)
t = df[0].to_numpy()
x_e = df[1].to_numpy()
y_e = df[2].to_numpy()
z_e = df[3].to_numpy()

df = pd.read_csv('coordinates_rk.txt', sep=' ', header=None)
t = df[0].to_numpy()
x_rk = df[1].to_numpy()
y_rk = df[2].to_numpy()
z_rk = df[3].to_numpy()


fig,ax = plt.subplots(1,2, figsize=(20,10))

fig.suptitle('Particle motion in Penning Trap')

ax[0].plot(x_e,y_e, '--', color='r', label='Euler')
ax[0].plot(x_rk,y_rk, '--', color='b', label='RK4')
ax[0].plot(x_an,y_an, '-', color='k', label='Analytical')
ax[0].plot(x_an[0], y_an[0], 'o', markersize=8, label='Starting position')
ax[0].plot(x_an[-1], y_an[-1], 'o', markersize=8, label='Final position')
ax[0].plot(x_rk[-1], y_rk[-1], 'o', markersize=8, label='Final position_rk')
ax[0].plot(x_e[-1], y_e[-1], 'o', markersize=8, label='Final position_eu')
ax[0].grid()
ax[0].set_ylabel(r'y [$\mu m$]')
ax[0].set_xlabel(r'x [$\mu m$]')
ax[0].legend()
ax[0].axis('equal')

ax[1].plot(t,z_e, '--', color='r', label='Euler')
ax[1].plot(t,z_rk, '--', color='b', label='RK4')
ax[1].plot(t,z_an, '--', color='k', label='Analytical')
ax[1].grid()
ax[1].set_ylabel(r'z [$\mu m$]')
ax[1].set_xlabel(r't [$\mu s$]')
ax[1].legend()

fig.tight_layout()
plt.savefig('1particle_motion.pdf', dpi=900)
plt.show()

# 2.

df = pd.read_csv('rk2particles_without.txt', sep=' ', header=None)
x_rk2 = df[1].to_numpy()
y_rk2 = df[2].to_numpy()
z_rk2 = df[3].to_numpy()
vx_rk2 = df[4].to_numpy()
vy_rk2 = df[5].to_numpy()
vz_rk2 = df[6].to_numpy()
x2_rk2 = df[7].to_numpy()
y2_rk2 = df[8].to_numpy()
z2_rk2 = df[9].to_numpy()
vx2_rk2 = df[10].to_numpy()
vy2_rk2 = df[11].to_numpy()
vz2_rk2 = df[12].to_numpy()

df = pd.read_csv('rk2particles_with.txt', sep=' ', header=None)
x_rk2_w = df[1].to_numpy()
y_rk2_w = df[2].to_numpy()
z_rk2_w = df[3].to_numpy()
vx_rk2_w = df[4].to_numpy()
vy_rk2_w = df[5].to_numpy()
vz_rk2_w = df[6].to_numpy()
x2_rk2_w = df[7].to_numpy()
y2_rk2_w = df[8].to_numpy()
z2_rk2_w = df[9].to_numpy()
vx2_rk2_w = df[10].to_numpy()
vy2_rk2_w = df[11].to_numpy()
vz2_rk2_w = df[12].to_numpy()


fig, ax = plt.subplots(1,2, figsize=(20,10))

fig.suptitle('Particle motions in the xy-plane')

ax[0].plot(x_rk2,y_rk2, '--', color='r', label='First particle')
ax[0].plot(x2_rk2,y2_rk2, '--', color='b', label='Second Particle')
ax[0].grid()
ax[0].set_ylabel(r'x [$\mu m$]')
ax[0].set_xlabel(r'y [$\mu m$]')
ax[0].legend()
ax[0].axis('equal')
ax[0].set_title('Without Particle interactions')

ax[1].plot(x_rk2_w,y_rk2_w, '--', color='r', label='First particle')
ax[1].plot(x2_rk2_w,y_rk2_w, '--', color='b', label='Second Particle')
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

ax[0].plot(x_rk2,vx_rk2, '--', color='r', label='First particle')
ax[0].plot(x2_rk2,vx2_rk2, '--', color='b', label='Second Particle')
ax[0].grid()
ax[0].set_ylabel(r'x [$\mu m$]')
ax[0].set_xlabel(r'$v_x$ [$\mu m / s$]')
ax[0].legend()
#ax[0].axis('equal')

ax[1].plot(y_rk2,vy_rk2, '--', color='r', label='First particle')
ax[1].plot(y2_rk2,vy2_rk2, '--', color='b', label='Second Particle')
ax[1].grid()
ax[1].set_ylabel(r'y [$\mu m$]')
ax[1].set_xlabel(r'$v_y$ [$\mu m / s$]')
ax[1].legend()
#ax[1].axis('equal')

ax[2].plot(z_rk2,vz_rk2, '--', color='r', label='First particle')
ax[2].plot(z2_rk2,vz2_rk2, '--', color='b', label='Second Particle')
ax[2].grid()
ax[2].set_ylabel(r'z [$\mu m$]')
ax[2].set_xlabel(r'$v_z$ [$\mu m / s$]')
ax[2].legend()
#ax[2].axis('equal')

fig.tight_layout()
plt.savefig('phase_spaces_without.pdf', dpi=900)
plt.show()


fig,ax = plt.subplots(1,3, figsize=(30,10))
fig.suptitle('Phase space in Penning Trap with particle interactions')

ax[0].plot(x_rk2_w,vx_rk2_w, '--', color='r', label='First particle')
ax[0].plot(x2_rk2_w,vx2_rk2_w, '--', color='b', label='Second Particle')
ax[0].grid()
ax[0].set_xlabel(r'x [$\mu m$]')
ax[0].set_ylabel(r'$v_x$ [$\mu m / s$]')
ax[0].legend()
#ax[0].axis('equal')

ax[1].plot(y_rk2_w,vy_rk2_w, '--', color='r', label='First particle')
ax[1].plot(y2_rk2_w,vy2_rk2_w, '--', color='b', label='Second Particle')
ax[1].grid()
ax[1].set_xlabel(r'y [$\mu m$]')
ax[1].set_ylabel(r'$v_y$ [$\mu m / s$]')
ax[1].legend()
#ax[1].axis('equal')

ax[2].plot(z_rk2_w,vz_rk2_w, '--', color='r', label='First particle')
ax[2].plot(z2_rk2_w,vz2_rk2_w, '--', color='b', label='Second Particle')
ax[2].grid()
ax[2].set_xlabel(r'z [$\mu m$]')
ax[2].set_ylabel(r'$v_z$ [$\mu m /s$]')
ax[2].legend()
#ax[2].axis('equal')

fig.tight_layout()
plt.savefig('phase_spaces_with.pdf', dpi=900)
plt.show()

#Plotting

from mpl_toolkits import mplot3d

fig = plt.figure(figsize=(20,10))
plt.suptitle('Particle motions')

ax = fig.add_subplot(1,2,1,projection='3d')

ax.set_title('No particle interactions')
ax.plot(x_rk2, y_rk2, z_rk2, label='Particle 1')
ax.plot(x2_rk2, y2_rk2, z2_rk2, label='Particle 2')
ax.set_xlabel(r'x [$\mu m$]')
ax.set_ylabel(r'y [$\mu m$]')
ax.set_zlabel(r'z [$\mu m$]')
ax.legend()

ax = fig.add_subplot(1,2,2, projection='3d')

ax.set_title('With particle interactions')
ax.plot(x_rk2_w, y_rk2_w, z_rk2_w, label='Particle 1')
ax.plot(x2_rk2_w, y2_rk2_w, z2_rk2_w, label='Particle 2')
ax.set_xlabel(r'x [$\mu m$]')
ax.set_ylabel(r'y [$\mu m$]')
ax.set_zlabel(r'z [$\mu m$]')
ax.legend()

plt.savefig('3d.pdf', dpi=900)
plt.show()



df = pd.read_csv('coordinates_rk0.100000.txt', sep=' ', header=None)
t1 = df[0].to_numpy()
rk1 = df[1].to_numpy()

df = pd.read_csv('coordinates_rk0.050000.txt', sep=' ', header=None)
t2 = df[0].to_numpy()
rk2 = df[1].to_numpy()

df = pd.read_csv('coordinates_rk0.010000.txt', sep=' ', header=None)
t3 = df[0].to_numpy()
rk3 = df[1].to_numpy()

df = pd.read_csv('coordinates_rk0.005000.txt', sep=' ', header=None)
t4 = df[0].to_numpy()
rk4 = df[1].to_numpy()

df = pd.read_csv('coordinates_rk0.001000.txt', sep=' ', header=None)
t5 = df[0].to_numpy()
rk5 = df[1].to_numpy()


df = pd.read_csv('coordinates_eu0.100000.txt', sep=' ', header=None)
eu1 = df[1].to_numpy()

df = pd.read_csv('coordinates_eu0.050000.txt', sep=' ', header=None)
eu2 = df[1].to_numpy()

df = pd.read_csv('coordinates_eu0.010000.txt', sep=' ', header=None)
eu3 = df[1].to_numpy()

df = pd.read_csv('coordinates_eu0.005000.txt', sep=' ', header=None)
eu4 = df[1].to_numpy()

df = pd.read_csv('coordinates_eu0.001000.txt', sep=' ', header=None)
eu5 = df[1].to_numpy()



def abs_error(x,y,z):
    return np.sqrt(x**2+y**2+z**2)

h = [0.1,0.05,0.01,0.005,0.001]
t = [t1,t2,t3,t4,t5]

r_rk = [rk1,rk2,rk3,rk4,rk5]

r_eu = [eu1,eu2,eu3,eu4,eu5]

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
