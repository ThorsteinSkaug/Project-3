import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Importing from txt file made in the cpp code
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


df = pd.read_csv('rk2particles_without.txt', sep=' ', header=None)
x_rk2 = df[1].to_numpy()
y_rk2 = df[2].to_numpy()
z_rk2 = df[3].to_numpy()
x2_rk2 = df[4].to_numpy()
y2_rk2 = df[5].to_numpy()
z2_rk2 = df[6].to_numpy()
vx_rk2 = df[7].to_numpy()
vy_rk2 = df[8].to_numpy()
vz_rk2 = df[9].to_numpy()
vx2_rk2 = df[10].to_numpy()
vy2_rk2 = df[11].to_numpy()
vz2_rk2 = df[12].to_numpy()


fig,ax = plt.subplots(1,3, figsize=(30,10))

fig.suptitle('Phase space in Penning Trap without particle interactions')

ax[0].plot(x_rk2,vx_rk2, '--', color='r', label='First particle')
ax[0].plot(x2_rk2,vx2_rk2, '--', color='b', label='Second Particle')
ax[0].grid()
ax[0].set_ylabel(r'x [$\mu m$]')
ax[0].set_xlabel(r'$v_x$ [$\mu m$]')
ax[0].legend()
ax[0].axis('equal')

ax[1].plot(y_rk2,vy_rk2, '--', color='r', label='First particle')
ax[1].plot(y2_rk2,vy2_rk2, '--', color='b', label='Second Particle')
ax[1].grid()
ax[1].set_ylabel(r'y [$\mu m$]')
ax[1].set_xlabel(r'$v_y$ [$\mu m$]')
ax[1].legend()
ax[1].axis('equal')

ax[2].plot(z_rk2,vz_rk2, '--', color='r', label='First particle')
ax[2].plot(z2_rk2,vz2_rk2, '--', color='b', label='Second Particle')
ax[2].grid()
ax[2].set_ylabel(r'z [$\mu m$]')
ax[2].set_xlabel(r'$v_z$ [$\mu m$]')
ax[2].legend()
ax[2].axis('equal')

fig.tight_layout()
plt.savefig('phase_spaces_without.pdf', dpi=900)
plt.show()


df = pd.read_csv('rk2particles_with.txt', sep=' ', header=None)
x_rk2_w = df[1].to_numpy()
y_rk2_w = df[2].to_numpy()
z_rk2_w = df[3].to_numpy()
x2_rk2_w = df[4].to_numpy()
y2_rk2_w = df[5].to_numpy()
z2_rk2_w = df[6].to_numpy()
vx_rk2_w = df[7].to_numpy()
vy_rk2_w = df[8].to_numpy()
vz_rk2_w = df[9].to_numpy()
vx2_rk2_w = df[10].to_numpy()
vy2_rk2_w = df[11].to_numpy()
vz2_rk2_w = df[12].to_numpy()

fig,ax = plt.subplots(1,3, figsize=(30,10))
fig.suptitle('Phase space in Penning Trap with particle interactions')

ax[0].plot(x_rk2_w,vx_rk2_w, '--', color='r', label='First particle')
ax[0].plot(x2_rk2_w,vx2_rk2_w, '--', color='b', label='Second Particle')
ax[0].grid()
ax[0].set_xlabel(r'x [$\mu m$]')
ax[0].set_ylabel(r'$v_x$ [$\mu m$]')
ax[0].legend()
ax[0].axis('equal')

ax[1].plot(y_rk2_w,vy_rk2_w, '--', color='r', label='First particle')
ax[1].plot(y2_rk2_w,vy2_rk2_w, '--', color='b', label='Second Particle')
ax[1].grid()
ax[1].set_xlabel(r'y [$\mu m$]')
ax[1].set_ylabel(r'$v_y$ [$\mu m$]')
ax[1].legend()
ax[1].axis('equal')

ax[2].plot(z_rk2_w,vz_rk2_w, '--', color='r', label='First particle')
ax[2].plot(z2_rk2_w,vz2_rk2_w, '--', color='b', label='Second Particle')
ax[2].grid()
ax[2].set_xlabel(r'z [$\mu m$]')
ax[2].set_ylabel(r'$v_z$ [$\mu m$]')
ax[2].legend()
ax[2].axis('equal')

fig.tight_layout()
plt.savefig('phase_spaces_without.pdf', dpi=900)
plt.show()

#Plotting

from mpl_toolkits import mplot3d

fig = plt.subplots(1,2, figsize=(20,10))
ax = plt.axes(projection='3d')

plt.suptitle('Particle motions')

ax[0].set_title('No particle interactions')
ax[0].plot3D(x_rk2, y_rk2, z_rk2, label='Particle 1')
ax[0].plot3D(x2_rk2, y2_rk2, z2_rk2, label='Particle 2')
ax[0].set_xlabel(r'x [$\mu m$]')
ax[0].set_ylabel(r'y [$\mu m$]')
ax[0].set_zlabel(r'z [$\mu m$]')

ax[1].set_title('With particle interactions')
ax[1].plot3D(x_rk2_w, y_rk2_w, z_rk2_w, label='Particle 1')
ax[1].plot3D(x2_rk2_w, y2_rk2_w, z2_rk2_w, label='Particle 2')
ax[1].set_xlabel(r'x [$\mu m$]')
ax[1].set_ylabel(r'y [$\mu m$]')
ax[1].set_zlabel(r'z [$\mu m$]')


"""
df = pd.read_csv('coordinates_rk0.1.txt', sep=' ', header=None)
t = df[0].to_numpy()
x_error1 = df[1].to_numpy()
y_error1 = df[2].to_numpy()
z_error1 = df[3].to_numpy()

df = pd.read_csv('coordinates_rk0.05.txt', sep=' ', header=None)
t = df[0].to_numpy()
x_error2 = df[1].to_numpy()
y_error2 = df[2].to_numpy()
z_error2 = df[3].to_numpy()

df = pd.read_csv('coordinates_rk0.01.txt', sep=' ', header=None)
t = df[0].to_numpy()
x_error3 = df[1].to_numpy()
y_error3 = df[2].to_numpy()
z_error3 = df[3].to_numpy()

df = pd.read_csv('coordinates_rk0.005.txt', sep=' ', header=None)
t = df[0].to_numpy()
x_error4 = df[1].to_numpy()
y_error4 = df[2].to_numpy()
z_error4 = df[3].to_numpy()

df = pd.read_csv('coordinates_rk0.001.txt', sep=' ', header=None)
t = df[0].to_numpy()
x_error5 = df[1].to_numpy()
y_error5 = df[2].to_numpy()
z_error5 = df[3].to_numpy()


plt.plot(t,)
"""
