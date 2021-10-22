import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Importing from txt file made in the cpp code
df = pd.read_csv('coordinates_euler.txt', sep=' ', header=None)
t = df[0].to_numpy()
x_e = df[1].to_numpy()
y_e = df[2].to_numpy()
z_e = df[3].to_numpy()


df = pd.read_csv('coordinates_euler.txt', sep=' ', header=None)
t = df[0].to_numpy()
x_rk = df[1].to_numpy()
y_rk = df[2].to_numpy()
z_rk = df[3].to_numpy()


df = pd.read_csv('coordinates_rk2particles_with.txt', sep=' ', header=None)
x_rk2 = df[1].to_numpy()
y_rk2 = df[2].to_numpy()
z_rk2 = df[3].to_numpy()
x2_rk2 = df[4].to_numpy()
y2_rk2 = df[5].to_numpy()
z2_rk2 = df[6].to_numpy()

df = pd.read_csv('coordinates_an.txt', sep=' ', header=None)
x_an = df[1].to_numpy()
y_an = df[2].to_numpy()
z_an = df[3].to_numpy()


#Plotting

fig,ax = plt.subplots(1,2, figsize=(20,10))

fig.suptitle('Paticle motion in Penning Trap')

ax[0].plot(x_e,y_e, '--', color='r', label='Euler')
ax[0].plot(x_rk,y_rk, '--', color='b', label='RK4')
#ax[0].plot(x_an,y_an, '-', color='k', label='Analytical
ax[0].grid()
ax[0].set_ylabel(r'y [$\mu m$]')
ax[0].set_xlabel(r'x [$\mu m$]')
ax[0].legend()
ax[0].axis('equal')

ax[1].plot(t,z_e, '--', color='r', label='Euler')
ax[1].plot(t,z_rk, '--', color='b', label='RK4')
#ax[1].plot(t,z_an, '-', color='k', label='Analytical')
ax[1].grid()
ax[1].set_ylabel(r'z [$\mu m$]')
ax[1].set_xlabel(r't [$\mu s$]')
ax[1].legend()

fig.tight_layout()
plt.savefig('coordinates.pdf', dpi=900)
plt.show()
