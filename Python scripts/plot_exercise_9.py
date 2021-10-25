import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib


#Set the size of the plots, text, and axes
small= 17
medium = 24
big = 30

plt.rc('font', size=small)          # controls default text sizes
plt.rc('axes', titlesize=medium)     # fontsize of the axes title
plt.rc('axes', labelsize=medium)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=medium)    # fontsize of the tick labels
plt.rc('ytick', labelsize=medium)    # fontsize of the tick labels
plt.rc('legend', fontsize=small)    # legend fontsize
plt.rc('figure', titlesize=big)  # fontsize of the figure title



#This function takes a filename, reads it and output the columns as a list of np.arrays
def read(filename):
    df = pd.read_csv(filename, sep=' ', header=None)
    var = []
    for i in range(len(df.columns)):
        var.append(df[i].to_numpy())
    return var

# 1.

#Read the files for the single particle cases
pos_an =  read('coordinates_an.txt')
t001 = pos_an[0]
pos_eu = read('coordinates_euler.txt')
pos_rk = read('coordinates_rk.txt')


#Plot the results for the single particle cases
fig,ax = plt.subplots(figsize=(10,10))

fig.suptitle('Particle motion in Penning Trap')

ax.plot(pos_eu[1],pos_eu[2], '--', color='r', label='Euler')
ax.plot(pos_rk[1],pos_rk[2], '--', color='b', label='RK4')
ax.plot(pos_an[1],pos_an[2], '-', color='k', label='Analytical')
ax.plot(pos_an[1][0], pos_an[2][0], 'o', markersize=8, label='Starting position')
ax.plot(pos_an[1][-1], pos_an[2][-1], 'o', markersize=8, label='Final position')
ax.plot(pos_rk[1][-1], pos_rk[2][-1], 'o', markersize=8, label='Final position_rk')
ax.plot(pos_eu[1][-1], pos_eu[2][-1], 'o', markersize=8, label='Final position_eu')
ax.grid()
ax.set_ylabel(r'y [$\mu m$]')
ax.set_xlabel(r'x [$\mu m$]')
ax.legend(loc=1)
ax.axis('equal')
fig.tight_layout()
plt.savefig('1particle_motion1.pdf', dpi=1200)
plt.show()

fig,ax = plt.subplots(figsize=(10,10))

ax.plot(t001,pos_eu[3], '--', color='r', label='Euler')
ax.plot(t001,pos_rk[3], '--', color='b', label='RK4')
ax.plot(t001,pos_an[3], '--', color='k', label='Analytical')
ax.grid()
ax.set_ylabel(r'z [$\mu m$]')
ax.set_xlabel(r't [$\mu s$]')
ax.legend()
fig.suptitle(r'Motion in $z$-direction')
fig.tight_layout()
plt.savefig('1particle_motion2.pdf', dpi=1200)
plt.show()


# 2.
#Read the files for the two particle cases
rk2_no =  read('rk2particles_without.txt')
rk2 =  read('rk2particles_with.txt')

#Plot the results for the two particle cases
fig, ax = plt.subplots(figsize=(10,10))

ax.plot(rk2_no[1],rk2_no[2], '--', color='r', label='First particle')
ax.plot(rk2_no[7],rk2_no[8], '--', color='b', label='Second Particle')
ax.set_ylabel(r'x [$\mu m$]')
ax.set_xlabel(r'y [$\mu m$]')
ax.legend()
ax.axis('equal')
ax.grid()
plt.suptitle('Particle motions without particle interactions')

fig.tight_layout()
plt.savefig('2particle_motion1.pdf', dpi=1200)
plt.show()

fig, ax = plt.subplots(figsize=(10,10))

ax.plot(rk2[1],rk2[2], '--', color='r', label='First particle')
ax.plot(rk2[7],rk2[8], '--', color='b', label='Second Particle')
ax.set_xlabel(r'x [$\mu m$]')
ax.set_ylabel(r'y [$\mu m$]')
ax.legend()
ax.grid()
ax.axis('equal')
plt.suptitle('Particle motions with particle interactions')

fig.tight_layout()
plt.savefig('2particle_motion2.pdf', dpi=1200)
plt.show()

# 3.

#Plot the phase space
fig,ax = plt.subplots(1,3, figsize=(30,10))

fig.suptitle('Phase space in Penning Trap without particle interactions')

for i in range(3):
    ax[i].plot(rk2_no[i+1],rk2_no[i+4], '--', color='r', label='First particle')
    ax[i].plot(rk2_no[i+7],rk2_no[i+10], '--', color='b', label='Second Particle')
    ax[i].grid()
    ax[i].legend()

ax[0].set_xlabel(r'x [$\mu m$]')
ax[0].set_ylabel(r'$v_x$ [$\mu m / \mu s$]')
ax[1].set_xlabel(r'y [$\mu m$]')
ax[1].set_ylabel(r'$v_y$ [$\mu m / \mu s$]')
ax[2].set_xlabel(r'z [$\mu m$]')
ax[2].set_ylabel(r'$v_z$ [$\mu m / \mu s$]')
fig.tight_layout()
plt.savefig('phase_spaces_without.pdf', dpi=1200)
plt.show()

fig,ax = plt.subplots(1,3, figsize=(30,10))

fig.suptitle('Phase space in Penning Trap with particle interactions')

for i in range(3):
    ax[i].plot(rk2[i+1],rk2[i+4], '--', color='r', label='First particle')
    ax[i].plot(rk2[i+7],rk2[i+10], '--', color='b', label='Second Particle')
    ax[i].grid()
    ax[i].legend()

ax[0].set_xlabel(r'x [$\mu m$]')
ax[0].set_ylabel(r'$v_x$ [$\mu m / \mu s$]')
ax[1].set_xlabel(r'y [$\mu m$]')
ax[1].set_ylabel(r'$v_y$ [$\mu m / \mu s$]')
ax[2].set_xlabel(r'z [$\mu m$]')
ax[2].set_ylabel(r'$v_z$ [$\mu m / \mu s$]')
fig.tight_layout()

plt.savefig('phase_spaces_with.pdf', dpi=1200)
plt.show()


#Plotting the 3D plots
from mpl_toolkits import mplot3d
plt.rc('xtick', labelsize=14)    # fontsize of the tick labels
plt.rc('ytick', labelsize=14)
plt.rc('axes', titlesize=21)     # fontsize of the axes title
plt.rc('axes', labelsize=21)

fig = plt.figure(figsize=(10,10))
plt.suptitle('Particle motions without interactions')
ax = plt.axes(projection='3d')

ax.plot(rk2_no[1], rk2_no[2], rk2_no[3], label='Particle 1')
ax.plot(rk2_no[7], rk2_no[8], rk2_no[9], label='Particle 2')
ax.set_xlabel(r'x [$\mu m$]')
ax.set_ylabel(r'y [$\mu m$]')
ax.set_zlabel(r'z [$\mu m$]')
ax.legend()

plt.tight_layout()
plt.savefig('3d1.pdf', dpi=1200)
plt.show()

fig = plt.figure(figsize=(10,10))
plt.suptitle('Particle motions with interactions')
ax = plt.axes(projection='3d')

ax.plot(rk2[1], rk2[2], rk2[3], label='Particle 1')
ax.plot(rk2[7], rk2[8], rk2[9], label='Particle 2')
ax.set_xlabel(r'x [$\mu m$]')
ax.set_ylabel(r'y [$\mu m$]')
ax.set_zlabel(r'z [$\mu m$]')
ax.legend()

plt.tight_layout()
plt.savefig('3d2.pdf', dpi=1200)
plt.show()

plt.rc('axes', titlesize=medium)     # fontsize of the axes title
plt.rc('axes', labelsize=medium)
plt.rc('xtick', labelsize=small)    # fontsize of the tick labels
plt.rc('ytick', labelsize=small)

#Read the error files
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

t = [rk_h1[0],rk_h05[0],rk_h01[0],rk_h005[0],rk_h001[0]] #List storing the different time arrays
r_rk = [rk_h1[1],rk_h05[1],rk_h01[1],rk_h005[1],rk_h001[1]] #List storing the results from runge-kutta
r_eu = [eu_h1[1],eu_h05[1],eu_h01[1],eu_h005[1],eu_h001[1]] #List storing the results from forward euler

h = [0.1,0.05,0.01,0.005,0.001] #List storing the different step size values

#Plot the error plots
fig,ax = plt.subplots(figsize=(10,10))
for i in range(5):
    ax.plot(t[i], r_rk[i],label=('h = %.3f' %(h[i])))
ax.legend()
plt.suptitle('Relative error from Runge-Kutta method')
ax.set_xlabel(r'Time [$\mu$ s]')
ax.set_ylabel(r'log(Relative error)')
ax.set_yscale('log')
ax.grid()
plt.tight_layout()
plt.savefig('logerror1.pdf', dpi=1200)


fig,ax = plt.subplots(figsize=(10,10))

for i in range(5):
    ax.plot(t[i], r_eu[i],label=('h = %.3f' %(h[i])))
ax.legend()
plt.suptitle('Relative error from Forward Euler method')
ax.set_xlabel(r'Time [$\mu$ s]')
ax.set_ylabel(r'log (Relative error)')
ax.set_yscale('log')
ax.grid()
plt.tight_layout()
plt.savefig('logerror2.pdf', dpi=1200)

fig,ax = plt.subplots(figsize=(10,10))
for i in range(5):
    ax.plot(t[i], r_rk[i],label=('h = %.3f' %(h[i])))
ax.legend()
plt.suptitle('Relative error from Runge-Kutta method')
ax.set_xlabel(r'Time [$\mu$ s]')
ax.set_ylabel(r'log(Relative error)')
ax.grid()
plt.tight_layout()
plt.savefig('error1.pdf', dpi=1200)

fig,ax = plt.subplots(figsize=(10,10))
for i in range(5):
    ax.plot(t[i], r_eu[i],label=('h = %.3f' %(h[i])))
ax.legend()
plt.suptitle('Relative error from Forward Euler method')
ax.set_xlabel(r'Time [$\mu$ s]')
ax.set_ylabel(r'log (Relative error)')
ax.grid()
plt.tight_layout()
plt.savefig('error2.pdf', dpi=1200)

#Print the mean error values for runge kutta
print('RK4')
for i in range(5):
    print(h[i], np.mean(r_rk[i]))

#Print the mean error values for forward euler
print('Euler')
for i in range(5):
    print(h[i], np.mean(r_eu[i]))


#Calculate r_err
delta_rk = [rk_h1[2],rk_h05[2],rk_h01[2],rk_h005[2],rk_h001[2]]
delta_eu = [eu_h1[2],eu_h05[2],eu_h01[2],eu_h005[2],eu_h001[2]]
r_err_rk = 0 #Storing the r_err for runge kutta
r_err_eu = 0 #Storing the r_err for forward euler
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


#Read and plot the resonance and resonance_zoom results
res = read('resonance.txt')

fig, ax = plt.subplots(figsize=(10,10))
ax.plot(res[0], res[1]/100, label='f=0.1')
ax.plot(res[0], res[2]/100, label='f=0.4')
ax.plot(res[0], res[3]/100, label='f=0.7')
ax.set_xlabel(r'$\omega_V$ [MHz]')
ax.set_ylabel('Num. particles', rotation=90)
plt.suptitle('Number of particles left in the box')
ax.legend()
ax.grid()
plt.tight_layout()
plt.savefig('resonance.pdf', dpi=1200)
plt.show()


res = read('resonance_zoom.txt')

fig, ax = plt.subplots(figsize=(10,10))
ax.plot(res[0], res[1]/100, label='Without particle interactions')
ax.plot(res[0], res[2]/100, label='With particle interactions')
ax.set_xlabel(r'$\omega_V$ [MHz]')
ax.set_ylabel('Num. particles', rotation=90)
plt.suptitle('Number of particles left in the box')
ax.legend()
ax.grid()
plt.tight_layout()
plt.savefig('resonance_zoom.pdf', dpi=1200)
plt.show()
