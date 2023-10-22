import numpy as np
import matplotlib.pyplot as plt

LA_Vp = np.genfromtxt('finalVp_1d_LA_dis.txt')
LA_Vs = np.genfromtxt('finalVs_1d_LA_dis.txt')
#Vs_amb_YC = np.genfromtxt('Vs_LA_YC.txt')
#NA_try=np.genfromtxt('aust.tvel')

ak135Vp = np.genfromtxt('finalVp_1d_AU_dis.txt')
ak135Vs = np.genfromtxt('finalVs_1d_AU_dis.txt')
#real= np.genfromtxt('real_itvel.txt')
###
#crst1_Vp = np.genfromtxt('crust1Vp.txt')
#crst1_Vs = np.genfromtxt('crust1Vs.txt')


# depth = LA_Vp[:,0]
# vp_la = LA_Vp[:,1]

# depth_o = our_1d[:,0]

fig = plt.figure(figsize=(5,6),dpi=250)

ax = plt.axes()

# Setting the background color of the plot
# using set_facecolor() method
# ax.set_facecolor("oldlace")
ax.set_facecolor("floralwhite")

# plt.style.use('ggplot')

#plt.plot(crst1_Vp[:,1],crst1_Vp[:,0], 'darkkhaki',label='Crust 1.0')
#plt.plot(crst1_Vs[:,1],crst1_Vs[:,0], 'darkkhaki',linestyle='-.')

#plt.plot(NA_try[:,2],NA_try[:,0], 'black',linestyle='-.')


plt.plot(LA_Vp[:,1],LA_Vp[:,0], 'indianred', label='AuSREM LE')
#plt.plot(AU_Vp[:,1],AU_Vp[:,0], 'darkblue', label='AuSREM 1D Aus')

plt.plot(LA_Vs[:,1],LA_Vs[:,0], 'indianred',linestyle='--')
#plt.plot(Vs_amb_YC[:,1],Vs_amb_YC[:,0], 'steelblue',linestyle='-.',label='ANT LE')

plt.plot(ak135Vp[:,1],ak135Vp[:,0], 'seagreen',label='AuSREM AU')
plt.plot(ak135Vs[:,1],ak135Vs[:,0], 'seagreen',linestyle='--')
######
kw={'markersize':7,'alpha':0.85,'markeredgecolor':'black','markeredgewidth':0.25}

#plt.plot(3.78,25,'o',markerfacecolor='darkkhaki',label='6.54, 3.78 km/s',**kw)
#plt.plot(6.56,25,'o',markerfacecolor='darkkhaki',**kw)

#plt.plot(3.76,25,'o',markerfacecolor='indianred',label='6.45, 3.76 km/s',**kw)
#plt.plot(6.45,25,'o',markerfacecolor='indianred',**kw)

#plt.plot(3.65,25,'o',markerfacecolor='seagreen',label='6.15, 3.65 km/s',**kw)
#plt.plot(6.15,25,'o',markerfacecolor='seagreen',**kw)




plt.xlabel('velocity (km/s)')
plt.ylabel('depth (km)')
plt.grid(color='black', linestyle='-', linewidth=.15,alpha=.15)
plt.gca().invert_yaxis()
plt.xticks(np.arange(1.5, 8, 1))
plt.yticks(np.arange(0, 50, 5))
plt.ylim((47,0))
plt.xlim((2,8.5))
plt.legend(fontsize="8")
plt.savefig('1dmodel_crust_.png',bbox_inches='tight', pad_inches=0.1)
plt.show()
