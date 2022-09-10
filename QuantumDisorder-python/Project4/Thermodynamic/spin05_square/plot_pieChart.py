import matplotlib.pyplot as plt
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"]})
import numpy as np
# Pie chart, where the slices will be ordered and plotted counter-clockwise:
mylabels = ['Strip', 'FM', 'AF']
mycolors = ["#EFE80E", "#F59407", "#2ECC71"]

sizes = np.array([33.3, 33.3, 33.3])
explode = (0, 0, 0)  # only "explode" the 2nd slice (i.e. 'Hogs')

fig, ax = plt.subplots(figsize=(6,6))
#ax1.pie(sizes, explode=explode, colors=mycolors, labels=mylabels, autopct='%1.1f%%',
#        shadow=True, startangle=30)
ax.pie(sizes, explode=explode, colors=mycolors,
        startangle=30,
        wedgeprops={"edgecolor":"k",'linewidth': 0.5, 'linestyle': 'solid', 'antialiased': True})
#wedgeprops={"edgecolor":"k",'linewidth': 2, 'linestyle': 'dashed', 'antialiased': True}

ax.text(-0.2,0.5,'Stripe',fontsize=24)
ax.text(0.3,-0.4,'AFM',fontsize=24)
ax.text(-0.6,-0.4,'FM',fontsize=24)

ax.text(-0.15,1.05,r'$J_2$',fontsize=18)
ax.text(1.0,-0.2,r'$J_1$',fontsize=18)

#---------------------
#:1
x,y=np.cos(np.deg2rad(0)), np.sin(np.deg2rad(0))
a_circle = plt.Circle((x,y), 0.02,color='r')
ax.add_artist(a_circle)
ax.text(x+0.05,y-0.07,r'Pb$_2$Cu(OH)$_4$Cl$_2$',color='k',fontsize=10)

#:2
x,y=np.cos(np.deg2rad(0)), np.sin(np.deg2rad(0))
a_circle = plt.Circle((x,y), 0.02,color='r')
ax.add_artist(a_circle)
ax.text(x+0.05,y+0.03,r'Cu(pz)$_2$(ClO$_4$)$_2$',color='k',fontsize=10)


#:3
x,y=np.cos(np.deg2rad(15)), np.sin(np.deg2rad(15))
a_circle = plt.Circle((x,y), 0.02, color='r')
ax.add_artist(a_circle)
ax.text(x+0.05,y-0.03,r'VOMoO$_4$',color='k',fontsize=10)

#:4
x,y=np.cos(np.deg2rad(18)), np.sin(np.deg2rad(18))
a_circle = plt.Circle((x,y), 0.02, color='r')
ax.add_artist(a_circle)
ax.text(x+0.05,y-0.0,r'PbVO$_3$',color='k',fontsize=10)



#:5
x,y=np.cos(np.deg2rad(55)), np.sin(np.deg2rad(55))
a_circle = plt.Circle((x,y), 0.02, color='r')
ax.add_artist(a_circle)
ax.text(x+0.1,y-0.06,r'Sr$_2$CuMoO$_6$',color='k',fontsize=10)


#:6
x,y=np.cos(np.deg2rad(65)), np.sin(np.deg2rad(65))
a_circle = plt.Circle((x,y), 0.02, color='r')
ax.add_artist(a_circle)
ax.text(x+0.05,y-0.02,r'Li$_2$VOGeO$_4$',color='k',fontsize=10)



#:7
x,y=np.cos(np.deg2rad(72)), np.sin(np.deg2rad(72))
a_circle = plt.Circle((x,y), 0.02, color='r')
ax.add_artist(a_circle)
ax.text(x+0.06,y+0.03,r'Sr$_2$CuWO$_6$',color='k',fontsize=10)




#:8
x,y=np.cos(np.deg2rad(85)), np.sin(np.deg2rad(85))
a_circle = plt.Circle((x,y), 0.02, color='r')
ax.add_artist(a_circle)
ax.text(x+0.035,y+0.045,r'Li$_2$VOSiO$_4$',color='k',fontsize=10)


#:9
x,y=np.cos(np.deg2rad(110)), np.sin(np.deg2rad(110))
a_circle = plt.Circle((x,y), 0.02, color='r')
ax.add_artist(a_circle)
ax.text(x-0.5,y+0.045,r'Na$_{1.5}$VOPO$_4$F$_{0.5}$',color='k',fontsize=10)


#:10
x,y=np.cos(np.deg2rad(120)), np.sin(np.deg2rad(120))
a_circle = plt.Circle((x,y), 0.02, color='r')
ax.add_artist(a_circle)
ax.text(x-0.5,y+0.035,r'Pb$_2$VO(PO$_4$)$_2$',color='k',fontsize=10)



#:11
x,y=np.cos(np.deg2rad(130)), np.sin(np.deg2rad(130))
a_circle = plt.Circle((x,y), 0.02, color='r')
ax.add_artist(a_circle)
ax.text(x-0.55,y+0.035,r'SrZnVO(PO$_4$)$_2$',color='k',fontsize=10)


#:12
x,y=np.cos(np.deg2rad(140)), np.sin(np.deg2rad(140))
a_circle = plt.Circle((x,y), 0.02, color='r')
ax.add_artist(a_circle)
ax.text(x-0.55,y+0.035,r'BaCdVO(PO$_4$)$_2$',color='k',fontsize=10)


#:13
x,y=np.cos(np.deg2rad(160)), np.sin(np.deg2rad(160))
a_circle = plt.Circle((x,y), 0.02, color='r')
ax.add_artist(a_circle)
ax.text(x-0.55,y+0.035,r'(CuCl)LaNb$_2$O$_7$',color='k',fontsize=10)

#------------------------



ax.quiver([-1, 0], [0, -1], [2.15, 0], [0, 2.15], width=0.0025,  linestyle= 'dashed',
            color=['k','k'], angles='xy', scale_units='xy', scale=1)
#ax.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

plt.savefig("pie.pdf",
               bbox_inches='tight',
               pad_inches=0)
plt.show()
