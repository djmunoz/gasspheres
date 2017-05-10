import numpy as np
import matplotlib.pyplot as plt
from pylab import gca
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import brentq
from scipy.optimize import fsolve
import matplotlib.cm as cm
from mpl_toolkits.axes_grid import make_axes_locatable 
import  matplotlib.axes as maxes
import sys
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import FancyArrowPatch

#############################################################################


#def plot_orbit(radius,incl,node):

def contours(phi,theta,spin_phi,spin_theta):
    ndotnspin = np.sin(spin_theta * np.pi/180) * np.cos(spin_phi * np.pi/180) * \
                np.sin(theta * np.pi/180) * np.cos(phi * np.pi/180) + \
                np.sin(spin_theta * np.pi/180) * np.sin(spin_phi * np.pi/180) * \
                np.sin(theta * np.pi/180) * np.sin(phi * np.pi/180) + \
                np.cos(spin_theta * np.pi/180) * np.cos(theta * np.pi/180)
    
    return  1.0/3 - np.arccos(ndotnspin)


if __name__ == "__main__":

    view_azim, view_lat,view_dist = 60, 20, 7.97

    spin_lambda, spin_i = 25, 80

    spin_theta= 180/np.pi * np.arccos(np.sin(spin_i * np.pi/180) * np.cos(spin_lambda * np.pi/180))
    spin_phi = 180/np.pi * np.arctan2(np.sin(spin_i * np.pi/180) * np.sin(spin_lambda * np.pi/180),np.cos(spin_i * np.pi/180))
    spin_view_azim, spin_view_lat = view_azim, view_lat

    print spin_theta,spin_phi

    
    fig = plt.figure(1, figsize=(4.0,4.0))
    fig.subplots_adjust(top=0.90,bottom=0.03,right=0.93,left=0.04, hspace=0.3,wspace=0.4)


    # Planet orbit
    ax3 = fig.add_subplot(111,projection='3d')
    fig.subplots_adjust(top=0.8,bottom=0.10,right=0.85,left=0.15)
    fig.subplots_adjust(top=1.0,bottom=0.0,right=1.0,left=0.0)
    ax3.set_axis_bgcolor('none')
    ax3.axis('off')
    phases = np.linspace(0,2*np.pi,200)
    radius = 1.4
    xorb,yorb,zorb = radius*np.cos(phases), radius*np.sin(phases),radius*np.zeros(phases.shape[0])    
    ax3.view_init(view_lat, view_azim)
    ax3.dist = view_dist
    ax3.plot(xorb,yorb,zorb,zorder=100,color='firebrick')
    ax3.set_xlim(-1,1.0)
    ax3.set_ylim(-1,1.0)
    ax3.set_zlim(-1,1.0)
    X03,X13 = ax3.get_xlim()
    Y03,Y13 = ax3.get_ylim()
    Z03,Z13 = ax3.get_zlim()


    # Arrows for coordinate triad
    arrow_begin,arrow_end = [0.7,-1.1,0.9],[0.95,-1.1,0.9]
    ax3.plot([arrow_begin[0],arrow_end[0]],[arrow_begin[1],arrow_end[1]],[arrow_begin[2],arrow_end[2]],
             color='k',lw=2.0)
    ax3.text(arrow_end[0],1.05*arrow_end[1],1.02*arrow_end[2],r'$\hat{\mathbf{x}}$')
    arrow_begin,arrow_end = [0.7,-1.1,0.9],[0.7,-0.75,0.9]
    ax3.plot([arrow_begin[0],arrow_end[0]],[arrow_begin[1],arrow_end[1]],[arrow_begin[2],arrow_end[2]],
             color='k',lw=2.0)
    ax3.text(arrow_end[0],0.9*arrow_end[1],arrow_end[2],r'$\hat{\mathbf{y}}$')
    arrow_begin,arrow_end = [0.7,-1.1,0.9],[0.7,-1.1,1.15]
    ax3.plot([arrow_begin[0],arrow_end[0]],[arrow_begin[1],arrow_end[1]],[arrow_begin[2],arrow_end[2]],
             color='k',lw=2.0)
    ax3.text(arrow_end[0],0.9*arrow_end[1],1.01*arrow_end[2],r'$\hat{\mathbf{z}}$')

    
    # Sphere
    ax1 = fig.add_axes([0.167, 0.162, 0.7, 0.7])
    
    lat_0,lon_0 = spin_view_lat, spin_view_azim
    map = Basemap(projection='ortho',lat_0=lat_0,lon_0=lon_0,resolution='l')
    #map.drawmeridians(np.arange(0,360,20), latmax=90,alpha=0.5,color='gray',linewidth=0.5,dashes=[3,1])
    #map.drawparallels(np.arange(-80,80,10),alpha=0.5,color='gray',linewidth=0.5,dashes=[3,1])
    
    x,y=map(np.linspace(-20,180,30), np.repeat(0,30))
    #map.plot(x,y,color='gray', lw=3.0,alpha=0.3)
    x,y=map(np.repeat(0,30),np.linspace(-90,90,30))
    #map.plot(x,y,color='gray', lw=3.0,alpha=0.3)
    
    
    nvals_y = 200; nvals_x = 200; deltax = 360.0/(nvals_x-1)
    yvals = (90.0-deltax*np.indices((nvals_y,nvals_x))[0,:,:])
    yvals[np.abs(yvals) < 5.0e-2] = 0
    xvals = (deltax*np.indices((nvals_y,nvals_x))[1,:,:])
    x, y = map(xvals, yvals)
    cont = np.vectorize(contours)(xvals,90.0-yvals,spin_phi,spin_theta)
    levels=np.linspace(cont.min()*0.995,cont.max()*0.995,30)
    cs = map.contour(x,y,cont,levels=levels,linewidths=0.5,colors='gray',linestyles='-')
    
    # And the arcs subtended by the observable angles
    map.drawgreatcircle(spin_phi, 90-spin_theta,0,0,del_s=50,color='orangered', lw=1.8)
    map.drawgreatcircle(90, 90-spin_theta*np.sin(spin_phi*np.pi/180),0,90,del_s=50,color='dodgerblue', lw=1.8)
    

    
    # Pole axis
    #x,y=map(spin_phi, 90-spin_theta)
    #llx2, lly2 = ax1.transAxes.inverted().transform(ax1.transData.transform((x,y)))
    #arrowx_begin,arrowx_end = llx2,llx2
    #arrowy_begin,arrowy_end = lly2,lly2+0.20
    #ax2 = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    #ax2.set_axis_bgcolor('none')
    #ax2.axis('off')
    #X0,X1 = ax2.get_xlim()
    #Y0,Y1 = ax2.get_ylim()
    #ax2.plot([arrowx_begin,arrowx_end],[arrowy_begin,arrowy_end],'k',transform=ax1.transAxes,lw=3.0)
    #ax2.plot([arrowx_end],[arrowy_end],'k',transform=ax1.transAxes,marker=[3,0,0],ms=10.0)
    #ax2.text(arrowx_end+0.03,arrowy_end-0.08,r'${\bf{S}}$',transform=ax1.transAxes,size=22)


    # 3-D plotting
    ax3 = fig.add_subplot(111,projection='3d')
    fig.subplots_adjust(top=1,bottom=0.,right=1.0,left=0.0)
    ax3.set_axis_bgcolor('none')

    # Pole axis
    arrow_begin,arrow_end = [np.sin(spin_theta*np.pi/180) * np.cos(spin_phi*np.pi/180),
                             np.sin(spin_theta*np.pi/180) * np.sin(spin_phi*np.pi/180),
                             np.cos(spin_theta*np.pi/180)],\
                             [1.5 * np.sin(spin_theta*np.pi/180) * np.cos(spin_phi*np.pi/180),
                              1.5 * np.sin(spin_theta*np.pi/180) * np.sin(spin_phi*np.pi/180),
                              1.5 * np.cos(spin_theta*np.pi/180)]

    ax3.plot([arrow_begin[0]],[arrow_begin[1]],[arrow_begin[2]],'o',color='gray',mew=0,markersize=5.0)
    ax3.plot([arrow_begin[0],arrow_end[0]],[arrow_begin[1],arrow_end[1]],[arrow_begin[2],arrow_end[2]],
             color='k',lw=3.0)
    ax3.plot([arrow_end[0]],[arrow_end[1]],[arrow_end[2]],'k',marker=[3,0,-spin_theta * np.cos(spin_phi*np.pi/180.0)],ms=10.0)
    ax3.text(arrow_end[0]+0.03,arrow_end[1]+0.2,arrow_end[2]-0.05,r'${\bf{S}}$',size=22)
    
    # More coordinate axes
    arrow_begin,arrow_end = [1,0,0],[1.7,0,0]
    ax3.plot([arrow_begin[0],arrow_end[0]],[arrow_begin[1],arrow_end[1]],[arrow_begin[2],arrow_end[2]],
             color='k',lw=3.0,alpha=0.4)
    arrow_begin,arrow_end = [0,1,0],[0,2.3,0]
    ax3.plot([arrow_begin[0],arrow_end[0]],[arrow_begin[1],arrow_end[1]],[arrow_begin[2],arrow_end[2]],
             color='k',lw=3.0,alpha=0.4)
    arrow_begin,arrow_end = [0,0,1],[0,0,1.7]
    ax3.plot([arrow_begin[0],arrow_end[0]],[arrow_begin[1],arrow_end[1]],[arrow_begin[2],arrow_end[2]],
             color='k',lw=3.0,alpha=0.4)
    angles_theta,angles_phi = np.linspace(0,np.pi/2,20), np.repeat(0,20)
    ax3.plot(np.cos(angles_phi) * np.sin(angles_theta),np.sin(angles_phi) * np.sin(angles_theta),np.cos(angles_theta),
             color='k',lw=1.5,alpha=0.4,ls='--')
    angles_theta,angles_phi = np.linspace(0,np.pi/2,20), np.repeat(np.pi/2,20)
    ax3.plot(np.cos(angles_phi) * np.sin(angles_theta),np.sin(angles_phi) * np.sin(angles_theta),np.cos(angles_theta),
            color='k',lw=1.5,alpha=0.4,ls='--')
    angles_theta,angles_phi = np.repeat(np.pi/2,20), np.linspace(0,np.pi/2,20)
    ax3.plot(np.cos(angles_phi) * np.sin(angles_theta),np.sin(angles_phi) * np.sin(angles_theta),np.cos(angles_theta),
             color='k',lw=1.5,alpha=0.4,ls='--')

    # And the observable angles
    #angles_theta,angles_phi = np.linspace(spin_theta*np.pi/180.0,np.pi/2,20), np.linspace(spin_phi*np.pi/180.0,0,20)
    #ax3.plot(np.cos(angles_phi) * np.sin(angles_theta),np.sin(angles_phi) * np.sin(angles_theta),np.cos(angles_theta),
    #         color='r',lw=1.5,ls='-')
    
    
    # Orbit again   
    ax3.axis('off')
    ax3.view_init(view_lat, view_azim)
    ax3.dist = view_dist
    ind = (phases > 0) & (phases < np.pi)
    xorb, yorb, zorb =  xorb[ind],yorb[ind],zorb[ind]
    ax3.plot(xorb,yorb,zorb,zorder=0,color='firebrick')
    ax3.set_xlim(X03,X13)
    ax3.set_ylim(Y03,Y13)
    ax3.set_zlim(Z03,Z13)


    
    
    
    fig.savefig("sphere_diagram.pdf")
    fig.clf()
    
