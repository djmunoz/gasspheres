import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import healpy as hp
from scipy.spatial import SphericalVoronoi
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from mpl_toolkits.mplot3d.art3d import Poly3DCollection




if __name__ == "__main__":

    #fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')
    
    #u = np.linspace(0, 2 * np.pi, 100)
    #v = np.linspace(0, np.pi, 100)
    
    #x = 10 * np.outer(np.cos(u), np.sin(v))
    #y = 10 * np.outer(np.sin(u), np.sin(v))
    #z = 10 * np.outer(np.ones(np.size(u)), np.cos(v))
    #z1 = z * np.cos(0.5*x)
    
    #N = z1 / z1.max()  # normalize 0..1
    #surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, facecolors=cm.jet(N), linewidth=0, antialiased=False, shade=False)
    #ax.set_axis_off()
    #m = cm.ScalarMappable(cmap=cm.jet)
    #m.set_array(z1)
    #plt.colorbar(m)
    #plt.show() 
    #exit()

    view_azim, view_lat,view_dist = 60, 20, 7.97
    spin_view_azim, spin_view_lat = view_azim, view_lat

    lat_0,lon_0 = spin_view_lat, spin_view_azim
    #m = Basemap(projection='ortho',lat_0=lat_0,lon_0=lon_0,resolution='l')
    m = Basemap(projection='hammer',lon_0=180)
    m.drawmeridians(np.arange(0,360,20),latmax=90,alpha=0.1,color='gray',linewidth=1.0,dashes=[10000,0.01])
    m.drawparallels(np.arange(-80,80,10),alpha=0.1,color='gray',linewidth=1.0,dashes=[10000,0.01])
    

    NSIDE = 4
    npix = hp.nside2npix(NSIDE)
    pixelcoords = hp.pix2ang(NSIDE,np.arange(npix),nest=True)
    xpix, ypix, zpix = np.cos(pixelcoords[1][:]) * np.sin(pixelcoords[0][:]), \
                       np.sin(pixelcoords[1][:]) * np.sin(pixelcoords[0][:]), \
                       np.cos(pixelcoords[0][:]) 

    x,y = m(pixelcoords[1]* 180.0/np.pi, pixelcoords[0] * 180.0/np.pi-90)
    m.scatter(x,y,color='k',s=10.0,zorder=1000)
    for pix in np.arange(npix):
        step = 20
        edges = hp.boundaries(NSIDE,pix,step=step,nest=True)
        edgetheta,edgephi = np.arcsin(edges[2]), np.arctan2(edges[1],edges[0])
        edgephi[edgephi < 0] = 2 * np.pi +  edgephi[edgephi < 0]
        edge_length = np.arccos(edges[0][::step][1:] * edges[0][::step][:-1] + \
                                edges[1][::step][1:] * edges[1][::step][:-1] + \
                                edges[2][::step][1:] * edges[2][::step][:-1])

        print edgephi * 180.0/np.pi,"hehe"
        try:
            lon = m.shiftdata(edgephi * 180.0/np.pi)
        except ValueError:
            lon = edgephi * 180.0/np.pi
        lat = edgetheta *180.0/np.pi
        xx,yy = m(lon,lat)
        m.plot(xx,yy,marker='.',markersize=1.0,color='b')
    
    
    plt.show()

    m = Basemap(projection='ortho',lat_0=lat_0,lon_0=lon_0,resolution='l')
    #m = Basemap(projection='hammer',lon_0=180)
    m.drawmeridians(np.arange(0,360,20),latmax=90,alpha=0.1,color='gray',linewidth=1.0,dashes=[10000,0.01])
    m.drawparallels(np.arange(-80,80,10),alpha=0.1,color='gray',linewidth=1.0,dashes=[10000,0.01])
    xx,yy = m(pixelcoords[1]* 180.0/np.pi, pixelcoords[0] * 180.0/np.pi-90)
    m.scatter(xx,yy,color='k',s=10.0,zorder=1000)

    sv = SphericalVoronoi(np.asarray([xpix,ypix,zpix]).T)
    sv.sort_vertices_of_regions()
    verttheta, vertphi = np.arcsin(sv.vertices[:,2]), np.arctan2(sv.vertices[:,1],sv.vertices[:,0])
    vertphi[vertphi < 0] = 2 * np.pi +  vertphi[vertphi < 0]
    
    for region in sv.regions:
        polygon = Poly3DCollection([sv.vertices[region]], alpha=1.0)
        
        verttheta, vertphi = np.arcsin(sv.vertices[region,2]), np.arctan2(sv.vertices[region,1],sv.vertices[region,0])
        vertphi[vertphi < 0] = 2 * np.pi +  vertphi[vertphi < 0]
        verttheta = np.append(verttheta,verttheta[0])
        vertphi = np.append(vertphi,vertphi[0])
        edge_length = np.arccos(sv.vertices[region,0][1:] * sv.vertices[region,0][:-1] + \
                                sv.vertices[region,1][1:] * sv.vertices[region,1][:-1] + \
                                sv.vertices[region,2][1:] * sv.vertices[region,2][:-1])
        print edge_length.shape
        plot_edge = np.ones(verttheta.shape[0]).astype(bool)
        plot_edge[1:-1] = edge_length < 3*np.sqrt(4.0/npix)
        xx,yy = m(vertphi[plot_edge] * 180.0/np.pi, verttheta[plot_edge] * 180.0/np.pi)
        m.plot(xx,yy,color='g')

    verttheta, vertphi = np.arcsin(sv.vertices[:,2]), np.arctan2(sv.vertices[:,1],sv.vertices[:,0])
    vertphi[vertphi < 0] = 2 * np.pi +  vertphi[vertphi < 0]
    xx,yy = m(vertphi * 180.0/np.pi, verttheta * 180.0/np.pi)
    m.plot(xx,yy,marker='.',linestyle='None',markersize=5.0,color='g')
    plt.show()


