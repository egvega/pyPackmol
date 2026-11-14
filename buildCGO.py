"""
PackmolVisualizer
Described at PyMOL wiki: http://www.pymolwiki.org/index.php/
Authors: Esteban Gabriel Vega-Hissi
email: egvega@gmail.com, egvega@unsl.edu.ar
Date: February 2019
License: GNU General Public License
Version 0.1
"""


import pymol
from pymol import cmd

import pymol.cgo as cgo
from pymol.cgo import *

import numpy as np


def crisscross(x,y,z,d=0.5,name="crisscross",color=(1.,1.,1.)):

    obj = [
    LINEWIDTH, 3,
    COLOR, *color,

    BEGIN, LINE_STRIP,
    VERTEX, float(x-d), float(y), float(z),
    VERTEX, float(x+d), float(y), float(z),
    END,

    BEGIN, LINE_STRIP,
    VERTEX, float(x), float(y-d), float(z),
    VERTEX, float(x), float(y+d), float(z),
    END,

    BEGIN, LINE_STRIP,
    VERTEX, float(x), float(y), float(z-d),
    VERTEX, float(x), float(y), float(z+d),
    END

    ]
    view = cmd.get_view()
    cmd.load_cgo(obj,name)
    cmd.set_view(view)

def cgo_cube(x, y, z, r, color=(1.,1.,1.)):
    #r *= 3 ** -.5
    xmax=x+r; ymax=y+r; zmax=z+r
    return cgo_prisma(x, y, z, xmax,ymax,zmax, color)
def cgo_prisma(x, y, z, xmax, ymax, zmax, color=(1.,1.,1.)):
    #r *= 3 ** -.5
    return [
        cgo.COLOR, *color,
        cgo.BEGIN, cgo.TRIANGLE_STRIP,
        cgo.NORMAL, 0., 0., 1.,
        cgo.VERTEX, xmax , ymax , zmax ,
        cgo.VERTEX, xmax , y    , zmax ,
        cgo.VERTEX, x    , ymax , zmax ,
        cgo.VERTEX, x    , y    , zmax ,
        cgo.END,
        cgo.BEGIN, cgo.TRIANGLE_STRIP,
        cgo.NORMAL, 1., 0., 0.,
        cgo.VERTEX, xmax , y    , z    ,
        cgo.VERTEX, xmax , ymax , z    ,
        cgo.VERTEX, xmax , y    , zmax ,
        cgo.VERTEX, xmax , ymax , zmax ,
        cgo.END,
        cgo.BEGIN, cgo.TRIANGLE_STRIP,
        cgo.NORMAL, 0., 1., 0.,
        cgo.VERTEX, xmax , ymax , z    ,
        cgo.VERTEX, x    , ymax , z    ,
        cgo.VERTEX, xmax , ymax , zmax ,
        cgo.VERTEX, x    , ymax , zmax ,
        cgo.END,
        cgo.BEGIN, cgo.TRIANGLE_STRIP,
        cgo.NORMAL, 0., 0., -1.,
        cgo.VERTEX, x    , y    , z    ,
        cgo.VERTEX, x    , ymax , z    ,
        cgo.VERTEX, xmax , y    , z    ,
        cgo.VERTEX, xmax , ymax , z    ,
        cgo.END,
        cgo.BEGIN, cgo.TRIANGLE_STRIP,
        cgo.NORMAL, -1., 0., 0.,
        cgo.VERTEX, x    , ymax , zmax ,
        cgo.VERTEX, x    , y    , zmax ,
        cgo.VERTEX, x    , ymax , z    ,
        cgo.VERTEX, x    , y    , z    ,
        cgo.END,
        cgo.BEGIN, cgo.TRIANGLE_STRIP,
        cgo.NORMAL, 0., -1., 0.,
        cgo.VERTEX, x    , y    , zmax ,
        cgo.VERTEX, xmax , y    , zmax ,
        cgo.VERTEX, x    , y    , z    ,
        cgo.VERTEX, xmax , y    , z    ,
        cgo.END,
    ]

def cgo_cube_lines(xmin, ymin, zmin, d, color=(1.,1.,1.)):
    xmax=xmin+d; ymax=ymin+d; zmax=zmin+d
    return cgo_prisma_lines(xmin, ymin, zmin, xmax,ymax,zmax, color)

def cgo_prisma_lines(xmin, ymin, zmin, xmax,ymax,zmax, color=(1.,1.,1.)):
    return [
        COLOR, *color,
        BEGIN, LINE_STRIP, 
        VERTEX, xmin, ymin, zmin,
        VERTEX, xmax, ymin, zmin,
        VERTEX, xmax, ymax, zmin,
        VERTEX, xmin, ymax, zmin,
        VERTEX, xmin, ymin, zmin,
        END,
        BEGIN, LINE_STRIP,
        VERTEX, xmin, ymin, zmax,
        VERTEX, xmax, ymin, zmax,
        VERTEX, xmax, ymax, zmax,
        VERTEX, xmin, ymax, zmax,
        VERTEX, xmin, ymin, zmax,
        END,
        BEGIN, LINES,        
        VERTEX, xmin, ymin, zmin,
        VERTEX, xmin, ymin, zmax,
        VERTEX, xmax, ymin, zmin,
        VERTEX, xmax, ymin, zmax,
        VERTEX, xmin, ymax, zmin,
        VERTEX, xmin, ymax, zmax,        
        VERTEX, xmax, ymax, zmin,
        VERTEX, xmax, ymax, zmax,        
        END,
    ]

#Packmol to Pymol CGO  t r a n s l a t e

def sphere_translate_v0(s_pack, color=(1.,1.,1.)):
    s_py=[COLOR, *color, SPHERE]+s_pack
    return s_py
#three orthogonal circles composed by VERTEX connected by lines
def sphere_translate_v1(s_pack, color=(1.,1.,1.)):
    x,y,z,r=s_pack
    s_py = [ COLOR, *color, BEGIN, LINES ]
    for i in range(180):
        s_py.append( VERTEX )
        s_py.append(r*np.cos(i) + x )
        s_py.append(r*np.sin(i) + y )
        s_py.append(z)
        s_py.append( VERTEX )
        s_py.append(r*np.cos(i+0.1) + x )
        s_py.append(r*np.sin(i+0.1) + y )
        s_py.append(z)

    for i in range(180):
        s_py.append( VERTEX )
        s_py.append(r*np.cos(i) + x )
        s_py.append(y)
        s_py.append(r*np.sin(i) + z )
        s_py.append( VERTEX )
        s_py.append(r*np.cos(i+0.1) + x )
        s_py.append(y)
        s_py.append(r*np.sin(i+0.1) + z )

    for i in range(180):
        s_py.append( VERTEX )
        s_py.append(x)
        s_py.append(r*np.cos(i) + y )
        s_py.append(r*np.sin(i) + z )
        s_py.append( VERTEX )
        s_py.append(x)
        s_py.append(r*np.cos(i+0.1) + y )
        s_py.append(r*np.sin(i+0.1) + z )
    s_py.append(END)
    return s_py
#based in last ellipsoid
def sphere_translate(s_pack, color=(1.,1.,1.)):
    s_py=ellipsoid_translate(s_pack[:3]+[1.,1.,1.]+[s_pack[-1]], color)
    return s_py

def ellipsoid_translate_v0(e_pack, color=(1.,1.,1.)):
    e_py=[COLOR, *color, ELLIPSOID]+e_pack[:3]+[e_pack[-1]]+[e_pack[3],0.,0.]+[0.,e_pack[4],0.]+[0.,0.,e_pack[5]]
    return e_py
def ellipsoid_translate(e_pack, color=(1.,1.,1.)):
    cx, cy, cz = e_pack[:3] # center
    a1,a2,a3 = e_pack[3]*e_pack[-1], e_pack[4]*e_pack[-1], e_pack[5]*e_pack[-1] #semiaxes length (in Angstrom)
    u1, u2, v1, v2 = np.radians(-90), np.radians(90), np.radians(-180), np.radians(180)
    u_segs, v_segs = 8, 8 #segments
    # increments
    dU = (u2 - u1) / u_segs
    dV = (v2 - v1) / v_segs

    e_py = [COLOR, *color, BEGIN, LINES]

    U = u1
    for Y in range(0, u_segs):
            V = v1
            for X in range(0, v_segs):
                x, y, z = [0] * 5, [0] * 5, [0] * 5 #
                for i in range(1, 5): # creates four sets of ellipsoid points
                    u=U + dU * (i == 2 or i == 3)
                    v=V + dV * (i == 3 or i == 4)
                    #coordinates of ellipsoid points, packmol ellipsoid definition do not forsee rotation
                    x[i] = a1 * np.cos(u) * np.cos(v) + cx
                    y[i] = a2 * np.cos(u) * np.sin(v) + cy
                    z[i] = a3 * np.sin(u) + cz
                for i in (1, 2): # extends the pymol cgo object list with the new points # 1, 2, 4, 2, 3, 4
                    e_py.extend([VERTEX, x[i], y[i], z[i]]) # color set at the beginning of the list
                V += dV #increment
            U += dU #increment
    e_py.append(END)
    return e_py

def prisma_translate(p_pack, color=(1.,1.,1.)):
    if len(p_pack)>4: #prisma
        x0,y0,z0,x1,y1,z1=p_pack
        return cgo_prisma_lines(x0,y0,z0,x1,y1,z1, color)
    else: #cube
        x0,y0,z0,d=p_pack
        return cgo_prisma_lines(x0,y0,z0,x0+d,y0+d,z0+d, color)

def cylinder_translate_old(c_pack, color=(1.,1.,1.)):
    x1,y1,z1,x2,y2,z2=c_pack[:6]
    dx=x2-x1; dy=y2-y1; dz=z2-z1
    li=np.sqrt(dx**2+dy**2+dz**2)
    t=c_pack[-1]/li
    x2=x1+t*dx
    y2=y1+t*dy
    z2=z1+t*dz
    c_py=[COLOR, *color, CYLINDER]+c_pack[:3]+[x2,y2,z2]+[c_pack[-2]] + [1.,0.,0.,0.,0.,1.]#[r1,g1,b1,r2,g2,b2]
    return c_py

#using circle and lines
def cylinder_translate(c_pack, color=(1.,1.,1.)):
    # packmol: base point (a1,b1,c1), axis direction (a2,b2,c2), radius d, length l
    point1 = np.array(c_pack[:3], dtype=np.float64)
    axis = np.array(c_pack[3:6], dtype=np.float64)
    radius, length = c_pack[-2], c_pack[-1]
    norm = np.linalg.norm(axis)
    if norm == 0:
        print('Wrong cylinder: axis direction is zero')
        return [COLOR, *color]
    k = axis / norm
    vector = length * k  # base to top
    # orthonormal basis of the circle plane, valid for any axis orientation
    e = np.zeros(3)
    e[np.argmin(np.abs(k))] = 1.
    v = np.cross(k, e)
    v /= np.linalg.norm(v)
    y = np.cross(k, v)

    theta = np.linspace(0, 2 * np.pi, 100) # 100 should be changed by a user parameter
    x_circle = point1[0] + radius * (np.cos(theta) * y[0] + np.sin(theta) * v[0])
    y_circle = point1[1] + radius * (np.cos(theta) * y[1] + np.sin(theta) * v[1])
    z_circle = point1[2] + radius * (np.cos(theta) * y[2] + np.sin(theta) * v[2])

    x_c2=x_circle+vector[0]
    y_c2=y_circle+vector[1]
    z_c2=z_circle+vector[2]

    #do the magic VERTEX
    c_py = [ COLOR, *color, BEGIN, LINES ]
    for i in range(len(x_circle)): #circle at base
        c_py.append( VERTEX )
        c_py.append(x_circle[i])
        c_py.append(y_circle[i])
        c_py.append(z_circle[i])
    for i in range(len(x_circle)): #circle at top
        c_py.append( VERTEX )
        c_py.append(x_c2[i])
        c_py.append(y_c2[i])
        c_py.append(z_c2[i])
    for i in range(0,len(x_circle),20): #borders, 5 lines, related to circle points
        c_py.append( VERTEX )
        c_py.append(x_circle[i])
        c_py.append(y_circle[i])
        c_py.append(z_circle[i])
        c_py.append( VERTEX )
        c_py.append(x_c2[i])
        c_py.append(y_c2[i])
        c_py.append(z_c2[i])

    c_py.extend([END])
    return c_py

def plane_translate(p_pack,iMesh=0.,jMesh=50., color=(1.,1.,1.)):
    #plane eq.  ax + by + cz - d = 0
    #i,j Mesh  used to define plane length ---> ? take from .inp other restraint like inside box
    a,b,c,d=p_pack
    #seed points
    x = np.linspace(iMesh, jMesh, 50)
    y = np.linspace(iMesh, jMesh, 50)
    X, Y = np.meshgrid(x, y)
    if c==0:
        if a!=0:
            Z=X
            X = (d - b * Y - c * Z) / a
        else:
            Z=Y
            Y = (d - a * X - c * Z) / b
    else:
        Z = (d - a * X - b * Y) / c
    ## plot the mesh. Each array is 2D, so we flatten them to 1D arrays
    #fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')
    #ax.plot(X.flatten(), Y.flatten(), Z.flatten(), 'bo ',ms=2.)
    p_py=[COLOR, *color]
    for x,y,z in zip(X.flatten(),Y.flatten(),Z.flatten()):
        p_py.extend([BEGIN,POINTS,VERTEX,x,y,z,END])
    return p_py

def plane_translate2(p_pack,iMesh=0.,jMesh=50.,iMesh2=0.,jMesh2=50., color=(1.,1.,1.)):
    #plane eq.  ax + by + cz - d = 0
    #i,j Mesh  used to define plane length ---> ? take from .inp other restraint like inside box
    a,b,c,d=p_pack
    #seed points
    x = np.linspace(iMesh, jMesh, 50)
    y = np.linspace(iMesh2, jMesh2, 50)
    X, Y = np.meshgrid(x, y)
    if c==0:
        if a!=0:
            Z=X
            X = (d - b * Y - c * Z) / a
        else:
            Z=Y
            Y = (d - a * X - c * Z) / b
    else:
        Z = (d - a * X - b * Y) / c
    ## plot the mesh. Each array is 2D, so we flatten them to 1D arrays
    #fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')
    #ax.plot(X.flatten(), Y.flatten(), Z.flatten(), 'bo ',ms=2.)
    p_py=[COLOR, *color]
    for x,y,z in zip(X.flatten(),Y.flatten(),Z.flatten()):
        p_py.extend([BEGIN,POINTS,VERTEX,x,y,z,END])
    return p_py

def xygauss_translate(xygs_pack, color=(1.,1.,1.)):
    a1,b1,a2,b2,c0,h=xygs_pack
    A = h
    x0 = a1
    y0 = b1

    sigma_X = a2
    sigma_Y = b2

    xg = np.linspace(a1-2*a2,a1+2*a2,num=100)
    yg = np.linspace(b1-2*b2,b1+2*b2,num=100)

    theta= np.pi

    X, Y = np.meshgrid(xg,yg)

    a = np.cos(theta)**2/(2*sigma_X**2) + np.sin(theta)**2/(2*sigma_Y**2);
    b = -np.sin(2*theta)/(4*sigma_X**2) + np.sin(2*theta)/(4*sigma_Y**2);
    c = np.sin(theta)**2/(2*sigma_X**2) + np.cos(theta)**2/(2*sigma_Y**2);

    aXXdet = np.array([a*(Xi-x0)**2 for Xi in X],float)
    bbXYdet = np.array([2*b*(Xi-x0)*(Y[ii]-y0) for ii,Xi in enumerate(X)],float)
    cYYdet = np.array([c*(Yi-y0)**2 for Yi in Y],float)
    Z = c0+np.array([A*np.exp( - (ai + bbXYdet[i] + cYYdet[i])) for i,ai in enumerate(aXXdet)],float);

    xygs_py=[COLOR, *color]
    for x,y,z in zip(X.flatten(),Y.flatten(),Z.flatten()):
        xygs_py.extend([BEGIN,POINTS,VERTEX,x,y,z,END])
    return xygs_py

#CGO dictionary of translate functions - associated keyword to function
CGO_DICT={'sphere':sphere_translate,'ellipsoid':ellipsoid_translate,'cylinder':cylinder_translate,
          'box':prisma_translate,'cube':prisma_translate,'plane':plane_translate2, 'xygauss':xygauss_translate}

def test():
    #testing.....
    sphere_packmol=[0., 0., 0., 1.5] #sample
    ellipsoid_packmol=[0.,0.,0.,1.,0.5,0.5,6.5] #sample
    cube_packmol=[0.,0.,0.,1.]
    prisma_packmol=[0.,0.,0.,1.,2.,1.]
    cylinder_packmol=[0.,0.,0.,2.,2.,2.,1.,5.]
    plane_packmol=[0., 0., 1., 2.]
    xygs_pack=[20.,20.,10.,10.,10.,30.]

    #SPHERE
    sphere_pymol=sphere_translate(sphere_packmol)
    cmd.load_cgo(sphere_pymol,'sphe')
    #ELLIPSOID
    ellipsoid_pymol=ellipsoid_translate(ellipsoid_packmol)
    cmd.load_cgo(ellipsoid_pymol,'ellip')
    #CUBE / PRISMA
    prisma_pymol=prisma_translate(cube_packmol)
    cmd.load_cgo(prisma_pymol,'prisma')
        #Inserted at the beginning, all faces same color and transparency
        #list.insert index, itm (first in last out)   # obj.insert(0,1.)
        #obj=[ALPHA,0.75,COLOR,0.,0.,0.1]+obj
    #CYLINDER
    cylinder_pymol=cylinder_translate(cylinder_packmol)
    cmd.load_cgo(cylinder_pymol,'cyl')
    #PLANE
    plane_pymol=plane_translate(plane_packmol)
    cmd.load_cgo(plane_pymol,'plane')
    plane_pymol=plane_translate([0.,0.,1.,12.])
    cmd.load_cgo(plane_pymol,'plane2')
    #XYGaussianSurface
    xygs_pymol=xygauss_translate(xygs_pack)
    cmd.load_cgo(xygs_pymol,'xygs')
