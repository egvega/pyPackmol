import numpy as np
from scipy.spatial import ConvexHull
from scipy.spatial.distance import cdist
from pymol import cmd, stored

#check if prism is completelly inside sphere
def totalInsideSphere(sphere,prism):
    prism=np.float64(prism)
    sphere=np.float64(sphere)
    r2=sphere[-1]**2
    l=prism[-3:]-prism[:3]
    if np.max(l)<sphere[-1]*2: #max axis lenght shorter than Diameter
        #print (prism[:3]-sphere[:3])**2
        #print r2, np.sum((prism[:3]-sphere[:3])**2)<=r2
        if np.sum((prism[:3]-sphere[:3])**2)<=r2 and np.sum((prism[-3:]-sphere[:3])**2)<=r2: #corners inside sphere
            return True
    # cube from sphere
    #cubeSph=sphere[:3]-sphere[-1]
    #cubeSph=np.append(cubeSph,sphere[:3]+sphere[-1])
    #t=cubeSph[:3]>=prism[:3]
    #t=np.append(t,cubeSph[-3:]<=prism[-3:])
    return False

#check if sphere is completelly inside prism
def totalInsidePrism(sphere,prism):
    prism=np.float64(prism)
    sphere=np.float64(sphere)
    # cube from sphere
    cubeSph=sphere[:3]-sphere[-1]
    cubeSph=np.append(cubeSph,sphere[:3]+sphere[-1])
    t=cubeSph[:3]>=prism[:3]
    t=np.append(t,cubeSph[-3:]<=prism[-3:])
    return np.sum(t)==6

def volCGO(element):    
    if element[0]=='sphere':
        vol=4/3.*np.pi*(np.float64(element[-1])**3)
        dim1=np.float64(element[-1])*2 #diameter
    elif element[0]=='ellipsoid':
        print(element[4:-1],'-',element[-1])
        t=np.prod(np.float64(element[4:-1]))*np.float64(element[-1])**3 # a*d * b*d * c*d
        vol=4/3.*np.pi*t
        dim1=np.max(np.float64(element[4:-1]))*np.float64(element[-1]) #major axis times volume multiplier
    elif element[0]=='cube':
        vol=np.float64(element[-1])**3
        dim1=np.float64(element[-1])*np.sqrt(3) #body diagonal
    elif element[0]=='box':        
        vol=np.prod(np.abs(np.float64(element[-3:])-np.float64(element[1:4])))
        dim1=np.sqrt(np.sum(np.square(np.float64(element[-3:])-np.float64(element[1:4])))) #body diagonal
    elif element[0]=='cylinder':
        vol=np.pi*np.float64(element[-2])**2*np.float64(element[-1]) # area of base (circle) * length
        dim1=np.float64(element[-1]) #length
        dim2=np.float64(element[-2]) #radius
    elif element[0]=='plane' or element[0]=='xygauss' or element[0]=='fixed':
        vol=0 #print('Plane has no volume! \n')
        dim1=0
    return (vol,dim1) # returns a tuple of (vol,max_dim)

def interVolPrisms(p1,p2):
    #min corner has to be first
    p1=np.float64(p1)
    p2=np.float64(p2)
    inter=1
    for i in range(3):
            p1[i]=max(p1[i],p2[i])
            p1[i+3]=min(p1[i+3],p2[i+3])
            inter*=max(0.,p1[i+3]-p1[i])
    return inter
#spheres intersection
def interVolSpheres(sp1,sp2):
    r=np.float64(sp1[3])
    R=np.float64(sp2[3])
    d=np.sqrt(np.sum((np.float64(sp2[:3])-np.float64(sp1[:3]))**2))
    if d<np.abs(R-r):
        v=4/3.*np.pi*(min(r,R)**3)
    elif d>=r+R:
        v=0
    else:
        v=np.pi*(R+r-d)**2
        v*=(d**2+2*d*r-3*r**2+2*d*R+6*r*R-3*R**2)/(12.*d)
    return v

def interVolSP(gridscale=1.,sphere=np.array([0.,0.,0.,1.]), prism=np.array([0.,0.,0.,1.,1.,1.])):
    prism=np.float64(prism)
    sphere=np.float64(sphere)

    # cube containing sphere to get min, max for bounding box
    cubeSph=sphere[:3]-sphere[-1]
    cubeSph=np.append(cubeSph,sphere[:3]+sphere[-1])
    cubeSphVol=np.prod(cubeSph[-3:]-cubeSph[:3])
    #print 'cubeSph volume ', cubeSphVol
    #print sphereVol/cubeSphVol
    
    #bounding box
    bbox=np.zeros(6)
    for i in range(3):
        bbox[i]=min(prism[i],cubeSph[i])
        bbox[i+3]=max(prism[i+3],cubeSph[i+3])
    bbox=np.reshape(bbox,(2,3))

    xn, yn, zn, = np.int_(np.ceil((bbox[1]-bbox[0])/gridscale))

    origin = bbox[0]
    
    values=np.zeros(xn*yn*zn)

    gridpoint = {} # for each grid point have a list of atoms it is in

    for i in range(xn): # looping over the cube about the center of the atom 
        x = (i * gridscale) + origin[0] # convert point to xyz position
        for j in range(yn):
            y = (j * gridscale) + origin[1] 
            for k in range(zn):
                z = (k * gridscale) + origin[2] 
                dist = distance([x,y,z],sphere[:3]) #distance between 2 vectors
                if (dist <= sphere[-1]): #point within atom radius
                    count = i*(yn*zn)+j*(zn) + k
                    if not (count in gridpoint): # if point not assigned, init dictionary key
                        gridpoint[count] = ['sp']
                    else: # point already in dic, shared by atoms
                        gridpoint[count].append('sp')
                p=np.array([x,y,z])
                if np.sum(p>=prism[:3])==3 and np.sum(p<=prism[-3:])==3:#point inside box
                    count = i*(yn*zn)+j*(zn) + k
                    if not (count in gridpoint): # if point not assigned, init dictionary key
                        gridpoint[count] = ['cub']
                    else: # point already in dic, shared by atoms
                        gridpoint[count].append('cub')

    interNpts = 0
    for key in gridpoint:
        if len(gridpoint[key])>1:
            interNpts+=1
    boxNpts = len(values)
    boxVol = xn*gridscale*yn*gridscale*zn*gridscale # volume of bbox

    interVol = float(interNpts)/float(boxNpts)*boxVol # final volume calculator cubic Angstroms

    return interVol


class Sphere:
    def __init__(self,x,y,z,radius):        
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.radius = float(radius)

def distance(v1,v2):
    if (len(v1)!=len(v2)):
       print("Vectors not same lenght!")
    v1=np.array(v1)
    v2=np.array(v2)
    return np.sqrt(sum((v2-v1)**2.0))


def calcVolume(mol='(all)',gridscale=0.25,tolerance=0.0): #tolerance for radius 
    # mol: selection
    # gridscale: the grid spacing (for all axis)
    gridscale=float(gridscale)
    tolerance=float(tolerance)
    # ### Generate sphere(atom) list from mol selection for current state
    sphs = []
    cmd.iterate_state(-1, mol, 'sphs.append(Sphere(x,y,z,vdw+tolerance))',space={'sphs':sphs,'Sphere':Sphere,'tolerance':tolerance}) #-1 means current state #1.06

    maxRadius=0.0
    for s in sphs:
        if s.radius>maxRadius:
            maxRadius=s.radius
    sphslen=len(sphs)
    if sphslen == 0:
        print("Warning: Sphere List is Emtpy")
    elif sphslen > 100: # actual atom number may be adjusted until faster method or includ faster method for proteins
        print("Warning: molecule too large. Calculation will take too long. Aborting...")
        return 0.0

    # ### Generate bounding box
    buff = gridscale # add/substract a buffer to make sure spheres are completly contained
    boundaries=np.array(cmd.get_extent(mol))
    boundaries[0]-=(maxRadius+buff) #min_x,min_y,min_z
    boundaries[1]+=(maxRadius+buff) #max_x,max_y,max_z
    #print "max corner = ", boundaries[1]
    #print "min corner = ", boundaries[0]
    
    # Define the bounding box.
    # max-min gets length, divided by scale gets points (needs round up integer value)
    xn, yn, zn, = np.int_(np.ceil((boundaries[1]-boundaries[0])/gridscale))    
    origin = boundaries[0] #[min_x, min_y, min_z]
    # Generate empty grid.
    values=np.zeros(xn*yn*zn)
    
    #print gridscale ,xn,yn,zn,origin
    
    # ### loop over each sphere. 
    dict_gridpoint = {} # for each grid point have a list of atoms it is in

    for atom_i,atom in enumerate(sphs):        
        radius = atom.radius

        grid_i = round((atom.x - origin[0] ) / gridscale) #get atom location in grid
        grid_j = round((atom.y - origin[1] ) / gridscale)
        grid_k = round((atom.z - origin[2] ) / gridscale)

        if (grid_i<0.0 or grid_j<0.0 or grid_k<0.0):
           print("ERROR zero grid lenght points")

        radius_gridpoint = np.ceil(radius / gridscale) + 1
        #print radius, gridscale, radius_gridpoint
        #print  grid_i,grid_j, grid_k

        ini_i = max([(grid_i - radius_gridpoint),0.0]) #get points boundaries of each atom in each axis
        end_i  = min([(grid_i + radius_gridpoint),xn])
        ini_j = max([(grid_j - radius_gridpoint),0.0])
        end_j  = min([(grid_j + radius_gridpoint),yn])
        ini_k = max([(grid_k - radius_gridpoint),0.0])
        end_k  = min([(grid_k + radius_gridpoint),zn])

        for i in range(int(ini_i),int(end_i)): # looping over the cube about the center of the atom 
            x = (i * gridscale) + origin[0] # convert point to xyz position
            for j in range(int(ini_j),int(end_j)):
                y = (j * gridscale) + origin[1] 
                for k in range(int(ini_k),int(end_k)):
                    z = (k * gridscale) + origin[2] 
                    dist = distance([x,y,z],[atom.x,atom.y,atom.z]) #distance between 2 vectors
                    if (dist <= radius): #point within atom radius
                        count = i*(yn*zn)+j*(zn) + k
                        if not (count in dict_gridpoint): # if point not assigned, init dictionary key
                           dict_gridpoint[count] = [atom_i]
                        else: # point already in dic, shared by atoms
                           dict_gridpoint[count].append(atom_i)

    molNpts = len(dict_gridpoint.keys()) #number of keys is the numper of grid points included into molecule volume
    boxNpts = len(values)
    boxVol = xn*gridscale*yn*gridscale*zn*gridscale # volume of prisma

    molVol = float(molNpts)/float(boxNpts)*boxVol # final volume calculator cubic Angstroms
    #print molNpts, boxNpts, boxVol
    #print molVol
    
    return molVol

def generate_sphere_points(n):
    """
    Returns an array of evenly distributed points on the 
    surface of a the unit sphere using the Golden-Section
    Spiral algorithm.
    
    Parameters
    ----------
    n : int. Number of points
    """
    
    golden_angle = np.pi * (3 - 5**0.5)

    theta = golden_angle * np.arange(n)
    z = np.linspace(1 - 1.0 / n, 1.0 / n - 1, n)

    radius = np.sqrt(1 - z * z)
    points = np.zeros((n, 3))
    # change from polar to cartesian coordinates
    points[:,0] = radius * np.cos(theta)
    points[:,1] = radius * np.sin(theta)
    points[:,2] = z
    
    return points

def get_Vpoints(coord, probe=1.4, n_points=1000):
    """
    Returns array of coordinates of all points of sasa method

    Parameters
    ----------
    coord : Array (n,3). Cartesian coordinates
    probe : float. Radisu of the solvent
    n : int. Number of points
    """

    points = generate_sphere_points(n_points)
    # compute the area each point represents
    Vcoord = []
    for i in range(len(coord)):
        # scale the points to the correct radius
        radius = stored.vdW[i] + probe
        points_scaled = coord[i] + points * radius
        Vcoord.extend(points_scaled.tolist())
    return Vcoord

def get_neighbors(coord, probe, k):
    """
    Returns list of indices of neighbors. 
    
    Parameters
    ----------
    coord : Array (n,3). Cartesian coordinates
    probe : float. Radius of the solvent
    k     : index of particle
    """

    neighbor_indices = []
    coord_k = coord[k]
    radius = stored.vdW[k] + probe * 2 # take into account probe radius of the possible neighbor 
    # avoid comparing an atom with itself
    indices = list(range(k))
    indices.extend(list(range(k+1, len(coord))))
    # compute distance all vs all
    dist = cdist(coord, coord, metric='euclidean')    
    #
    for idx in indices:
        if dist[idx, k] < radius + stored.vdW[idx]:
            neighbor_indices.append(idx)
    return neighbor_indices

def get_surface_points(coord, probe=1.4, n_points=1000): # !! for volume calculation probe=0.0
    """
    Returns surface points of the molecular solvent-accessible-surface area
    
    Parameters
    ----------
    coord : Array (n,3). Cartesian coordinates
    probe : float. Radisu of the solvent
    n : int. Number of points
    """

    points = generate_sphere_points(n_points)
    # compute the area each point represents
    surfPoints = np.empty((0,3))
    
    for i in range(len(coord)):
        # scale the points to the correct radius
        radius = stored.vdW[i] + probe
        points_scaled = coord[i] + points * radius
        # get all the neighbors of the i residue
        neighbors = get_neighbors(coord, 1.4, i) #idx array of neighbors
        # compute the distance between points and neighbors
        d_matrix = cdist(points_scaled, coord[neighbors],
                     metric='euclidean')
        # create a matrix and store the vdW radii for each neighbor
        nb_matrix = np.zeros((n_points, len(neighbors)))
        for nb_i, nb in enumerate(neighbors):
            nb_matrix[:,nb_i] =  stored.vdW[nb]
        # compute the number of buried points, we have to be carefull
        # because we have counted how many times a point is buried
        # and we only need how many points are buried
        
        buried = d_matrix < nb_matrix + probe # bool array of buried points according criteria
        # buried points may be "inside" more than one atom, then
        # sum over atoms
        # if sum equals 0 we got a surface point 
        surfP=points_scaled[np.sum(buried,axis=1)==0] # pass as index array
        surfPoints=np.concatenate((surfPoints,surfP),axis=0)
    return surfPoints

class Body:
    def __init__(self,x,y,z):
        global chull
        from . import chull # legacy hull, only used by calcVolume2
        self.body=[]
        for i in range(len(x)):
            self.body.append(chull.Vector(float(x[i]),float(y[i]),float(z[i])))
    def getVolume(self):
        ref=chull.Vector(0,0,0)
        #create the ConvexHull
        h=chull.Hull(self.body)
        vT=0
        #loop over faces and sum up tetrahedron volume
        for fa in h.faces:
            vT+=h.Volume(fa,ref)/6
        return vT

def calcVolume2(mol='(all)',probe=1.4,n_points=100): #
    '''
    Returns molecular volume using solvent accesible surface and qhull. 
    
    Parameters
    ----------
    mol     : string. Pymol object name of Molecule
    probe   : float. Radius of the solvent
    n_points: integer. number of points of each sphere (more > precision > computation time) 
    '''
    # get molecule pymol object model and x,y,z coordinates
    model = cmd.get_model(mol)
    xyz = np.array(model.get_coord_list())

    # TODO cancel when molecular size too large; switch to: get outermost atoms, then process only those without Vdw.
    if len(xyz) < 200:
        # generate an array of van der Waals radius of each atom
        stored.vdW = []
        cmd.iterate(mol, 'stored.vdW.append(vdw)')
        # compute surface points
        surface_points=get_surface_points(xyz,probe,n_points)
    else: #bigger molecules
        hull=ConvexHull(xyz)        
        surface_points=xyz[hull.vertices] #pass index of points to get coords
    # calculate volume
    b=Body(surface_points[:,0],surface_points[:,1],surface_points[:,2])
    return b.getVolume()


def calcVolume3(mol='(all)',probe=1.4,n_points=100): #
    '''
    Returns molecular volume using solvent accesible surface and scipy ConvexHull.

    Parameters
    ----------
    mol     : string. Pymol object name of Molecule
    probe   : float. Radius of the solvent
    n_points: integer. number of points of each sphere (more > precision > computation time)
    '''
    # get molecule pymol object model and x,y,z coordinates
    model = cmd.get_model(mol)
    xyz = np.array(model.get_coord_list())

    stored.vdW = []
    cmd.iterate(mol, 'stored.vdW.append(vdw)')

    points=get_Vpoints(xyz, probe, n_points)

    hull=ConvexHull(points)
    return hull.volume
