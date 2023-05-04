# This source file is part of the chiraltube project, made by Jose M. de Albornoz Caratozzolo and Felipe Cervantes Sodi
# This source code is released under the 3-Clause BSD License, see "LICENSE.txt".

import sys
import math
import random as rd
import statistics as stat


opts = [opt for opt in sys.argv[1:] if opt.startswith("-")]
args = [arg for arg in sys.argv[1:] if not arg.startswith("-")]

ERROR=0.025

class point3d:
    '''
    3 dimensional point used for atoms (and vectors), with an
    identifier for the element.

    Examples
    --------
    >>> my_atom = point3d(0.0, 1.0, 0.5, 'H')
    '''
    def __init__(self, x, y, z, ele):
        self.x=x #X, Y, Z coordinates
        self.y=y
        self.z=z
        self.ele=ele #element
    def __str__(self):
        return "{}, {}, {}, {}".format(self.x, self.y, self.z, self.ele)
    def __repr__(self):
        return "({}, {}, {}, {})".format(self.x, self.y, self.z, self.ele)
    def __eq__(self, other): #equality, element is not taken into account
        return self.x==other.x and self.y==other.y and self.z==other.z
    def __add__(self, other): #sum
        return point3d(self.x+other.x, self.y+other.y, self.z+other.z, other.ele)
    def __mul__(self, other): #multiplication
        if isinstance(self, point3d):
            if isinstance(other, point3d):
                return point3d(self.x*other.x, self.y*other.y, self.z*other.z, self.ele)
            else:
                return point3d(self.x*other, self.y*other, self.z*other, self.ele)
        elif isinstance(other, point3d):
            return point3d(self*other.x, self*other.y, self*other.z, other.ele)
        else:
            return point3d(self*other.x, self*other.y, self*other.z, other.ele)
    def __sub__(self,other): #substraction
        return point3d(self.x-other.x, self.y-other.y, self.z-other.z, other.ele)
    def mag(self):#magnitude
        '''Returns the magnitude of the vector defined by the X, Y, Z coordinates of the atom'''
        return math.sqrt(self.x**2+self.y**2+self.z**2)


def distance(p1, p2):
    '''
    Returns the distance between two point3d's

    Parameters
    ----------
    p1: point3d
        First atom/vector
    p2: point3d
        Second aom/vector

    Returns
    -------
    float
        The Euclidean distance between the positions of p1 and p2

    See Also
    --------
    point3d: 3 dimensional point used for atoms (and vectors)

    Examples
    --------
    >>> p1 = point3d(0.0, 1.0, 0.5, 'Na')
    >>> p2 = point3d(-1.0, 0.0, 1.5, 'Cl')
    >>> distance(p1,p2)
    1.7320508075688772
    
    '''
    return math.sqrt((p2.x-p1.x)**2 + (p2.y-p1.y)**2 + (p2.z-p1.z)**2)

def testgraph(Arr):
    '''
    VERY basic visualization for testing purposes using matplotlib.

    Parameters
    ----------
    Arr: List of point3d's
        List of all the atoms that you want to graph

    '''
    import matplotlib.pyplot as plt
    from mpl_toolkits import mplot3d
    x=[ele.x for ele in Arr]
    y=[ele.y for ele in Arr]
    z=[ele.z for ele in Arr]
    C=[]
    for ele in Arr:
        if ele.ele=='S':
            C.append('b')
        elif ele.ele=='C':
            C.append('g')
        else:
            C.append('r')
    fig=plt.figure()
    ax=plt.axes(projection='3d')
    ax.scatter3D(x,y,z, c=C)
    plt.show()
    

def arr_initial(n, m, x, y):
    '''
    Returns initial array of atoms

    Repeats the Unit Cell given in the input file to create a 2d array of
    unit cells, with the right dimensions from which the nanoribbon/nanotube
    will be made.

    Parameters
    ----------
    n: int
        The first chiral index
    m: int
        The second chiral index
    x: int
        The first translational index
    y: int
        The second translational index

    Returns
    -------
    List of point3d
        Array of atoms arranged in a 2d plane.

    See Also
    --------
    robtainxy: Gives a list of possible x,y pairs from a specific n,m pair

    Example
    -------
    >>> n=5
    >>> m=7
    >>> theta, Res = robtainxy(n,m)
    >>> x, y, error, p_error=Res[0]
    >>> x=round(x)
    >>> y=round(y)
    >>> Arr=arr_initial(n, m, x, y)
    '''
    y_tot=n+x+1
    x_inf=-y-1
    x_sup=m+1
    Arr=[]
    py=0
    pz=0
    for i in range(y_tot+1):
        px=a2*x_inf+a1*i
        for _ in range(x_sup-x_inf+1):
            for atom in UnitCell:
                Arr.append(px+atom)            
            px=px+a2
    return Arr

def robtainxy(n, m, error=ERROR, MAX=100):
    '''
    Returns chiral angle theta, and an array with the X, Y pairs
    and their respective errors corresponding to the n,m indices,
    up to an error "error" and looks for x,y<"MAX".

    Parameters
    ----------
    n: int
        The first chiral index
    m: int
        The second chiral index
    error: float
        The allowed error when looking for integer x,y pairs. Its in
        units of the cell vectors.
    MAX: int
        The maximum range of search for x and y. Default is 100

    Returns
    -------
    theta: float
        The chiral angle defined by the indices n and m
    Res: List of floats
        A list with all the results of x, y pairs (not rounded to
        nearest integer), along with their respective errors and
        percentage errors.

    See Also
    --------
    obtainxy: Prints theta and the X, Y pairs instead of returning them

    Examples
    --------
    >>> n=5
    >>> m=7
    >>> theta, Res = robtainxy(n,m)
    >>> x, y, error, p_error=Res[0]
    >>> x=round(x)
    >>> y=round(y)
    >>> coordinatesNT_xyz(n,m,x,y, arch_out="./OutputFile.xyz")
    '''
    Res=[]
    if n==0:
        theta=0
        for x in range(1,MAX):
            y=x*a1.mag()*math.cos(phi)/a2.mag()
            if abs(round(y)-y)<error:
                Res.append((x, y, abs(round(y)-y), (abs(round(y)-y)/y)*100))
        for y in range(1,MAX):
            x=y*a2.mag()/(a1.mag()*math.cos(phi))
            if abs(round(x)-x)<error:
                Res.append((x, y, abs(round(x)-x), (abs(round(x)-x)/x)*100))
    elif m==0:
        theta=phi
        for x in range(1,MAX):
            y=x*a1.mag()/(a2.mag()*math.cos(phi))
            if abs(round(y)-y)<error:
                Res.append((x, y, abs(round(y)-y), (abs(round(y)-y)/y)*100))
        for y in range(1,MAX):
            x=y*a2.mag()*math.cos(phi)/a1.mag()
            if abs(round(x)-x)<error:
                Res.append((x, y, abs(round(x)-x), (abs(round(x)-x)/x)*100))
    else:
        theta= math.acos((a2.x*(a2*m+a1*n).x + a2.y*(a2*m+a1*n).y+a2.z*(a2*m+a1*n).z)/(a2.mag()*(a2*m+a1*n).mag()))
        for x in range(1, MAX):
            #y=((n*a1.mag()*a1.mag())/(m*a2.mag()*a2.mag()))*(math.sin(math.pi-2*phi+theta)/math.sin(theta))*x
            y=((n*a1.mag()*a1.mag())/(m*a2.mag()*a2.mag()))*(math.sin(phi-theta)*math.cos(phi-theta))*x/(math.sin(theta)*math.cos(theta))
            if abs(round(y)-y)<error:
                Res.append((x, y, abs(round(y)-y), (abs(round(y)-y)/y)*100))
        for y in range(1, MAX):
            #x=((m*a2.mag()*a2.mag())/(n*a1.mag()*a1.mag()))*(math.sin(theta)/math.sin(math.pi-2*phi+theta))*y
            x=((m*a2.mag()*a2.mag())/(n*a1.mag()*a1.mag()))*(math.sin(theta)*math.cos(theta))*y/(math.sin(phi-theta)*math.cos(phi-theta))
            if abs(round(x)-x)<error:
                Res.append((x, y, abs(round(x)-x), (abs(round(x)-x)/x)*100))
    Res.sort(key=sortfuncT)
    return theta, Res

def obtainxy(n, m, error=ERROR, MAX=100):
    '''
    Prints chiral angle theta, and an array with the X, Y pairs
    and their respective errors corresponding to the n,m indices,
    up to an error "error" and looks for x,y<"MAX".

    Parameters
    ----------
    n: int
        The first chiral index
    m: int
        The second chiral index
    error: float
        The allowed error when looking for integer x,y pairs. Its in
        units of the cell vectors.
    MAX: int
        The maximum range of search for x and y. Default is 100

    See Also
    --------
    robtainxy: Returns theta and the X, Y pairs instead of printing them
    '''
    theta, Res= robtainxy(n,m, error, MAX)
    print("Chiral angle: ", math.degrees(theta),  "\n\nResults up to an error <", error)
    print("X,    Y,    error, percentage error")
    for ele in Res:
        x,y,err,por=ele
        print("{}, {}, {}, {}%".format(round(x),round(y), err, por))
        
def obtainpq(n, m, x,y, MAX=100):
    '''
    Returns the symmetry vector R indices: P and Q as a tuple

    Parameters
    ----------
    n: int
        The first chiral index
    m: int
        The second chiral index
    x: int
        The first translational index
    y: int
        The second translational index

    Returns
    -------
    p: int
        First index of the symmetry vector
    q: int
        Second index of the symmetry vector

    See Also
    --------
    robtainxy: Gives the chiral angle and a list of
               possible x,y pairs from a specific n,m pair
    '''
    N=x*m + y*n
    for p in range(MAX):
        q=(1-y*p)/x
        if int(q)==q and (0< m*p - n*q <= N):
            break
    return p, int(q)

def rotate(Arr, theta):
    '''
    Rotates every point in the array (phi-theta) radians counter-clockwise.
    Returns the rotated array.

    Parameters
    ----------
    Arr: list of point3d
        List of atoms in a 2D plane that will be rotated.
    theta: float
        Chiral angle in radians

    Returns
    -------
    list of point3d
        Rotated array of atoms

    See Also
    --------
    robtainxy: Gives the chiral angle and a list of
               possible x,y pairs from a specific n,m pair

    Examples
    --------
    >>> n=5
    >>> m=7
    >>> theta, Res = robtainxy(n,m)
    >>> x, y, error, p_error=Res[0]
    >>> x=round(x)
    >>> y=round(y)
    >>> Arr=arr_initial(n, m, x, y)
    >>> Arr= rotate(Arr, theta)
    '''
    for ele in Arr:
        tmpx=ele.x
        tmpy=ele.y
        ele.x=tmpx*math.cos(phi-theta)-tmpy*math.sin(phi-theta)
        ele.y=tmpx*math.sin(phi-theta)+tmpy*math.cos(phi-theta)
    return Arr

def eliminate(Arr, n, m, x, y):
    '''
    Removes those atoms in the array "Arr" which do not belong in the
    nanoribbon. Returns the new array, the Y distance (|T|) and the X
    distance (|C|).

    Parameters
    ----------
    Arr: list of point3d
        Array of atoms in a rotated, 2D plane
    n: int
        The first chiral index
    m: int
        The second chiral index
    x: int
        The first translational index
    y: int
        The second translational index

    Returns
    -------
    new_Arr: list of point3d
        Array of atoms that do belong in the nanoribbon
    disy: float
        The Y distance (height, |T|) of the nanoribbon
    disx: float
        The X distance (length, |C|) of the nanoribbon

    See Also
    --------
    rotate: Rotates a 2d array of atoms (phi - theta) radians counter-clockwise.

    Examples
    --------
    >>> n=5
    >>> m=7
    >>> theta, Res = robtainxy(n,m)
    >>> x, y, error, p_error=Res[0]
    >>> x=round(x)
    >>> y=round(y)
    >>> Arr=arr_initial(n, m, x, y)
    >>> Arr= rotate(Arr, theta)
    >>> Arr, disy, disx = eliminate(Arr, n, m, x, y)
    '''
    new_Arr=[]
    disy= distance(point3d(0,0,0,0), (a1*x-a2*y))
    disx= distance(point3d(0,0,0,0), (a1*n+a2*m))
    for ele in Arr:
        if ele.y>=-0.000001 and ele.y<(disy-0.000001) and ele.x>=-0.000001 and ele.x<(disx-0.000001):
            new_Arr.append(ele)
    return new_Arr, disy, disx

def coordinates_xyz(n, m, x, y, prt=1, arch_out=False, grf=False):
    '''
    Does all the steps necessary to get a nanoribbon with chiral indices
    (n,m) and translational indices (x,y). If "prt"=1, it prints each
    completed step. It automatically prints to file "arch_out",
    defaults to "./Ribbon_n-m.xyz". Output file is in special .xyz format.

    Parameters
    ----------
    n: int
        The first chiral index
    m: int
        The second chiral index
    x: int
        The first translational index
    y: int
        The second translational index

    Optional
    --------    
    prt: int
        If set to 1 then it prints each completed step in the process
    arch_out: 
        Name of the output file, if False, it defaults to "./Ribbon_n-m.xyz"
    grf: int
        If set to 1 then it also graphs the result using testgraph

    See Also
    --------
    robtainxy: Gives the chiral angle and a list of
               possible x,y pairs from a specific n,m pair
    testgraph: VERY basic visualization for testing purposes using matplotlib.
    
    Examples
    --------
    >>> n=5
    >>> m=3
    >>> theta, Res = robtainxy(n,m)
    >>> x, y, error, p_error=Res[0]
    >>> x=round(x)
    >>> y=round(y)
    >>> coordinates_xyz(n,m,x,y, arch_out="./OutputFile.xyz")
    '''
    if not arch_out:
        arch_out= "./Ribbon_{}-{}.xyz".format(n,m)
    theta, Res = robtainxy(n, m)
    if prt==1:
        print("\nChiral angle obtained correctly: Theta=", math.degrees(theta), "\n")
    Arr=arr_initial(n, m, x, y)
    if prt==1:
        print("Initial array correctly created\n")
    Arr= rotate(Arr, theta)
    if prt==1:
        print("Rotation succesfully completed\n")
    Arr, disy, disx = eliminate(Arr, n, m, x, y)
    for ele in Arr:
        displaced = point3d(0,-disy, 0,'')+ ele
        for ele2 in Arr:
            if distance(displaced, ele2)<0.2:
                Arr.remove(ele)
    for ele in Arr:
        displaced = point3d(-disx,0,0,'')+ ele
        for ele2 in Arr:
            if distance(displaced, ele2)<0.2:
                Arr.remove(ele)
    for ele in Arr:
        displaced = point3d(-disx,-disy,0,'')+ ele
        for ele2 in Arr:
            if distance(displaced, ele2)<0.2:
                Arr.remove(ele)
    if prt==1:
        print("Ribbon created correctly, {} atoms in total.".format(len(Arr)))
    with open(arch_out, mode = 'w') as outfile:
        print(len(Arr), file=outfile)
        print("Nanoribbon (n,m)=({},{}). Y={} A, X={} A".format(n, m, disy, disx), file=outfile)
        for ele in Arr:
            print("{}    {}    {}    {}".format(ele.ele, ele.x, ele.y, ele.z), file=outfile)
    if grf:
        testgraph(Arr)

def coordinates_VASP(n, m, x, y, prt=1, arch_out=False, grf=False):
    '''
    Does all the steps necessary to get a nanoribbon with chiral indices
    (n,m) and translational indices (x,y). If "prt"=1, it prints each
    completed step. It automatically prints to file "arch_out",
    defaults to "./POSCAR". Output file is POSCAR VASP format.

    Parameters
    ----------
    n: int
        The first chiral index
    m: int
        The second chiral index
    x: int
        The first translational index
    y: int
        The second translational index

    Optional
    --------    
    prt: int
        If set to 1 then it prints each completed step in the process
    arch_out: 
        Name of the output file, if False, it defaults to "./Ribbon_n-m.xyz"
    grf: int
        If set to 1 then it also graphs the result using testgraph

    See Also
    --------
    robtainxy: Gives the chiral angle and a list of
               possible x,y pairs from a specific n,m pair
    testgraph: VERY basic visualization for testing purposes using matplotlib.
    
    Examples
    --------
    >>> n=5
    >>> m=3
    >>> theta, Res = robtainxy(n,m)
    >>> x, y, error, p_error=Res[0]
    >>> x=round(x)
    >>> y=round(y)
    >>> coordinates_VASP(n,m,x,y, arch_out="./OutputFile.xyz")
    '''
    if not arch_out:
        arch_out= "POSCAR"
    theta, Res = robtainxy(n, m)
    if prt==1:
        print("\nChiral angle obtained correctly: Theta=", math.degrees(theta), "\n")
    Arr=arr_initial(n, m, x, y)
    if prt==1:
        print("Initial array correctly created\n")
    Arr= rotate(Arr, theta)
    if prt==1:
        print("Rotation succesfully completed\n")
    Arr, disy, disx = eliminate(Arr, n, m, x, y)
    for ele in Arr:
        displaced = point3d(0,-disy, 0,'')+ ele
        for ele2 in Arr:
            if distance(displaced, ele2)<0.2:
                Arr.remove(ele)
    for ele in Arr:
        displaced = point3d(-disx,0,0,'')+ ele
        for ele2 in Arr:
            if distance(displaced, ele2)<0.2:
                Arr.remove(ele)
    for ele in Arr:
        displaced = point3d(-disx,-disy,0,'')+ ele
        for ele2 in Arr:
            if distance(displaced, ele2)<0.2:
                Arr.remove(ele)
    if prt==1:
        print("Ribbon created correctly, {} atoms in total.".format(len(Arr)))
    with open(arch_out, mode = 'w') as outfile:
        print("Nanoribbon (n,m)=({},{})".format(n, m), file=outfile)
        print(1, file=outfile)
        print("{} 0.0 0.0".format(disx), file=outfile)
        print("0.0 {} 0.0".format(disy), file=outfile)
        print("0.0 0.0 {}".format(a3.z), file=outfile)
        unique_elements = []
        for atom in Arr:
            if atom.ele not in unique_elements:
                unique_elements.append(atom.ele)
        Arr_divided=[[atom for atom in Arr if atom.ele==element] for element in unique_elements]
        ele_names=""
        ele_numbers=""
        for i in range(len(unique_elements)):
            ele_names=ele_names+unique_elements[i]+' '
            ele_numbers = ele_numbers + str(len(Arr_divided[i]))+' '
        print(ele_names.strip(), file=outfile)
        print(ele_numbers.strip(), file=outfile)
        print('Cartesian', file=outfile)
        for ele in Arr_divided:
            for atom in ele:
                print("{} {} {}".format(atom.x, atom.y, atom.z), file=outfile)
    if grf:
        testgraph(Arr)
        
def nanotube(Arr, n, m, z_med=False):
    '''
    Rolls the ribbon into a nanotube with chiral indices (n,m).
    "z_med" determines the central plane, defaults to a_3/2.
    Returns the new array and the NT radius.

    Parameters
    ----------
    Arr: list of point3d
        Array of atoms in a nanoribbon
    n: int
        First chiral index
    m: int
        Second chiral index
    z_med: float
        Determines the central plane in the z-direction. Defaults to a_3/2.

    Returns
    -------
    Arr: list of point3d
        Array of atoms folded into a nanotube
    radio: float
        Radius of the newly formed nanotube (measured at the central plane).

    Examples
    --------
    >>> n=5
    >>> m=7
    >>> theta, Res = robtainxy(n,m)
    >>> x, y, error, p_error=Res[0]
    >>> x=round(x)
    >>> y=round(y)
    >>> Arr=arr_initial(n, m, x, y)
    >>> Arr= rotate(Arr, theta)
    >>> Arr, disy, disx = eliminate(Arr, n, m, x, y)
    >>> Arr, radio=nanotube(Arr, n, m)
    '''
    x_max= distance(point3d(0,0,0,0), (a1*n+a2*m))
    radio= x_max/(2*math.pi)
    if z_med==False:
        z_med= (a3.z)/2        
    for ele in Arr:
        gamma= ele.x/radio
        tmpy=ele.y
        radio_rel= ele.z-z_med
        ele.x= (radio+radio_rel)*math.cos(gamma)
        ele.y= (radio+radio_rel)*math.sin(gamma)
        ele.z=tmpy
    return Arr, radio

def coordinatesNT_xyz(n, m, x, y, rep=0,prt=1, arch_out=False, grf=False):
    '''
    Does all the steps necessary to get a nanotube with chiral indices (n,m)
    and translational indices (x,y), repeated along the NT axis "repeat" times.
    If "prt"=1, it prints the completed steps. It automatically prints to file
    "arch_out", defaults to "./Nanotube_n-m.xyz". Output file is in special
    .xyz format.

    Parameters
    ----------
    n: int
        The first chiral index
    m: int
        The second chiral index
    x: int
        The first translational index
    y: int
        The second translational index

    Optional
    --------
    rep: int
        Number of times the nanotube is repeated along its axis.
        Defaults to 0 (no repetition)
    prt: int
        If set to 1 then it prints each completed step in the process
    arch_out: 
        Name of the output file, if False, it defaults to "./Nanotbe_n-m.xyz"
    grf: int
        If set to 1 then it also graphs the result using testgraph

    See Also
    --------
    coordinates_xyz: Does the same but only makes the nanoribbon
    robtainxy: Gives the chiral angle and a list of
               possible x,y pairs from a specific n,m pair
    repeat: Repeats the array N times along the NT axis
    testgraph: VERY basic visualization for testing purposes using matplotlib.

    Examples
    --------
    >>> n=5
    >>> m=3
    >>> theta, Res = robtainxy(n,m)
    >>> x, y, error, p_error=Res[0]
    >>> x=round(x)
    >>> y=round(y)
    >>> coordinatesNT_xyz(n,m,x,y, arch_out="./OutputFile.xyz")
    '''
    extra=""
    theta, Res = robtainxy(n, m)
    if prt==1:
        print("\nChiral angle obtained correctly: Theta=", math.degrees(theta), "\n")
    Arr=arr_initial(n, m, x, y)
    if prt==1:
        print("Initial array correctly created\n")
    Arr= rotate(Arr, theta)
    if prt==1:
        print("Rotation succesfully completed\n")
    Arr, disy, disx = eliminate(Arr, n, m, x, y)
    if prt==1:
        print("Ribbon created correctly, {} atoms in total.\n".format(len(Arr)))
    Arr, radio=nanotube(Arr, n, m)
    for ele in Arr:
        displaced = point3d(0,0,-disy,'')+ ele
        for ele2 in Arr:
            if distance(displaced, ele2)<0.2:
                Arr.remove(ele)
    if prt==1:
        print("Nanotube created correctly, {} atoms in total.\n".format(len(Arr)))
    if rep!=0:
        Arr=repeat(Arr, disy, rep)
        extra="x{}".format(rep+1)
        print("Nanotube repeated correctly {} times, {} atoms in total.\n".format(rep, len(Arr)))
    if not arch_out:
        arch_out= "Nanotube_{}-{}x{}.xyz".format(n,m, extra)    
    with open(arch_out, mode = 'w') as outfile:
        print(len(Arr), file=outfile)
        print("NT (n,m)=({},{}). z={} A, r={} A. T(x,y)=({},{})".format(n, m, disy*(rep+1), radio, x,y), file=outfile)
        for ele in Arr:
            print("{}    {}    {}    {}".format(ele.ele, ele.x, ele.y, ele.z), file=outfile)
        print("alat", file=outfile)
        print("1.0", file=outfile)
        print("supercell", file=outfile)
        print("{}\t0.0\t0.0".format(2*radio), file=outfile)
        print("0.0\t{}\t0.0".format(2*radio), file=outfile)
        print("0.0\t0.0\t{}".format(disy*(rep+1)), file=outfile)
        print("cartesian coordinates", file=outfile)
    if grf:
        testgraph(Arr)

def coordinatesNT_VASP(n, m, x, y, rep=0,prt=1, arch_out=False, grf=False):
    '''
    Does all the steps necessary to get a nanotube with chiral indices (n,m)
    and translational indices (x,y), repeated along the NT axis "repeat" times.
    If "prt"=1, it prints the completed steps. It automatically prints to file
    "arch_out", defaults to "./POSCAR". Output file is POSCAR VASP
    format.

    Parameters
    ----------
    n: int
        The first chiral index
    m: int
        The second chiral index
    x: int
        The first translational index
    y: int
        The second translational index

    Optional
    --------
    rep: int
        Number of times the nanotube is repeated along its axis.
        Defaults to 0 (no repetition)
    prt: int
        If set to 1 then it prints each completed step in the process
    arch_out: 
        Name of the output file, if False, it defaults to "./Nanotbe_n-m.xyz"
    grf: int
        If set to 1 then it also graphs the result using testgraph

    See Also
    --------
    coordinates_VASP: Does the same but only makes the nanoribbon
    robtainxy: Gives the chiral angle and a list of
               possible x,y pairs from a specific n,m pair
    repeat: Repeats the array N times along the NT axis
    testgraph: VERY basic visualization for testing purposes using matplotlib.

    Examples
    --------
    >>> n=5
    >>> m=3
    >>> theta, Res = robtainxy(n,m)
    >>> x, y, error, p_error=Res[0]
    >>> x=round(x)
    >>> y=round(y)
    >>> coordinatesNT_VASP(n,m,x,y, arch_out="./OutputFile.xyz")
    '''
    theta, Res = robtainxy(n, m)
    if prt==1:
        print("\nChiral angle obtained correctly: Theta=", math.degrees(theta), "\n")
    Arr=arr_initial(n, m, x, y)
    if prt==1:
        print("Initial array correctly created\n")
    Arr= rotate(Arr, theta)
    if prt==1:
        print("Rotation succesfully completed\n")
    Arr, disy, disx = eliminate(Arr, n, m, x, y)
    if prt==1:
        print("Ribbon created correctly, {} atoms in total.\n".format(len(Arr)))
    Arr, radio=nanotube(Arr, n, m)
    for ele in Arr:
        displaced = point3d(0,0,-disy,'')+ ele
        for ele2 in Arr:
            if distance(displaced, ele2)<0.2:
                Arr.remove(ele)
    if prt==1:
        print("Nanotube created correctly, {} atoms in total.\n".format(len(Arr)))
    if rep!=0:
        Arr=repeat(Arr, disy, rep)
        print("Nanotube repeated correctly {} times, {} atoms in total.\n".format(rep, len(Arr)))
    if not arch_out:
        arch_out= "POSCAR"    
    with open(arch_out, mode = 'w') as outfile:
        print("NT (n,m)=({},{}). T(x,y)=({},{})".format(n, m, x,y), file=outfile)
        print(1, file=outfile)
        print("{} 0.0 0.0".format(2*radio), file=outfile)
        print("0.0 {} 0.0".format(2*radio), file=outfile)
        print("0.0 0.0 {}".format(disy*(rep+1)), file=outfile)
        unique_elements = []
        for atom in Arr:
            if atom.ele not in unique_elements:
                unique_elements.append(atom.ele)
        Arr_divided=[[atom for atom in Arr if atom.ele==element] for element in unique_elements]
        ele_names=""
        ele_numbers=""
        for i in range(len(unique_elements)):
            ele_names=ele_names+unique_elements[i]+' '
            ele_numbers = ele_numbers + str(len(Arr_divided[i]))+' '
        print(ele_names.strip(), file=outfile)
        print(ele_numbers.strip(), file=outfile)
        print('Cartesian', file=outfile)
        for ele in Arr_divided:
            for atom in ele:
                print("{} {} {}".format(atom.x, atom.y, atom.z), file=outfile)
    if grf:
        testgraph(Arr)
        
def repeat(Arr, disy, times):
    '''
    Repeats the array in the z direction "times" times. The height of
    the array is "disy". Returns the new array.

    Parameters
    ----------
    Arr: list of point3d
        Array of atoms to be repeated
    disy: float
        Height of the array (length in the z direction)
    times: int
        Times the array will be repeated (0 means no repetition, the array is
        unchanged)
        
    Returns
    -------
    list of point3d
        New array of atoms that has been repeated the corresponding amount
        of times.

    Examples
    --------
    >>> n=5
    >>> m=7
    >>> theta, Res = robtainxy(n,m)
    >>> x, y, error, p_error=Res[0]
    >>> x=round(x)
    >>> y=round(y)
    >>> Arr=arr_initial(n, m, x, y)
    >>> Arr= rotate(Arr, theta)
    >>> Arr, disy, disx = eliminate(Arr, n, m, x, y)
    >>> Arr, radio=nanotube(Arr, n, m)
    >>> rep=2
    >>> Arr=repeat(Arr, disy, rep)
    '''
    Copia=[]
    for ele in Arr:
        Copia.append(ele)
    py=0
    for _ in range(times):
        for ele in Copia:
            Arr.append(point3d(0,0,py+disy,'na')+ele)
        py=py+disy
    return Arr


def searchnum(limit, top=10, error=0.025, r=0):
    '''
    Searches for achiral nanotubes with a number of atoms <= "limit", up
    to n,m < top. With an error = "error". If r=0, it prints the results,
    if r=1 it returns an array with the results.

    Parameters
    ----------
    limit: int
        Max number of atoms that will be searched for
    top: int
        Maximum value for n and m while searching for NTs. Defaults to 10
    error: float
        Error used for getting the x, y pairs. Defaults to 0.025.
    r: int
        If r=0, results are printed. If r=1, results are returned.

    Returns
    -------
    nRes: list of float
        List containing the results: No. of atoms, n,m indices,
        x,y indices, error and percentage error.
    '''
    nRes=[]
    for n in range(1,top):
        for m in range(1,top):
            t, Res = robtainxy(n,m, error=error)
            x,y,err,perr=Res[0]
            x=round(x)
            y=round(y)
            Arr, disy,disx,rad =do_everything(n,m,x,y)
            if len(Arr)<=limit:
                nRes.append((len(Arr),n,m,x,y,err,perr))
    if r==0:
        print("no., n, m, X,    Y,    Error, Percentage error")
        for ele in nRes:
            num, n, m, x, y, err, por=ele
            print("{}, {}, {}, {}, {}, {}, {}".format(num, n, m, x,y, err, por))
    if r==1:
        return nRes


def printcoords(Arr, file_name, comment):
    '''
    Prints the coordinates of an array of point3d to a file "file_name",
    with a comment line "comment". Output is in a simple .xyz format.
    Mainly for testing purposes.

    Parameters
    ----------
    Arr: list of point3d
        Array of atoms to print in a file
    file_name: string
        Name of output file
    comment: string
        Comment to write in a designated comment line on the output file.
    '''
    with open(file_name, mode = 'w') as outfile:
        print(len(Arr), file=outfile)
        print(comment, file=outfile)
        for ele in Arr:
            print("{}    {}    {}    {}".format(ele.ele, ele.x, ele.y, ele.z), file=outfile)


def read_arch(arch, filetype=False):
    '''
    Reads the input file in .in, .xyz (special), or POSCAR VASP formats.
    Returns the number of atoms in the unit cell, the three unit vectors in
    an array and the atoms in the unit cell as an array.

    Parameters
    ----------
    arch: string
        Name of input file
    filetype:
        If set to 'VASP', it will read the file as a POSCAR VASP file.
        Otherwise it detects if its a .xyz or a .in file

    Returns
    -------
    nat: int
        Number of atoms in unit cell
    A: list of point3d
        List with the three unit vectors of the unit cell
    UnitCell: list of point3d
        List containing the atoms with their respective positions in the
        unit cell.

    Examples
    --------
    >>> nat, A, UnitCell = read_arch("./MoS2UnitCell.in")

    >>> nat, A, UnitCell = read_arch("./POSCAR", filetype='VASP')
    '''
    A=[]
    UnitCell=[]
    if filetype=='VASP':
        with open(arch, mode='r') as infile:
            infile.readline()
            scale = float(infile.readline().strip())
            for _ in range(3):
                valores= infile.readline().strip().split()
                A.append(point3d(float(valores[0])*scale, float(valores[1])*scale, float(valores[2])*scale, 'Vector'))
            line = infile.readline().strip()
            element_numbers = []
            try:
                for number in line.split():
                    tmp = int(number)
                print("\n***ERROR*** No line with species names. Please specify the species.\n")
                sys.exit()
            except ValueError:
                element_names = line.split()
            for number in infile.readline().strip().split():
                element_numbers.append(int(number))
            nat=sum(element_numbers)
            line=infile.readline().strip()
            if (line.startswith("s") or line.startswith("S")):
                line = infile.readline().strip()
            if  line == 'Direct' or line=='direct':
                for ele_no in range(len(element_names)):
                    for _ in range(element_numbers[ele_no]):
                        valores = infile.readline().strip().split()
                        position = a1*valores[0] + a2*valores[1] + a3*valores[2]
                        UnitCell.append(point3d(position.x, position.y, position.z, element_names[ele_no]))
            elif line == 'Cartesian' or line=='cartesian':
                for ele_no in range(len(element_names)):
                    for _ in range(element_numbers[ele_no]):
                        valores = infile.readline().strip().split()
                        UnitCell.append(point3d(float(valores[0]), float(valores[1]), float(valores[2]), element_names[ele_no]))     
            else:
                print("\n***ERROR*** Coordinates can only be in 'Direct' format or 'Cartesian' format.\n {} format is not recognized\n".format(line))
                sys.exit()
            return nat, A, UnitCell
                    
    if arch.endswith(".in"):
        with open(arch, mode='r') as infile:
            for line in infile:
                line=line.strip()
                if line.startswith("nat"):
                    nat= line.replace("nat","").replace("=","")
                    nat=int(nat.strip().strip(','))
                elif line.startswith("CELL_PARAMETERS"):
                    for _ in range(3):
                        valores= infile.readline().strip().split()
                        A.append(point3d(float(valores[0]), float(valores[1]), float(valores[2]), 'Vector'))
                elif line.startswith("ATOMIC_POSITIONS"):
                    units = line.replace("ATOMIC_POSITIONS", "").strip().strip('{}()')
                    if units=='angstrom' or units=='Angstrom':
                        for _ in range(nat):
                            valores = infile.readline().strip().split()
                            UnitCell.append(point3d(float(valores[1]), float(valores[2]), float(valores[3]), valores[0]))
                    elif units=='crystal' or units =='Crystal':
                        print("\n***ERROR*** The use of input files with 'crystal' units has not been implemented. Please convert the units to angstroms.\n")
                        sys.exit()
                    else:
                        print("\n***ERROR*** The units '{}' are not known or have not been implemented.".format(units))
                        sys.exit()
    elif arch.endswith(".xyz"):
        with open(arch, mode='r') as infile:
            for line in infile:
                if line.strip().endswith("coordinates"):
                    units = line.replace("coordinates", "").strip()
                    if units == "cartesian" or units== "Cartesian":
                        break
                    else:
                        print("\n***ERROR*** The units '{}' are not known or have not been implemented. Please convert to cartesian coordinates.".format(units))
                        sys.exit()
            infile.seek(0)
            nat=int(infile.readline().strip())
            line=infile.readline()
            for _ in range(nat):
                valores = infile.readline().strip().split()
                UnitCell.append(point3d(float(valores[1]), float(valores[2]), float(valores[3]), valores[0]))
            for line in infile:
                if line.startswith("supercell"):
                    for _ in range(3):
                        valores= infile.readline().strip().split()
                        A.append(point3d(float(valores[0]), float(valores[1]), float(valores[2]), 'Vector'))
    else:
        with open(arch, mode='r') as infile:
            print("***ERROR*** The use of .{} files has not been implemented.".format(arch.split(".")[1]))
    return nat, A, UnitCell

def shift_cell(A, UnitCell):
    '''
    Shifts the unit cell so that it is centered around the a_3/2 plane.

    Parameters
    ----------
    A: list of point3d
        List containing the unit vectors of the cell.
    UnitCell: list of point3d
        List containing all atoms in the unit cell.

    Returns
    -------
    UnitCell: list of point3d
        List containing all atoms in the unit cell shifted to be centered
        around the a3/2 plane.
    '''
    Alturas=[]
    for ele in UnitCell:
        if ele.z not in Alturas:
            Alturas.append(ele.z)
    rango= max(Alturas)-min(Alturas)
    #prom=stat.mean(Alturas)
    punto_medio=min(Alturas)+rango/2
    altura_tot=A[2].z
    medio= altura_tot/2
    #distance=medio-prom
    diferencia=medio-punto_medio
    for ele in UnitCell:
        ele.z=ele.z+diferencia
    for ele in UnitCell:
        if ele.z>altura_tot or ele.z<0:
            print("\n***ERROR*** Problems with the z axis of the Unit Cell (a3 Vector might be too small).")
            sys.exit()
    return UnitCell

def sortfunc(ele):
    '''
    Returns cartesian distance between two elements.

    Parameters
    ----------
    ele: list of float
        Tuple containing four elements, the first two are x and y which will
        define the distance

    Returns
    -------
    float
        Cartesian distance between x and y

    See Also
    --------
    sortfuncT: Returns distance scaled by the unit cell vectors
    
    '''
    x, y, err, perr=ele
    return math.sqrt(x*x + y*y)

def sortfuncT(ele):
    '''
    Returns distance between two elements, scaled by the unit cell vectors.

    Parameters
    ----------
    ele: list of float
        Tuple containing four elements, the first two are x and y which will
        define the distance

    Returns
    -------
    float
        Distance scaled by unit cell vectors between x and y.

    See Also
    --------
    sortfunc: Returns distance between two elements.
    
    '''
    x,y,err,perr=ele
    return distance(point3d(0,0,0,0), a1*x+ a2*y)

def rsearchnm(radius, error=0.1, MAX=100):
    '''
    Looks for nanotubes with radius close to "radius" up to an
    error "error" which defaults to 0.1. Returns an array with
    the results. Goes up to n,m < MAX.

    Parameters
    ----------
    radius: float
        objective radius to look for
    error: float
        Allowed error when looking for nanotube with the desired radius.
        Defauts to 0.1
    MAX: int
        Defaults to 100. Maximum value of n and m when looking for NTs.

    Returns
    -------
    list of floats
        List with results. Contains the n,m indices, the error and the
        relative error        
    '''
    Res = []
    for n in range(MAX):
        for m in range(MAX):
            rad_prueba= distance(point3d(0,0,0,0), (a1*n+a2*m))/(2*math.pi)
            if abs(rad_prueba-radius)<error:
                Res.append((n, m, abs(rad_prueba-radius), abs(rad_prueba-radius)/radius))
    Res.sort(key=sortfuncT)    
    return Res

def searchnm(radius, error=0.1, MAX=100):
    '''
    Looks for nanotubes with radius close to "radius" up to an
    error "error" which defaults to 0.1. Prints an array with
    the results. Goes up to n,m < MAX. 

    Parameters
    ----------
    radius: float
        objective radius to look for
    error: float
        Allowed error when looking for nanotube with the desired radius.
        Defauts to 0.1
    MAX: int
        Defaults to 100. Maximum value of n and m when looking for NTs.    
    '''
    Res = rsearchnm(radius, error=error, MAX=MAX)
    if len(Res)==0:
        print("No nanotube found with a radius {} +- {}. Try increasing the accepted error.".format(radius, error))
    else:
        print("Results for radius={} +- {}\n\nn,m,\tError,\tRelative error".format(radius, error))
        for ele in Res:
            print("{},{},\t{},{}".format(ele[0], ele[1], round(ele[2], 6), round(ele[3],6)))
    
def do_everything(n,m,x,y, nt=True):
    '''
    For testing purposes. If nt = True, makes a nanotube with C(n,m)
    and T(x,y). If nt=False, makes only the nanoribbon. Returns the array
    of atoms, the height of the NT, its circumference, and its radius.

    Parameters
    ----------
    n: int
        The first chiral index
    m: int
        The second chiral index
    x: int
        The first translational index
    y: int
        The second translational index
    nt: boolean
        If False, only makes the nanoribbon. If True, it makes the whole NT.

    Returns
    -------
    Arr: list of point3d
        Array containing the Nanotube or nanoribbon atoms.
    disy: float
        Height of the Nanotube / nanoribbon.
    disx: float
        Length of nanoribbon (or circumference of NT)
    radio: float
        Radius of the NT. Defaults to 0 if the user made a nanoribbon.
    '''
    theta, Res = robtainxy(n,m)
    Arr, disy, disx = eliminate(rotate(arr_initial(n,m,x,y), theta), n,m,x,y)
    radio=0
    if nt:
        Arr, radio = nanotube(Arr, n,m)
    return Arr, disy, disx, radio       
    
def main():
    ribbon=False
    r=0
    Error=ERROR
    
    if "-nr" in opts:
        ribbon=True
        
    for word in opts:
        if word.startswith("-E"):
            Error= float(word.strip().replace("-E", ""))
            
    if "-radius" in opts:
        rad=float(input("\nDesired radius: "))
        if Error == ERROR:
            searchnm(rad)
        else:
            searchnm(rad, error=Error)
        sys.exit()
        
    if "-diameter" in opts:
        rad=float(input("\nDesired diameter: "))
        rad=rad/2
        if Error == ERROR:
            searchnm(rad)
        else:
            searchnm(rad, error=Error/2)
        sys.exit()
        
    for word in opts:
        if word.startswith("-r"):
            r= int(word.strip().replace("-r", ""))
            
    n= int(input("\nValue of n: "))
    m= int(input("\nValue of m: "))
    
    if "-o" in opts:
        obtainxy(n,m,error=Error)
        sys.exit()
        
    try:
        arch_out=args[1]
    except IndexError:
        print("\n***No output file specified. \n\t A generic name will be used.***\n\n")
        arch_out=False
    
    theta, Res = robtainxy(n, m, error=Error)
    try:
        x, y, err, perr=Res[0]
    except IndexError:
        print("\n***ERROR*** No integer solutions found for translational indices (x,y) in the range x,y<100 \nPlease try another pair of chiral indices.")
        sys.exit
    x= round(x)
    y= round(y)

    if "-VASPout" in opts:
        if ribbon:
            coordinates_VASP(n,m,x,y,arch_out=arch_out, grf=False)
        else:
            coordinatesNT_VASP(n,m,x,y,rep=r, arch_out=arch_out, grf=False)
    else:
        if ribbon:
            coordinates_xyz(n,m,x,y,arch_out=arch_out, grf=False)
        else:
            coordinatesNT_xyz(n,m,x,y,rep=r, arch_out=arch_out, grf=False)
        



if __name__== "__main__":
    if ("-h" in opts) or ("--h" in opts) or ("-help" in opts):
        print("\nChiraltube prints atomic coordinates for different structures in special .xyz files for their usage in other visualization, simulation or calculation software.")
        print("This source code is released under the 3-Clause BSD License, see 'LICENSE.txt'")
        print("Created by José M. de Albornoz-Caratozzolo and Felipe Cervantes-Sodi.")
        print("\nUsage: \n\t\t python3 ~/chiraltube.py <input file> [<output file>] [<options>]")
        print("\n If an output file is not specified, an automatically generated name will be used: './Ribbon_n-m.xyz' or './Nanotube_n-m.xyz', accordingly.")
        print("The input file must contain the unit cell from which you want to create a nanotube or nanoribbon. It must be in special .xyz format or in the .in format used for Quantum Espresso.")
        print("NOTE: For the moment, Chiraltube only accepts cartesian coordinates in angstroms.")
        print("\nOPTIONS: (by order of program execution)")
        print("\tOption \t\tResult\n")
        print("\t-h / -help / --h \tPrints this help message\n")
        print("\t-VASPin \tDenotes that the input file is in the POSCAR VASP format.\n")
        print("\t-p \t\tPrints atomic positions and unit vectors of the unit cell.\n\t\t\tUseful for checking if the program read \n\t\t\tthe input file correctly.\n")
        print("\t-nr \t\tGives coordinates for the nanoribbon instead of the nanotube.\n")        
        print("\t-E<num> \tChanges the acceptable error when looking for\n\t\t\tinteger solutions to the translational indices (x,y)\n\t\t\tExample: '-E0.05' changes the acceptable error to 0.05.\n\t\t\tAlso changes the error if looking for nanotubes by radius.\n")
        print("\t-radius \tLooks for nanotubes that have a certain radius.\n\t\t\tThe program will ask for the radius desired.\n")
        print("\t-diameter \tSame as above but with diameter.")
        print("\t-r<num> \tRepeats the nanotube <num> times along its axis.\n\t\t\t<num> must be an integer.\n\t\t\tExample: '-r3' would yield a nanotube repeated 3 times.\n")
        print("\t-o \t\tPrints the (x,y) pairs that are solutions for\n\t\t\tthe (n,m) provided.\n")
        print("\t-VASPout \tGives the output file in the POSCAR VASP format.\n")
        sys.exit()
    try:
        arch_in=args[0]
    except IndexError:
        print("\n***ERROR*** No input file specified.\nType 'python3 chiraltube.py -h' for more help.")
        sys.exit()
    filetype=False
    if "-VASPin" in opts:
        filetype="VASP"
    try: nat, A, UnitCell =read_arch(arch_in, filetype=filetype)
    except FileNotFoundError:
        print("\n***ERROR*** Couldn't find the file '{}'\n".format(arch_in))
        sys.exit()
    UnitCell=shift_cell(A, UnitCell)
    print("\nUnit cell read succesfully from input file {}\n".format(arch_in))
    if "-p" in opts:    
        print("\nnat=", nat)
        print("\nCell Parameters:\n a1={}\n a2={}\n a3={}".format(A[0], A[1], A[2]))
        print("\nUnit Cell:")
        for ele in UnitCell:
            print("", ele)
        
        
    a1= A[0]
    a2= A[1]
    a3= A[2]
    a2= a2*(-1) 
    for i in range(len(UnitCell)):
        UnitCell[i]=a2+UnitCell[i] #done for the program calculations
    phi= math.acos((a1.x*a2.x+a1.y*a2.y+a1.z*a2.z)/(a1.mag()*a2.mag()))


    if "-p" in opts:
        sys.exit()
    main()

    
        
