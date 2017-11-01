from subprocess import call
pi = 3.14159265359
call(["ls", "-l"])
import glob
from Tkinter import *
#import gui
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import random
import sys
import os
import numpy as np
import matplotlib.cm as cm
#import tkinter as tk
import tkFileDialog as filedialog
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
# implement the default mpl key bindings
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
#import enthought.mayavi.mlab as mylab  #bead models
import sys
import plotly as py
py.tools.set_credentials_file(username='Molodenskiy', api_key='Jc4uP32ntpjUyxtMuWhl')
import plotly.plotly as py
import plotly.tools as plotly_tools
from plotly.graph_objs import *
from mpl_toolkits.mplot3d import proj3d
if sys.version_info[0] < 3:
    import Tkinter as Tk
else:
    import tkinter as Tk
import time

class Arrow3D(FancyArrowPatch):

    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        FancyArrowPatch.draw(self, renderer)

start = time.time()

# Function make ATOM structure from a pdb string
def mk_atom(pdbstring):

    atomserialnumber = pdbstring[7:11]
    atomname = pdbstring[13:16]
    alternativelocationindicator = pdbstring[17]
    residuename = pdbstring[18:20]
    chainidentifier = pdbstring[21]
    residuesequencenumber = pdbstring[23:26]
    codeforinsectionsofresidues = pdbstring[27]
    xcoord = pdbstring[31:38]
    ycoord = pdbstring[39:46]
    zcoord = pdbstring[47:54]
    occupancy = pdbstring[55:60]
    tempfactor = pdbstring[61:66]
    segmentidentifier = pdbstring[73:76]
    elementsymbol = pdbstring[77:78]
    return (ATOM(pdbstring, atomserialnumber, atomname,alternativelocationindicator,residuename,chainidentifier,residuesequencenumber,
                 codeforinsectionsofresidues,xcoord,ycoord,zcoord,occupancy,tempfactor,segmentidentifier,elementsymbol))

# File save dialog
def file_save(text):
    f = Tk.filedialog.asksaveasfile(mode='w', defaultextension=".txt")
    if f is None: # asksaveasfile return `None` if dialog closed with "cancel".
        return
    f.write(text)
    f.close()
# File load dialog
def file_load():
    file_path = filedialog.askopenfilename()
    with open(file_path) as f:
        pdblist = f.readlines()
        return pdblist

#class ATOM contains all info about the atom from pdb string (starting with ATOM) and returns each value by response
class ATOM:
    atomserialnumber = -1
    atomname = ""
    alternativelocationindicator = ""
    residuename = ""
    chainidentifier = ""
    residuesequencenumber = -1
    codeforinsectionsofresidues = ""
    xcoord = -1
    ycoord = -1
    zcoord = -1
    occupancy = -1
    tempfactor = -1
    segmentidentifier = ""
    elementsymbol = ""
    #backup
    __initstring = ""
    weights = {'H' : 1.008,'C' : 12.011,'N' : 14.007,'O' : 15.999,'S' : 32.060}
    __weight = 1

    def __init__(self,initstring, atomserialnumber, atomname,alternativelocationindicator,residuename,chainidentifier,residuesequencenumber,
                 codeforinsectionsofresidues,xcoord,ycoord,zcoord,occupancy,tempfactor,segmentidentifier,elementsymbol):
        self._initstring = initstring
        self._atomserialnumber = atomserialnumber
        self._atomname = atomname
        self._alternativelocationindicator = alternativelocationindicator
        self._residuename = residuename
        self._chainidentifier = chainidentifier
        self._residuesequencenumber = residuesequencenumber
        self._codeforinsectionsofresidues = codeforinsectionsofresidues
        self._xcoord = xcoord
        self._ycoord = ycoord
        self._zcoord = zcoord
        self._occupancy = occupancy
        self._tempfactor = tempfactor
        self._segmentidentifier = segmentidentifier
        self._elementsymbol = elementsymbol
        self._weight = self.weights[atomname[0]]

    #simple overload
    #def __init__(self,initstring, element,  residuecount, xcoord, ycoord, zcoord, occupancy):
    #    self.__initstring = initstring
    #    self.__element = element
    #    #self._chain = chain
    #    self._atomcount = residuecount
    #    self.__xcoord = xcoord
    #    self.__ycoord = ycoord
    #    self.__zcoord = zcoord
    #    self.__occupancy = occupancy
    #    if element[0] != ('H' or 'C' or 'N' or 'O' or 'S'):
    #        print("Bad atom! ")
    #        print (initstring)
    #    self.__weight = self.weights[element[0]]
    #def set_name(self, name):
    #    self.__name = name


    def get_pdb_string(self):
        return self._initstring

    def get_atomserialnumber(self):
        return self._atomserialnumber

    def get_atomname(self):
        return self._atomname

    def get_alternativelocationindicator(self):
        return self._alternativelocationindicator

    def get_residuename(self):
        return self._residuename

    def get_chainidentifier(self):
        return self._chainidentifier

    def get_residuesequencenumber(self):
        return int(self._residuesequencenumber)

    def get_codeforinsectionsofresidues(self):
        return self._codeforinsectionsofresidues

    def get_x_coord(self):
        return float(self._xcoord)

    def get_y_coord(self):
        return float(self._ycoord)

    def get_z_coord(self):
        return float(self._zcoord)

    def get_occupancy(self):
        return float(self._occupancy)

    def get_temp_factor(self):
        return float(self._tempfactor)


    def get_segmentidentifier(self):
        return self._segmentidentifier

    def get_elementsymbol(self):
        return self._elementsymbol

    def get_weight(self):
        return float(self._weight)

    #polymorphism?
    def get_type(self):
        return("ATOM")

    def file_save(self):
        f = tkFileDialog.asksaveasfile(mode='w', defaultextension=".txt")
        if f is None:  # asksaveasfile return `None` if dialog closed with "cancel".
            return
        text2save = str(text.get(1.0, END))  # starts from `1.0`, not `0.0`
        f.write(text2save)
        f.close()

# draw atoms in centre of mass and draw principal axes of inertia
def draw_atoms(atoms, eigenVectors, centerofmass):
    x = []
    y = []
    z = []
    Nm = []
    colors = []
    colors_dic = {'H' : 'white', 'C' : 'black', 'N' : 'blue', 'O': 'red', 'S': 'yellow'}
    for a in atoms:
        x.append(float(a.get_x_coord()) - centerofmass[0])
        y.append(float(a.get_y_coord()) - centerofmass[1])
        z.append(float(a.get_z_coord()) - centerofmass[2])
        Nm.append(a.get_atomname()[0])
        colors.append(colors_dic[a.get_atomname()[0]])
    # bead models (alternative)
    #mylab.points3d(x, y, z, 20)
    #mylab.show()
    length = max(x) - min(x)
    if (max(y) - min (y))> length: length = max(y) - min(y)
    if (max(z) - min(z)) > length: length = max(z) - min(z)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x, y, z, c=[ colors_dic[x] for x in Nm], s=5)
    ax.set_xlabel('X, Angstrem')
    ax.set_ylabel('Y, Angstrem')
    ax.set_zlabel('Z, Angstrem')

    eigenVectors = np.transpose(eigenVectors) * 0.5 * length

    for vec in eigenVectors:
        drawvec = Arrow3D([0, vec[0, 0]], [0, vec[0, 1]], [0, vec[0, 2]],
                      mutation_scale=20, lw=3, arrowstyle="-|>", color='r')
    # adding the arrow to the plot
        ax.add_artist(drawvec)
    return

def move_to_center_mass(atoms, centerofmass):
    a = []
    for r in range(len(atoms)):  # len(atoms) raws
        a.append([])  # create empty raw
        a[r].append(float(atoms[r].get_x_coord()) - centerofmass[0])
        a[r].append(float(atoms[r].get_y_coord()) - centerofmass[1])
        a[r].append(float(atoms[r].get_z_coord()) - centerofmass[2])
    a = np.array(a)
    return a

#calculation of weight center in certain range (in residues)
def center_mass(atoms):
    CMx = 0
    M = 0
    CMy = 0
    CMz = 0
    position = []
    for atom in atoms:
        x_coord = atom.get_x_coord()
        y_coord = atom.get_y_coord()
        z_coord = atom.get_z_coord()
        weight = atom.get_weight()
        CMx += x_coord * weight
        CMy += y_coord * weight
        CMz += z_coord * weight
        M += weight
    position.append(CMx / M)
    position.append(CMy / M)
    position.append(CMz / M)
    return position

def inertia_tensor_calc(atoms):
    Jxx = 0
    Jyy = 0
    Jzz = 0
    Jxy = 0
    Jyz = 0
    Jxz = 0
    for atom in atoms:
        weight = atom.get_weight()
        x = atom.get_x_coord()
        y = atom.get_y_coord()
        z = atom.get_z_coord()
        Jxx += weight*(y*y + z*z)
        Jyy += weight * (x * x + z * z)
        Jzz += weight * (x * x + y * y)
        Jxy += weight * (-x*y)
        Jxz += weight * (-x*z)
        Jyz += weight * (-y*z)
    tens_inertia = np.matrix([[Jxx, Jxy, Jxz], [Jxy, Jyy, Jyz], [Jxz, Jyz, Jzz]])
    return tens_inertia

def inertia_tensor_calc_matrix(matrix, weight):
    num_rows, num_cols = matrix.shape
    x = []
    y = []
    z = []
    for a in range(num_rows):
        x.append(matrix[a,0])
        y.append(matrix[a,1])
        z.append(matrix[a,2])

    Jxx = 0
    Jyy = 0
    Jzz = 0
    Jxy = 0
    Jyz = 0
    Jxz = 0
    for i in range(num_rows):
        Jxx += weight[i]*(y[i]*y[i] + z[i]*z[i])
        Jyy += weight[i] * (x[i] * x[i] + z[i] * z[i])
        Jzz += weight[i] * (x[i] * x[i] + y[i] * y[i])
        Jxy += weight[i] * (-x[i]*y[i])
        Jxz += weight[i] * (-x[i]*z[i])
        Jyz += weight[i] * (-y[i]*z[i])
    tens_inertia = np.matrix([[Jxx, Jxy, Jxz], [Jxy, Jyy, Jyz], [Jxz, Jyz, Jzz]])
    return tens_inertia

def rmsd(V, W):
    """
    Calculate Root-mean-square deviation from two sets of vectors V and W.

    Parameters
    ----------

    V : array
        (N,D) matrix, where N is points and D is dimension.
    W : array
        (N,D) matrix, where N is points and D is dimension.

    Returns
    -------
    rmsd : float
        Root-mean-square deviation

    """
    D = len(V[0])
    N = len(V)
    M = len(W)
    if M < N: N = M
    rmsd = 0.0
    for v, w in zip(V, W):
        rmsd += sum([(v[i] - w[i])**2.0 for i in range(D)])
    return np.sqrt(rmsd/(N))

def kabsch_rotate(P, Q):
    """
    Rotate matrix P unto matrix Q using Kabsch algorithm.

    Parameters
    ----------
    P : array
        (N,D) matrix, where N is points and D is dimension.
    Q : array
        (N,D) matrix, where N is points and D is dimension.

    Returns
    -------
    P : array
        (N,D) matrix, where N is points and D is dimension,
        rotated

    """
    U = kabsch(P, Q)

    # Rotate P
    P = np.dot(P, U)
    return P

def kabsch(P, Q):
    """
    The optimal rotation matrix U is calculated and then used to rotate matrix
    P unto matrix Q so the minimum root-mean-square deviation (RMSD) can be
    calculated.

    Using the Kabsch algorithm with two sets of paired point P and Q, centered
    around the center-of-mass. Each vector set is represented as an NxD
    matrix, where D is the the dimension of the space.

    The algorithm works in three steps:
    - a translation of P and Q
    - the computation of a covariance matrix C
    - computation of the optimal rotation matrix U

    http://en.wikipedia.org/wiki/Kabsch_algorithm

    Parameters
    ----------
    P : array
        (N,D) matrix, where N is points and D is dimension.
    Q : array
        (N,D) matrix, where N is points and D is dimension.

    Returns
    -------
    U : matrix
        Rotation matrix (D,D)

    """

    # Computation of the covariance matrix
    C = np.dot(np.transpose(P), Q)

    # Computation of the optimal rotation matrix
    # This can be done using singular value decomposition (SVD)
    # Getting the sign of the det(V)*(W) to decide
    # whether we need to correct our rotation matrix to ensure a
    # right-handed coordinate system.
    # And finally calculating the optimal rotation matrix U
    # see http://en.wikipedia.org/wiki/Kabsch_algorithm
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]

    # Create Rotation matrix U
    U = np.dot(V, W)

    return U

def surface_plot(X,Y,Z,**kwargs):
    """ WRITE DOCUMENTATION
    """
    xlabel, ylabel, zlabel, title = kwargs.get('X, Ang',""), kwargs.get('Y Ang',""), kwargs.get('Z Ang',""), kwargs.get('Molecule title',"")
    fig = plt.figure()
    fig.patch.set_facecolor('white')
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X,Y,Z)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    ax.set_title(title)
    plt.show()
    plt.close()

def draw_and_calc(atoms, start_residue, end_residue):
    newatoms = []
    weights = []
    #subdivide atoms
    for atom in atoms:
        current_rnumb = atom.get_residuesequencenumber()
        if (start_residue <= current_rnumb <  end_residue):
            newatoms.append(atom)
            weights.append(atom.get_weight())
    #position of center mass
    cm_position = center_mass(newatoms)
    #calculation of inertia tensor
    i_tensor = inertia_tensor_calc(newatoms)
    # calculating eigen values and eigen vectors
    eigenValues, eigenVectors = np.linalg.eig(i_tensor)
    # sorting
    idx = eigenValues.argsort()[::-1]
    eigenValues = eigenValues[idx]
    eigenVectors = eigenVectors[:, idx]
    #draw subplot in center of mass and eigen vectors
    draw_atoms(newatoms, eigenVectors, cm_position)
    # A in format (x1,y1,z1),...,(xn,yn,zn)
    A = move_to_center_mass(newatoms, cm_position)
    return A, weights

def draw_atoms_from_matrix(A):
    num_rows, num_cols = A.shape
    x = []
    y = []
    z = []
    for a in range(num_rows):
        x.append(A[a,0])
        y.append(A[a,1])
        z.append(A[a,2])

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x, y, z, c= 'black', s=5)

    ax.set_xlabel('X, Angstrem')
    ax.set_ylabel('Y, Angstrem')
    ax.set_zlabel('Z, Angstrem')
    #plt.show()
    return

def draw_atoms_from_2matrix(A, weights1, B, weights2):
    num_rows1, num_cols1 = A.shape
    num_rows2, num_cols2 = B.shape
    JA = inertia_tensor_calc_matrix(A, weights1)
    JB = inertia_tensor_calc_matrix(B, weights2)
    #finding eigen values and eigen vectors
    eigenValues1, eigenVectors1 = np.linalg.eig(JA)
    eigenValues2, eigenVectors2 = np.linalg.eig(JB)

    # sorting
    idx1 = eigenValues1.argsort()[::-1]
    eigenValues1 = eigenValues1[idx1]
    eigenVectors1 = eigenVectors1[:, idx1]
    idx2 = eigenValues2.argsort()[::-1]
    eigenValues2 = eigenValues2[idx2]
    eigenVectors2 = eigenVectors2[:, idx2]

    x1 = []
    y1 = []
    z1 = []
    x2 = []
    y2 = []
    z2 = []
    for a in range(num_rows1):
        x1.append(A[a,0])
        y1.append(A[a,1])
        z1.append(A[a,2])
    for b in range(num_rows2):
        x2.append(B[b,0])
        y2.append(B[b,1])
        z2.append(B[b,2])
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x1, y1, z1, c= 'red', s=5)
    ax.scatter(x2, y2, z2, c= 'black', s=5)
    ax.set_xlabel('X, Angstrem')
    ax.set_ylabel('Y, Angstrem')
    ax.set_zlabel('Z, Angstrem')
    #draw eigen vectors
    length1 = max(x1) - min(x1)
    if (max(y1) - min(y1)) > length1: length1 = max(y1) - min(y1)
    if (max(z1) - min(z1)) > length1: length1 = max(z1) - min(z1)
    length2 = max(x2) - min(x2)
    if (max(y2) - min(y2)) > length2: length2 = max(y2) - min(y2)
    if (max(z2) - min(z2)) > length2: length2 = max(z2) - min(z2)

    eigenVectors1=np.transpose(eigenVectors1)*0.5*length1
    eigenVectors2 = np.transpose(eigenVectors2)*0.5 * length2

    for vec in eigenVectors1:
        drawvec = Arrow3D([0, vec[0,0]], [0, vec[0,1]], [0, vec[0,2]],
                          mutation_scale=20, lw=3, arrowstyle="-|>", color='b')
        # adding the arrow to the plot
        ax.add_artist(drawvec)
        #ax.plot([0, 0.5*length1*a[0,0]], [0, 0.5*length1*a[0,1]], zs=[0, 0.5*length1*a[0,2]], color = 'blue')
    for vec in eigenVectors2:
        #ax.plot([0, 0.5*length2*b[0,0]], [0, 0.5*length2*b[0,1]], zs=[0, 0.5*length2*b[0,2]], color = 'green')
        drawvec = Arrow3D([0, vec[0, 0]], [0, vec[0, 1]], [0, vec[0, 2]],
                          mutation_scale=20, lw=3, arrowstyle="-|>", color='g')
        # adding the arrow to the plot
        ax.add_artist(drawvec)
    #plt.show()
    return

def draw_atoms_from_2matrix_subJ(A, weights1, B, weights2):
    num_rows1, num_cols1 = A.shape
    num_rows2, num_cols2 = B.shape
    JA = inertia_tensor_calc_matrix(A, weights1)
    JB = inertia_tensor_calc_matrix(B, weights2)
    #finding eigen values and eigen vectors
    eigenValues1, eigenVectors1 = np.linalg.eig(JA)
    eigenValues2, eigenVectors2 = np.linalg.eig(JB)

    # sorting
    idx1 = eigenValues1.argsort()[::-1]
    eigenValues1 = eigenValues1[idx1]
    eigenVectors1 = eigenVectors1[:, idx1]
    idx2 = eigenValues2.argsort()[::-1]
    eigenValues2 = eigenValues2[idx2]
    eigenVectors2 = eigenVectors2[:, idx2]
    # superimpose their principal axes
    #A, eigenVectors1 = superimpose_two_inertia_tensors_very_simple(A , eigenVectors1, eigenVectors2)
    #rotate manually
    thickness = 5
    A = rotation_z(A, 0)
    A = rotation_y(A, 90)
    A = rotation_x(A, 0)
    print("RMSD = " + str(rmsd(A,B)))
    x1 = []
    y1 = []
    z1 = []
    x2 = []
    y2 = []
    z2 = []
    for a in range(num_rows1):
        x1.append(A[a,0])
        y1.append(A[a,1])
        z1.append(A[a,2])
    for b in range(num_rows2):
        x2.append(B[b,0])
        y2.append(B[b,1])
        z2.append(B[b,2])
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x1, y1, z1, c= 'red', s=thickness)
    ax.scatter(x2, y2, z2, c= 'black', s=thickness)
    ax.set_xlabel('X, Angstrem')
    ax.set_ylabel('Y, Angstrem')
    ax.set_zlabel('Z, Angstrem')
    #draw eigen vectors
    length1 = max(x1) - min(x1)
    if (max(y1) - min(y1)) > length1: length1 = max(y1) - min(y1)
    if (max(z1) - min(z1)) > length1: length1 = max(z1) - min(z1)
    length2 = max(x2) - min(x2)
    if (max(y2) - min(y2)) > length2: length2 = max(y2) - min(y2)
    if (max(z2) - min(z2)) > length2: length2 = max(z2) - min(z2)

    eigenVectors1=np.transpose(eigenVectors1)*0.5*length1
    eigenVectors2 = np.transpose(eigenVectors2)*0.5 * length2

    for vec in eigenVectors1:
        drawvec = Arrow3D([0, vec[0,0]], [0, vec[0,1]], [0, vec[0,2]],
                          mutation_scale=20, lw=3, arrowstyle="-|>", color='b')
        # adding the arrow to the plot
        ax.add_artist(drawvec)
        #ax.plot([0, 0.5*length1*a[0,0]], [0, 0.5*length1*a[0,1]], zs=[0, 0.5*length1*a[0,2]], color = 'blue')
    for vec in eigenVectors2:
        #ax.plot([0, 0.5*length2*b[0,0]], [0, 0.5*length2*b[0,1]], zs=[0, 0.5*length2*b[0,2]], color = 'green')
        drawvec = Arrow3D([0, vec[0, 0]], [0, vec[0, 1]], [0, vec[0, 2]],
                          mutation_scale=20, lw=3, arrowstyle="-|>", color='g')
        # adding the arrow to the plot
        ax.add_artist(drawvec)
    #plt.show()
    return

def superimpose_two_inertia_tensors(A, eigenVectors1, B, eigenVectors2):
    #superimposing axes using Kabsh algorithm
    ROT = kabsch (eigenVectors1, eigenVectors2)
    ROT2Euler(ROT)
    RMSD_min = rmsd(np.array(eigenVectors1), np.array(eigenVectors2))
    ev1 = eigenVectors1.copy()
    ev1 = ev1[:, [0, 1, 2]]
    RMSD = rmsd(np.array(ev1), np.array(eigenVectors2))
    if RMSD < RMSD_min:
        RMSD_min =  RMSD
        ev1_final = ev1
    ev1 = eigenVectors1.copy()
    ev1 = ev1[:, [0, 2, 1]]
    RMSD = rmsd(np.array(ev1), np.array(eigenVectors2))
    if RMSD < RMSD_min:
        ev1_final = ev1
        RMSD_min =  RMSD
    ev1 = eigenVectors1.copy()
    ev1 = ev1[:, [1, 0, 2]]
    RMSD = rmsd(np.array(ev1), np.array(eigenVectors2))
    if RMSD < RMSD_min:
        ev1_final = ev1
        RMSD_min =  RMSD
    ev1 = eigenVectors1.copy()
    ev1 = ev1[:, [1, 2, 0]]
    RMSD = rmsd(np.array(ev1), np.array(eigenVectors2))
    if RMSD < RMSD_min:
        ev1_final = ev1
        RMSD_min =  RMSD
    ev1 = eigenVectors1.copy()
    ev1 = ev1[:, [2, 1, 0]]
    RMSD = rmsd(np.array(ev1), np.array(eigenVectors2))
    if RMSD < RMSD_min:
        ev1_final = ev1
        RMSD_min =  RMSD
    ev1 = eigenVectors1.copy()
    ev1 = ev1[:, [2, 0, 1]]
    RMSD = rmsd(np.array(ev1), np.array(eigenVectors2))
    if RMSD < RMSD_min:
        ev1_final = ev1
        RMSD_min =  RMSD
    ROT = kabsch (ev1_final, eigenVectors2)
    ROT2Euler(ROT)
    return  np.dot(A, ROT), np.dot(eigenVectors1, ROT)

def superimpose_two_inertia_tensors_simple(A, eigenVectors1, B, eigenVectors2):
    #superimposing axes using Kabsh algorithm
    ROT = kabsch (eigenVectors1, eigenVectors2)
    ROT2Euler(ROT)
    # something wrong with eigenvectors1
    return  np.dot(A, ROT), np.dot(eigenVectors1, ROT)

def superimpose_two_inertia_tensors_very_simple(A, eigenVectors1, eigenVectors2):
    ev1_r = np.transpose(eigenVectors1)
    ev2_r = np.transpose(eigenVectors2)
    #return  np.dot(A, ev1_r), np.dot(eigenVectors1, ev2_r)
    #converting arbitrary rotational matrix into Euler angles
    # something wrong with eigenvectors1
    return  np.dot(np.dot(A, ev1_r),ev2_r), np.dot(np.dot(eigenVectors1, ev2_r), ev2_r)

def rotate(A, ROT):
    return  np.dot(A, ROT)

def ROT2Euler(ROT):
    ROT00 = float(ROT[0][0])
    ROT10 = float(ROT[0][1])
    ROT20 = float(ROT[0][2])

    ROT01 = float(ROT[1][0])
    ROT11 = float(ROT[1][1])
    ROT21 = float(ROT[1][2])

    ROT02 = float(ROT[2][0])
    ROT12 = float(ROT[2][1])
    ROT22 = float(ROT[2][2])

    Tettax = math.atan2(ROT21, ROT22)
    Tettay = math.atan2(-ROT20, math.sqrt(ROT21*ROT21 + ROT22*ROT22))
    Tettaz = math.atan2(ROT10, ROT00)
    print("Euler angles (yaw, pitch, roll):")
    print("Tetta X = " + str(Tettax * (180 / pi)))
    print("Tetta Y = " + str(Tettay * (180 / pi)))
    print("Tetta Z = " + str(Tettaz * (180 / pi)))
    return Tettax, Tettay, Tettaz

def Euler2ROT(Tettax, Tettay, Tettaz):
    X = [[1, 0, 0], [0, math.cos(Tettax),-math.sin(Tettax)], [0, math.sin(Tettax), math.cos(Tettax)]]
    Y = [[math.cos(Tettay), 0, math.sin(Tettay)], [0 ,1, 0], [-math.sin(Tettay), 0, math.cos(Tettay)]]
    Z = [[math.cos(Tettaz), -math.sin(Tettaz), 0],[math.sin(Tettaz), math.cos(Tettaz), 0], [0, 0, 1]]

    return np.transpose(np.dot(Z,np.dot(Y, X)))
# finds one-to-one correspondence between points
def general_kabsh(P, Q, iteration):
    num_rows_Q, num_cols_Q = Q.shape
    if iteration >= num_rows_Q-1:
        return Q
    best_index = -1
    RMSD_min = rmsd(P, kabsch_rotate(Q, P))
    print("iteration = " + str(iteration) + " row = 0  RMSD = " + str(RMSD_min))
    for i in range(iteration, num_rows_Q):
        R = Q.copy()
        # raws exchange between row with index (iteration) all rows below
        R[iteration, 0], R[i,0] = R[i, 0], R[iteration,0]
        RMSD = rmsd(P, kabsch_rotate(R, P))
        if RMSD < RMSD_min:
            #print("iteration = " + str(iteration) + " row = " + str(i) + " RMSD = " + str(RMSD))
            RMSD_min = RMSD
            best_index = i
    if best_index == -1:
        iteration += 1
        Q = general_kabsh(P, Q, iteration)
        return Q
    else:
        #print ("Best index = " + str(best_index))
        #print ("RMSD old = " + str(rmsd(P, kabsch_rotate(Q, P))))
        Q[[iteration, 0]], Q[[best_index, 0]] = Q[[best_index, 0]], Q[[iteration, 0]]
        iteration +=1
        Q = general_kabsh(P, Q, iteration)
        return Q
        #print ("RMSD new = " + str(rmsd(P, kabsch_rotate(Q, P))))

def rotation_x(points, angle):
    angle = angle*(3.1416/180)
    r = np.array([[1, 0, 0], [0, math.cos(angle), -math.sin(angle)],
                  [0, math.sin(angle), math.cos(angle)]] )
    #p = points.T
    newpoints = np.dot(points,r)
    return np.around(newpoints)

def rotation_y(points, angle):
    angle = angle*(3.1416/180)
    r = np.array([[math.cos(angle), 0, math.sin(angle)], [0, 1, 0],
                  [-math.sin(angle), 0, math.cos(angle)]] )
    #p = points.T
    newpoints = np.dot(points,r)
    return np.around(newpoints)

def rotation_z(points, angle):
    angle = angle*(3.1416/180)
    r = np.array([[math.cos(angle), -math.sin(angle), 0], [math.sin(angle), math.cos(angle), 0],
                  [0, 0, 1]] )
    #p = points.T
    newpoints = np.dot(points,r)
    return np.around(newpoints)



def atom_inside_sphere(atom1, atom2, R):
    distance = math.sqrt((atom1.get_x_coord() - atom2.get_x_coord())**2+(atom1.get_y_coord() - atom2.get_y_coord())**2+
                         (atom1.get_z_coord() - atom2.get_z_coord())**2)
    if distance > R: return False
    if distance <= R: return True

def plot_arr(atoms):
    hits = []
    numbers = []
    count = len(atoms)
    i = 0
    for atom1 in atoms:
        num = 0
        for atom2 in atoms:
            if atom_inside_sphere(atom1, atom2, R):
                num+=1
        hits.append(num)
        numbers.append(atom1.get_atomserialnumber())
        print(str(i) + ' of  ' + str(count))
        i+=1

    #plt.hist(hits)
    plt.figure()
    plt.subplot()
    plt.plot(numbers, hits, 'r-')
    plt.title("Atoms in sphere with R =  " + str(R) + " Angstrems")
    plt.xlabel("Atom number")
    plt.ylabel("Count of atoms")

    #fig = plt.gcf()
    #(username = 'Molodenskiy', api_key = 'Jc4uP32ntpjUyxtMuWhl')
    #py.plot_mpl(fig, filename='2y25')

def subdivide_chains (atoms, chain_identifiers):
    chainatoms = []
    for id in chain_identifiers:
        chainatoms.append([a for a in atoms if a.get_chainidentifier() == id])
    return chainatoms

def subdivide_chains_term(atoms, chain_borders):
    number_of_chains = len(chain_borders) - 1
    allatoms = []
    i = 0
    while i < number_of_chains:
        locatoms = []
        left_border = chain_borders[i]
        right_border = chain_borders[i+1]
        for atom in atoms:
            if atom.get_atomserialnumber() > left_border and atom.get_atomserialnumber() < right_border:
                locatoms.append(atom)
        allatoms.append(locatoms)
        i+=1
    return allatoms

def atoms2pdb(atoms):
    pdb = ""
    for atom in atoms:
        pdb += atom.get_pdb_string()

    return pdb

def pdb2atoms(filename):
    with open(filename, "r") as myfile:
        pdblist = myfile.readlines()
    atoms = []
    for x in pdblist:
        if x.startswith("ATOM"):
            atoms.append(mk_atom(x))
    return atoms
def angle(vec1,vec2):
    return np.arccos(np.clip(np.dot(vec1, np.transpose(vec2)), -1.0, 1.0))

def tilt_twist(atoms1, atoms2):
    # calculation of inertia tensor
    i_tensor1 = inertia_tensor_calc(atoms1)
    i_tensor2 = inertia_tensor_calc(atoms2)
    # calculating eigen values and eigen vectors
    eigenValues1, eigenVectors1 = np.linalg.eig(i_tensor1)
    eigenValues2, eigenVectors2 = np.linalg.eig(i_tensor2)
    # sorting
    idx1 = eigenValues1.argsort()[::-1]
    eigenValues1 = eigenValues1[idx1]
    eigenVectors1 = eigenVectors1[:, idx1]
    idx2 = eigenValues2.argsort()[::-1]
    eigenValues2 = eigenValues2[idx2]
    eigenVectors2 = eigenVectors2[:, idx2]
    #calculate tilt and twist angle
    tilt = math.acos(np.dot(eigenVectors1[0], np.transpose(eigenVectors2[0])))
    twist = math.acos(np.dot(eigenVectors1[2], np.transpose(eigenVectors2[2])))
    return (180 / pi)*tilt, (180 / pi)*twist
def tilt_c(atoms1, atoms2):
    # calculation of inertia tensor
    i_tensor1 = inertia_tensor_calc(atoms1)
    i_tensor2 = inertia_tensor_calc(atoms2)
    # calculating eigen values and eigen vectors
    eigenValues1, eigenVectors1 = np.linalg.eig(i_tensor1)
    eigenValues2, eigenVectors2 = np.linalg.eig(i_tensor2)
    # sorting
    idx1 = eigenValues1.argsort()[::-1]
    eigenValues1 = eigenValues1[idx1]
    eigenVectors1 = eigenVectors1[:, idx1]
    idx2 = eigenValues2.argsort()[::-1]
    eigenValues2 = eigenValues2[idx2]
    eigenVectors2 = eigenVectors2[:, idx2]
    #calculate tilt angle
    a, b = eigenVectors1[0], eigenVectors2[0]
    tilt = angle(a,b)
    #calculate rotational matrix
    v = np.cross(a, b)[0][:]
    c = np.dot(a, np.transpose(b))[0]
    s = np.linalg.norm(v)
    I = np.identity(3)
    vXStr = '{} {} {}; {} {} {}; {} {} {}'.format(0, -v[2], v[1], v[2], 0, -v[0], -v[1], v[0], 0)
    k = np.matrix(vXStr)
    r = I + k + np.multiply(np.dot(k,k), ((1 - c) / (s ** 2)))
    #tilt = math.acos(np.vdot(eigenVectors1[0], eigenVectors2[0]))
    return (180 / pi)*tilt, r
def twist_c(atoms1, atoms2, r):
    # calculation of inertia tensor
    i_tensor1 = inertia_tensor_calc(atoms1)
    i_tensor2 = inertia_tensor_calc(atoms2)
    # calculating eigen values and eigen vectors
    eigenValues1, eigenVectors1 = np.linalg.eig(i_tensor1)
    eigenValues2, eigenVectors2 = np.linalg.eig(i_tensor2)
    # sorting
    idx1 = eigenValues1.argsort()[::-1]
    eigenValues1 = eigenValues1[idx1]
    eigenVectors1 = eigenVectors1[:, idx1]
    idx2 = eigenValues2.argsort()[::-1]
    eigenValues2 = eigenValues2[idx2]
    eigenVectors2 = eigenVectors2[:, idx2]
    #rotate first body to align z
    a = eigenVectors1[2]
    a = r.dot(np.transpose(a))
    #calculate twist angle
    twist = angle(np.transpose(a), eigenVectors2[2])
    #twist = math.acos(np.dot(eigenVectors1[2], np.transpose(eigenVectors2[2])))
    return (180 / pi)*twist
def domains_parsing(chain_identifiers):
    domains_bord = raw_input("For the chain " + str(chain_identifiers[i]) + " (AAAA-BBBB,CCCC-DDDD,... )")
    left_borders = []
    right_borders = []
    doms = domains_bord.split(',')
    for dom in doms:
        borders = dom.split('-')
        if (len(borders) != 2):
            print("Bad domain indication!")
            domains_parsing(chain_identifiers)

        left_borders.append(int(borders[0]))
        right_borders.append(int(borders[1]))
        if (int(borders[0]) > int(borders[1])):
            print ("Left border should be less than right border!")
            domains_parsing(chain_identifiers)

    return left_borders, right_borders
#///////////////////////START PROGRAM////////////////////////////////////////////////////////////

# open and read pdb file
init_pdb_name = raw_input("Please enter a pdb file name ")
with open (init_pdb_name, "r") as myfile:
    pdblist = myfile.readlines()

# create list of atoms
atoms = []
# non-inclusive
chain_borders = []
chain_identifiers = []
chain_borders.append(0)
id = ''
for x in pdblist:
    if x.startswith("ATOM"):
        atoms.append(mk_atom(x))
        if  atoms[-1].get_chainidentifier()!= id: chain_identifiers.append(atoms[-1].get_chainidentifier())
        id = atoms[-1].get_chainidentifier()
    if x.startswith("TER"):
        chain_borders.append(x[7:11])
if len(atoms) < 1:
    print ("Wrong format!")
    exit()

number_of_chains = len(chain_identifiers)
if number_of_chains > 0:
    allchains_atoms = subdivide_chains(atoms, chain_identifiers)
else:
    allchains_atoms = []
    allchains_atoms.append(atoms)
    number_of_chains = 1

print ("The file has " + str(number_of_chains) + " chains:")
for ch in chain_identifiers:
    print (ch)

wth = raw_input("Do you want to use automatic domains subdivision? (y/n) ")
filenames = []

if wth == 'y':
    # subdivision on chains
    i=0
    for at in allchains_atoms:
        pdb_chain = atoms2pdb(at)
        filename = str(init_pdb_name[0:-4]) + "." + chain_identifiers[i]
        file = open(filename, 'w')
        i+=1
        file.write(pdb_chain)
        file.close()
        call(["parcoor", filename, "-p", filename])
        os.remove(filename)
        #save in filenames all found subdomains
        os.chdir("./")
        for file in glob.glob(filename+ "*"):
            filenames.append(file)
            print(file + " is written")
    #Sorting
    #filenames = sorted(filenames, key=lambda name: int(name[-6:-4]))
    #for f in filenames: print(f + " is written")

if wth == 'n':
    print("Please enter the borders between domains manually (in residues) ")
    # subdivision on chains
    i=0
    for at in allchains_atoms:
        print ("Chain " +str(chain_identifiers[i]) + " has residues in range " + str(at[0].get_residuesequencenumber())+
                             "-"+str(at[-1].get_residuesequencenumber()))
        left_borders, right_borders = domains_parsing(chain_identifiers)
        for left_border, right_border in zip(left_borders, right_borders):
            atomsinrange = []
            for atom in at:
                    residue_number = atom.get_residuesequencenumber()
                    if (residue_number >= left_border and residue_number<= right_border): atomsinrange.append(atom)
            if len(atomsinrange) == 0:
                print("No atoms in "+ str(left_border) +"-"+str(right_border) +" region!")
            else:
                pdb_domain = atoms2pdb(atomsinrange)
                filename = "o." + str(init_pdb_name) + ".chain" + str(chain_identifiers[i]) + ".res" \
                       + str(left_border) +"-" +str(right_border)+".pdb"
                filenames.append(filename)
                file = open(filename, 'w')
                file.write(pdb_domain)
                file.close()
                print("File " + str(filename) + " is written")
        i += 1
    #Sorting
    #filenames = sorted(filenames, key=lambda name: int(name[0]))
    #for f in filenames: print(f + " is written")

# now turn of supcomb!
eulersstring = ""
print ("Begining of superposition...")
i = 0
while i < (len(filenames) - 1):
    first_file = filenames[i]
    second_file = filenames[i+1]
    print ("Superposition of " + first_file + " and " + second_file)
    outname = "o."+ second_file + ".aligned2." + first_file
    call(["supcomb", first_file, second_file, "-o", outname])
    # Read and transform matrix into Euler angles
    L = []
    st = ""
    for index, line in enumerate(open(outname)):
        if index <= 19: continue
        elif index ==20: st = line[12:48]
        elif index >= 25: break

        else:
            L.append(line[40:72].split())  # split on whitespace and append value from third columns to list.
    L = np.array(L)
    Tettax,Tettay,Tettaz = ROT2Euler(L)
    #ROT = Euler2ROT(Tettax, Tettay, Tettaz)
    eulersstring += "Superposition of " + first_file + " and " + second_file + ":\n"
    eulersstring += "Tetta z = " + str(Tettaz * (180 / pi)) + "\nTetta y = " + str(Tettay * (180 / pi))+"\nTetta x = "+str(Tettax * (180 / pi)) + "\n"
    # tilt and rotational matrix to align z and z'
    tilt, r = tilt_c(pdb2atoms(second_file), pdb2atoms(outname))
    twist = twist_c(pdb2atoms(second_file), pdb2atoms(outname), r)
    eulersstring += "Tilt = " +str(tilt) +  ", Twist = " + str(twist) + "\n"
    print("************************************************")
    print("Twist = " + str(float(twist)))
    print("Tilt = " + str(float((tilt))))
    eulersstring += st + "\n\n"
    #print("from supcomb")
    #print (L)
    #print("after double transformation")
    #print(ROT)
    #####################################################################################################
    i+=1
file = open("Results.txt", 'w')
file.write(eulersstring)
file.close()
end = time.time()
print (" Time consumed : " + str(end - start) + " sec")


#######################################################################################3