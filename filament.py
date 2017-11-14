from subprocess import call
pi = 3.14159265359
call(["ls", "-l"])
import sys
import glob
import math
import os
import numpy as np
import time
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
    weights = {'H' : 1.008,'C' : 12.011,'N' : 14.007,'O' : 15.999,'S' : 32.060, 'P' : 30.974}
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
        #print("Atom: " + atomname.replace(" ","")[0])
        self._weight = self.weights[atomname.replace(" ","")[0]]


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

#calculation of inertia tensor from set of atoms
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

#calculation of transformational matrix from Euler angles (x-y-z)
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

#calculation of Euler angles (x-y-z) from transformational matrix
def Euler2ROT(Tettax, Tettay, Tettaz):
    X = [[1, 0, 0], [0, math.cos(Tettax),-math.sin(Tettax)], [0, math.sin(Tettax), math.cos(Tettax)]]
    Y = [[math.cos(Tettay), 0, math.sin(Tettay)], [0 ,1, 0], [-math.sin(Tettay), 0, math.cos(Tettay)]]
    Z = [[math.cos(Tettaz), -math.sin(Tettaz), 0],[math.sin(Tettaz), math.cos(Tettaz), 0], [0, 0, 1]]

    return np.transpose(np.dot(Z,np.dot(Y, X)))

#creates n lists of atoms, where n - number of chains
def subdivide_chains (atoms, chain_identifiers):
    chainatoms = []
    for id in chain_identifiers:
        chainatoms.append([a for a in atoms if a.get_chainidentifier() == id])
    return chainatoms

#reverse transformation from list of atoms to pdb file
def atoms2pdb(atoms):
    pdb = ""
    for atom in atoms:
        pdb += atom.get_pdb_string()

    return pdb

#transforms pdb file to list of atoms
def pdb2atoms(filename):
    with open(filename, "r") as myfile:
        pdblist = myfile.readlines()
    atoms = []
    for x in pdblist:
        if x.startswith("ATOM"):
            atoms.append(mk_atom(x))
    return atoms

#calculates angle between two vectors
def angle(vec1,vec2):
    return np.arccos(np.clip(np.dot(vec1, np.transpose(vec2)), -1.0, 1.0))

#calculation of tilt angle between Z1 and Z2 and of transformational matrix of this operation
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
    return (180 / pi)*tilt, r

#calculation of twist angle between X1 and X2', where X2' is X2 after alignment of axes Z1 and Z2
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
    return (180 / pi)*twist

#dialog with user and borders parsing
def domains_parsing(chain_identifiers):
    suc = True
    left_borders = []
    right_borders = []
    while suc:
        suc = False
        domains_bord = raw_input("For the chain " + str(chain_identifiers[i]) + " (AAAA-BBBB,CCCC-DDDD,... )")
        domains_bord = domains_bord.replace(" ", "")
        doms = domains_bord.split(',')
        for dom in doms:
            borders = dom.split('-')
            if (len(borders) != 2):
                print("Bad domain indication!")
                suc = True
                continue
            try:
                left = int(borders[0])
                right = int(borders[1])
            except ValueError:
                print("Please, use numbers only!")
                suc = True
                continue
            if (left > right):
                print ("Left border should be less or equal to right border!")
                suc = True
                continue
            left_borders.append(int(borders[0]))
            right_borders.append(int(borders[1]))
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
print ("Beginning of superposition...")
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
    ###################################################################################
    #ROT = Euler2ROT(Tettax, Tettay, Tettaz)
    #print("from supcomb")
    #print (L)
    #print("after double transformation")
    #print(ROT)
    ###################################################################################
    i+=1
file = open("Results.txt", 'w')
file.write(eulersstring)
file.close()
end = time.time()
print (" Time consumed : " + str(end - start) + " sec")

##############################END PROGRAM##############################################
