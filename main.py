from pymol import cmd, stored
from pymol.cgo import *
from colour import Color
import numpy as np
import cclib 
import os

def color_gradient(vals):
    c1 = Color('blue') #Blue
    c2 = Color('red') #Red
    min_vals,max_vals = min(vals), max(vals)
    ix = [int((x-min_vals)/(max_vals-min_vals)*100) for x in vals]
    crange = list(c1.range_to(c2,101))
    mycolor = []
    for x in ix:
        myc = crange[x].get_rgb()
        mycolor.append([COLOR, myc[0], myc[1],myc[2]])
    return mycolor
    
def load_dots(Coordinates, name="ply", dotSize=0.1):
    Coordinates = np.array(Coordinates)
    Coordinates = Coordinates[Coordinates[:,3].argsort()]
    verts = [[float(x[0]), float(x[1]), float(x[2])] for x in Coordinates]
    obj = []

    colors = color_gradient([float(x[3]) for x in Coordinates])

    for v_ix in range(len(verts)):
        vert = verts[v_ix]
        obj.extend(colors[v_ix])
        obj.extend([SPHERE, vert[0], vert[1], vert[2], dotSize])
    cmd.load_cgo(obj, name, 1.0)

def load_structure(filename):
    parser = cclib.io.ccopen(filename)
    data = parser.parse()
    data.writexyz('Coordinates_from_gaussian_file.xyz')
    cmd.load('Coordinates_from_gaussian_file.xyz')
    os.remove('Coordinates_from_gaussian_file.xyz')

def load_gaussian_esp(Gaussian_file):
    f = open(Gaussian_file, 'r'); f = f.readlines(); f = [i.split('\n')[0] for i in f ]
    ESP_coordinates = []
    
    for line in f:
        if 'ESP Fit Center' not in line:
            continue
        line = [i for i in line.split(' ') if len(i)>0]
        ESP_coordinates.append(line[6:9])
    
    counter = 0
    for line in f:
        if (' Fit ' in line) and ('ESP' not in line):            
            line = [i for i in line.split(' ') if len(i)>0]
            ESP_coordinates[counter].append(line[2])
            counter = counter + 1
    cmd.set('connect_cutoff',0.3)
    load_structure(Gaussian_file)
    load_dots(ESP_coordinates, name='QM ESP')

def Read_RESP_charges_file(RESP_charges_file):

    code = {"1" : "H", "2" : "He", "3" : "Li", "4" : "Be", "5" : "B", \
            "6" : "C", "7" : "N", "8" : "O", "9" : "F", "10" : "Ne", \
            "11" : "Na" , "12" : "Mg" , "13" : "Al" , "14" : "Si" , "15" : "P", \
            "16" : "S" , "17" : "Cl" , "18" : "Ar" , "19" : "K" , "20" : "Ca", \
            "21" : "Sc" , "22" : "Ti" , "23" : "V" , "24" : "Cr" , "25" : "Mn", \
            "26" : "Fe" , "27" : "Co" , "28" : "Ni" , "29" : "Cu" , "30" : "Zn", \
            "31" : "Ga" , "32" : "Ge" , "33" : "As" , "34" : "Se" , "35" : "Br", \
            "36" : "Kr" , "37" : "Rb" , "38" : "Sr" , "39" : "Y" , "40" : "Zr", \
            "41" : "Nb" , "42" : "Mo" , "43" : "Tc" , "44" : "Ru" , "45" : "Rh", \
            "46" : "Pd" , "47" : "Ag" , "48" : "Cd" , "49" : "In" , "50" : "Sn", \
            "51" : "Sb" , "52" : "Te" , "53" : "I" , "54" : "Xe" , "55" : "Cs", \
            "56" : "Ba" , "57" : "La" , "58" : "Ce" , "59" : "Pr" , "60" : "Nd", \
            "61" : "Pm" , "62" : "Sm" , "63" : "Eu" , "64" : "Gd" , "65" : "Tb", \
            "66" : "Dy" , "67" : "Ho" , "68" : "Er" , "69" : "Tm" , "70" : "Yb", \
            "71" : "Lu" , "72" : "Hf" , "73" : "Ta" , "74" : "W" , "75" : "Re", \
            "76" : "Os" , "77" : "Ir" , "78" : "Pt" , "79" : "Au" , "80" : "Hg", \
            "81" : "Tl" , "82" : "Pb" , "83" : "Bi" , "84" : "Po" , "85" : "At", \
            "86" : "Rn" , "87" : "Fr" , "88" : "Ra" , "89" : "Ac" , "90" : "Th", \
            "91" : "Pa" , "92" : "U" , "93" : "Np" , "94" : "Pu" , "95" : "Am", \
            "96" : "Cm" , "97" : "Bk" , "98" : "Cf" , "99" : "Es" ,"100" : "Fm", \
            "101": "Md" ,"102" : "No" ,"103" : "Lr" ,"104" : "Rf" ,"105" : "Db", \
            "106": "Sg" ,"107" : "Bh" ,"108" : "Hs" ,"109" : "Mt" ,"110" : "Ds", \
            "111": "Rg" ,"112" : "Uub","113" : "Uut","114" : "Uuq","115" : "Uup", \
            "116": "Uuh","117" : "Uus","118" : "Uuo"}

    f = open(RESP_charges_file,'r'); f = f.readlines(); f = [line.strip('\n').split() for line in f]

    molecule_mass = []

    Molecule_number = False
    Molecule_counter = -1
    Atom_Coordinates = False
    Atom_name = False

    for line in f:
        if 'TITLE' in line: Molecule_number = True; continue
        if ['atm.no','X','Y','Z'] == line: Atom_Coordinates = True; continue
        if "q(opt)" in line: Atom_name = True; continue

        if Molecule_number:
            if len(line) == 0: Molecule_number = False; continue
            Molecule_counter += 1
            molecule_mass.append([Molecule_counter,[]])
            Atom_counter = 0
        
        if Atom_Coordinates: 
            if len(line) == 0: Atom_Coordinates = False;continue
            molecule_mass[Molecule_counter][1].append([[],line[1],line[2],line[3]])
            continue
        
        if Atom_name:
            if len(line) == 0: Atom_name = False; continue            
            molecule_mass[Molecule_counter][1][Atom_counter][0] = code[line[1]]
            Atom_counter += 1
            continue
    return molecule_mass

def Read_RESP_outfile(RESP_outfile):
    f = open(RESP_outfile,'r'); f = f.readlines(); f = [line.strip('\n').split() for line in f]
    
    ESP_mass = [[],[],[]] #QM/PC/Diff
    ESP_counter = -1
    ESP_bool = False

    for line in f:
        if 'esp_qm' in line: ESP_bool = True; ESP_counter +=1; ESP_mass[0].append([]);ESP_mass[1].append([]);ESP_mass[2].append([]); continue
        if 'molecule' in line: ESP_bool = False; continue
        if ESP_bool: 
            ESP_mass[0][ESP_counter].append([line[0],line[1],line[2],line[3]])
            ESP_mass[1][ESP_counter].append([line[0],line[1],line[2],line[4]])
            ESP_mass[2][ESP_counter].append([line[0],line[1],line[2],line[5]])
    
    return ESP_mass

def load_RESP_charges(files):
    RESP_charges_file, RESP_outfile, molecule = files.split()
    XYZ = Read_RESP_charges_file(RESP_charges_file)
    ESP_mass = Read_RESP_outfile(RESP_outfile)
        
    if len(molecule.split(':')) == 2:        
        for mol in range(int(molecule.split(':')[0]), int(molecule.split(':')[1])+1):  

            if len(ESP_mass[0]) < mol-1:
                break
            
            load_dots(ESP_mass[0][int(mol)], name='Molecule '+str(mol)+' QM ESP')
            load_dots(ESP_mass[1][int(mol)], name='Molecule '+str(mol)+' PC ESP')
            load_dots(ESP_mass[2][int(mol)], name='Molecule '+str(mol)+' Diff ESP')
            
            f = open('Molecule_'+str(mol)+'.xyz','a')
            f.write(str(len(XYZ[int(mol)][1]))+'\n')
            f.write('Coordinates_'+str(mol)+'\n')
            for atom in range(len(XYZ[int(mol)][1])):
                f.write("{:4} {:11.6f} {:11.6f} {:11.6f}\n".format(XYZ[int(mol)][1][atom][0],float(XYZ[int(mol)][1][atom][1]),float(XYZ[int(mol)][1][atom][2]),float(XYZ[int(mol)][1][atom][3])))
            f.write('\n')
            f.close()
            cmd.set('connect_cutoff',2)
            cmd.load('Molecule_'+str(mol)+'.xyz',format='xyz')
            os.remove('Molecule_'+str(mol)+'.xyz')            

    else:   

        load_dots(ESP_mass[0][int(molecule)], name='Molecule '+str(molecule)+' QM ESP')
        load_dots(ESP_mass[1][int(molecule)], name='Molecule '+str(molecule)+' PC ESP')
        load_dots(ESP_mass[2][int(molecule)], name='Molecule '+str(molecule)+' Diff ESP')

        f = open('Molecule_'+str(molecule)+'.xyz','a')
        f.write(str(len(XYZ[int(molecule)][1]))+'\n')
        f.write('Coordinates_'+str(molecule)+'\n')

        for atom in range(len(XYZ[int(molecule)][1])):
                f.write("{:4} {:11.6f} {:11.6f} {:11.6f}\n".format(XYZ[int(molecule)][1][atom][0],float(XYZ[int(molecule)][1][atom][1]),float(XYZ[int(molecule)][1][atom][2]),float(XYZ[int(molecule)][1][atom][3])))
        f.write('\n')
        f.close()
        cmd.set('connect_cutoff',2)
        cmd.load('Molecule_'+str(molecule)+'.xyz','Molecule_'+str(molecule)+'.xyz')
        os.remove('Molecule_'+str(molecule)+'.xyz') 



cmd.extend("load_gaussian_esp", load_gaussian_esp)
cmd.extend("load_RESP_files", load_RESP_charges)