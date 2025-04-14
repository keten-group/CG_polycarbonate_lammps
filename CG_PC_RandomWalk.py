"""
*************************************************************
*** Please cite ref [1] this script is used in your work. ***
*************************************************************

Author: Heather White, Department of Mechanical Engineering, Northwestern University
Acknowledgements: Zhenghao Wu @ Xi'an Jiaotong-Liverpool University, Zhaoxu Meng @
  Clemson University, Wenjie Xia @ Iowa State University
Title: CG_PC_RandomWalk_ToPublish.py
Programming language: Python 3.9.6
Date: Last updated April 13, 2025

Description
* Generates a LAMMPS input file containing coordinates and bonds for the 
  coarse-grained (CG) model of polycarbonate (PC)
* CG PC bead placement is based on a self-avoiding random walk with minimum angle 
  (~84 degrees)
* Intended for use with the CG PC model developed in refs [1,2]
* Published alongside ref [1]

Notes
* The generated PC melt is has periodic boundary conditions in x and y, but not z. 
  The z boundary can be made perdiodic through use of LAMMPS boundary settings 
  during relaxation.
* PC CG bead types are 3, 4, 5, and 6. Types 1 and 2 are an artifact of development. 
  This was left in place for consistency with LAMMPS commands provided in other
  documents.

References
[1] White, H., Fermen-Coker, M., Chen, W., Keten, S. (2025) "Characterizing the 
Mechanical Response of a Polycarbonate Coarse-Grained Model Developed with Energy 
Renormalization". Macromolecules. (DOI to be established on publication.)
[2] Xia, W., Hansoge, N., Xu, W., Phelan, F. R., Keten, S., Douglas, J. (2019) "Energy 
renormalization for coarse-grained polymers having different segmental structures". 
Science Advances, 5(4), eaav4683. DOI:10.1126/sciadv.aav4683

********** BEGIN COPYRIGHT INFORMATION **********

Copyright 2025 Heather L. White

Permission is hereby granted, free of charge, to any person obtaining a copy of this 
software and associated documentation files (the “Software”), to deal in the Software 
without restriction, including without limitation the rights to use, copy, modify, merge,
publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons 
to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or 
substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
DEALINGS IN THE SOFTWARE.

********** END COPYRIGHT INFORMATION **********
"""

### Imports

import math
import random
import numpy as np

### User inputs (basic)
print_info = True                     # Toggle for printing change placement info to terminal
side_length = 100                     # Desired side length (x=y, Angstroms) of the polycarbonate block
height = 100                          # Polymer thickness (z, Angstroms)
clength = 20                          # Chain length (# monomers) for PC chains
volume_multiplier = 6.0               # Adjust to change density where higher values result
                                      #   in lower density and easier bead placement during 
                                      #   the random walk. A value of 6.0 is recommended.

### Parameters for PC

massA = 76                                # Bead A - phenylene (g/mol)
massB = 42                                # Bead B - isopropylidene (g/mol)
massC = 76                                # Bead C - phenylene (g/mol)
massD = 60                                # Bead D - carbonate (g/mol)
monomer_mass = (massA + massB + massC     # Total mass of a monomer (g/mol)
                + massD)
rho = 1.3                                 # Desired polymer density, g/cm^3. This value is a
                                          #   little higher than typical (1.2) to account for
                                          #   the zspacing

alength = clength*4                       # Number of beads per polymer chain

### calculate the volume of the polymer phase and corresponding number of chains

Ntotal_ap = height*side_length*side_length/(monomer_mass/       # Number of monomers in each chunk of polymer needed
                    6.02E23/rho*1E24)              #   to reach target density
nchains = math.ceil(Ntotal_ap/clength)             # Total number of chains in each chunk of polymer
Ntotal = nchains*clength                           # Final number of monomers in each chunk of polymer

b_len_pc = 3.16                            # Bond length PC
min_angle_distance = 4.44 #3.75                               # Minimum distance between atom 1 and atom 3
                                          #   (provides minimum angle, 3.75 = ~84 deg)
                                          #   (4.44 = ~100 deg)

### box boundaries

# Increase volume for random walk 9.29.23
lx = side_length*(volume_multiplier**1/3)
ly = side_length*(volume_multiplier**1/3)
lz = height*(volume_multiplier**1/3)

zspacing = 5                              # Extra wiggle room
comp_spacing = lz + zspacing              # Height of polymer + wiggle room

xlo = -lx/2.0; xhi = lx/2.0
ylo = -ly/2.0; yhi = ly/2.0
zlo = -lz/2.0; zhi = lz/2.0

### Random walk PC_4bead model ###

# Determine the box dimensions with the assembly size and the extra box space
minx = round(xlo, 3); maxx = round(xhi, 3)
miny = round(ylo, 3); maxy = round(yhi, 3)
x_box_length = 2*maxx; y_box_length = 2*maxy

# Generate coordinates for random walk
# The while loop continuously generates new random coordinates
# Three checks to pass: z boundary, angle, and overlap
random.seed()
coor_unwrap_all_chain = np.empty([0,3])      # Storage for coordinates for all chains
coord_wrap_all_chain = np.empty([0,3])      # Storage for wrapped coordinates for all chains

# Loop through each chain
for jj in range(0,nchains):                 
  
  if print_info == True:
    print("Creating " + str(jj+1) + " of " + str(nchains) + " chains...")
  
  # Reset boundary checks and attempt count
  coor_unwrap_this_chain = np.empty([0,3])          # Storage for coordinates within a single chain
  coord_wrap_this_chain = np.empty([0,3])                  # Storage for wrapped coordinates 
  x_boundary_check = 0
  y_boundary_check = 0
  attempts = 0
  
  # Loop through each monomer within each chain
  beadnum = -1
  while beadnum < alength-1:
    beadnum += 1

    # Pick new coordinates until loop is broken
    while True:  
      
      # Manage chain restarts
      if beadnum > 0 and attempts > 0 and (attempts%1000) == 0:
        print("\t", attempts, "attempts at placing bead", beadnum, "of", alength)
      if attempts > 2000:
        print("\tToo many attempts. Retrying...")
        x_boundary_check = 0; y_boundary_check = 0
        attempts = 0
        beadnum = -1
        coor_unwrap_this_chain = np.empty([0,3])          # Storage for coordinates within a single chain
        coord_wrap_this_chain = np.empty([0,3])           # Storage for wrapped coordinates 
        break
        
      # If this is the first bead, generate its position randomly in the box    
      if beadnum == 0: 
                                    
        xij = minx + random.random()*x_box_length          
        yij = miny + random.random()*y_box_length
        zij = zlo + random.random()*lz

        xij_wrapped = xij
        yij_wrapped = yij
        zij_wrapped = zij

      # If not first bead, generate coordinates relative to the previous bead
      else:                   
                                               
        # Get a random direction in space (-1,1)                       
        dx = 2.0*random.random() - 1.0      
        dy = 2.0*random.random() - 1.0
        dz = 2.0*random.random() - 1.0
        r1 = math.sqrt(dx**2 + dy**2 + dz**2)
        
        # Normalize direction to unit vectors
        dx = dx/r1                          
        dy = dy/r1
        dz = dz/r1

        # Tentative position of the next bead in each chain relative to previous bead
        xij = coor_unwrap_this_chain[beadnum-1,0] + dx*b_len_pc    
        yij = coor_unwrap_this_chain[beadnum-1,1] + dy*b_len_pc
        zij = coor_unwrap_this_chain[beadnum-1,2] + dz*b_len_pc

        # Determine wrapped coordinates (z is nonperiodic)
        x_boundary_check = math.floor((xij + x_box_length/2) / x_box_length)
        xij_wrapped = xij - x_boundary_check*x_box_length
        y_boundary_check = math.floor((yij + y_box_length/2) / y_box_length)
        yij_wrapped = yij - y_boundary_check*y_box_length
        
      new_bead_coord_unwrap = np.array([[xij, yij, zij]])  
      new_bead_coord_wrap = np.array([[xij_wrapped, yij_wrapped, zij]])

      ### CHECK 1: Z BOUNDARY ###
      # If the z (nonperiodic) box boundaries are breached, get new tentative coordinates
      if (zij<zlo) or (zij>zhi):
        attempts += 1
        continue             
      
      ### CHECK 2: ANGLE ###
      # If the bead is at least number 3 in the chain, check the angle it is creating
      if beadnum >= 2: 
        r_squared = sum([(new_bead_coord_unwrap[0,dim] - coor_unwrap_this_chain[beadnum-2,dim])**2 for dim in range(3)])
        # If the angle between beads in a chain is big enough, angle check passes
        if math.sqrt(r_squared) < min_angle_distance:  
          attempts += 1
          continue
                
      ### CHECK 3: OVERLAP ###
      # Make sure wrapped coordinates don't overlap beads
      
      continue_while_loop_flag = False
      # Loop through wrapped coordinates of other chains
      for index in range(0,len(coord_wrap_all_chain)):
        
        wrapped_other_chain_bead_coords_to_check = [(new_bead_coord_wrap[0,dim] - coord_wrap_all_chain[index, dim])**2 for dim in range(3)]
        r_sq_other_chains = math.sqrt(sum(wrapped_other_chain_bead_coords_to_check))

        if r_sq_other_chains < 6.13: 
          continue_while_loop_flag = True
          break # break for loop

      # Loop through wrapped coordinates of this chain, excluding the previously placed two beads
      for index in range(0,len(coord_wrap_this_chain)-3):
        
        wrapped_this_chain_bead_coords_to_check = [(new_bead_coord_wrap[0,dim] - coord_wrap_this_chain[index, dim])**2 for dim in range(3)]
        distr_this_chain = math.sqrt(sum(wrapped_this_chain_bead_coords_to_check))

        if distr_this_chain < b_len_pc: 
          continue_while_loop_flag = True
          break # break for loop
  
      if continue_while_loop_flag == True:
          attempts += 1
          continue
    
      # Reset attempts and accept new placement
      attempts = 0
      break          
    # End while loop

    # Add bead coordinates to this chain's beads
    if beadnum != -1: # If not restarting chain
      coor_unwrap_this_chain = np.append(coor_unwrap_this_chain, new_bead_coord_unwrap, axis=0)
      coord_wrap_this_chain = np.append(coord_wrap_this_chain, new_bead_coord_wrap, axis=0)

  # If not restarting chain
  if beadnum != -1:
    # Add the coordinates for the chain to the total coordinate list for the PC chunk
    coor_unwrap_all_chain=np.append(coor_unwrap_all_chain, coor_unwrap_this_chain, axis=0)
    coord_wrap_all_chain=np.append(coord_wrap_all_chain, coord_wrap_this_chain, axis=0)
    beadnum = -1

###Generate the coordinates for one layer of polymer

coords_p = np.zeros([0,3])                  # New matrix for first chunk of polymer coordinates

 # Add polymer coordinates from random walk to new matrix
for i in range(0,len(coor_unwrap_all_chain)):                 
  coords_p = np.append(coords_p, [coor_unwrap_all_chain[i,:]], axis=0)

### Coordinates for all polymer chunks

# New matrix for first chunk of polymer coordinates
coords = np.zeros([0,3])     

# Loop through each layer                     
for i in range(0,1):                          
  
  # x-coordinates
  coords1 = np.array([coords_p[:,0]]).T 
  # y-coordinates          
  coords2 = np.array([coords_p[:,1]]).T            
  # z-coordinates                         
  coords3 = np.array([coords_p[:,2]]).T
  # Concatentate x, y, and z coordinates
  coords123 = np.concatenate((coords1, coords2, coords3), axis = 1)
  # Append to list of coordinates for all polymer chunks
  coords = np.concatenate((coords, coords123), axis=0)                 

### Atom index

monomers_chain = clength                          # Number of monomers per polymer chain
chains = nchains                                  # Number of chains per chunk
atoms_p = chains*monomers_chain*4                 # Total number of polymer beads (all chunks)

# Total number of beads 
num = np.arange(1,atoms_p+1)                     # Array for numbering all the beads (g+p)

### Mol index

molecule = np.empty([1,chains*monomers_chain*4])  # Empty array for numbering the molecules
for i in range(1,chains+1):                  # Loop through all the polymer chains
  for j in range(1,monomers_chain*4+1):           # Loop through all the monomers in each chain
    molecule[0,j+(i-1)*monomers_chain             # Assign all beads in a chain a single
              *4-1] = i

### Charge (all 0)

charge = np.empty([1,atoms_p])                    # Empty array for polymer bead charges
for i in range(1,atoms_p+1):                      # Set all polymer bead charges to zero
  charge[0, i-1] = 0

### Atom type

type = np.empty([1,atoms_p])                      # Empty array for polymer bead types

for i in range(0,atoms_p,4):                      # Alternate bead types between 3, 4, 5, 6
  type[0,i] = 3
for i in range(1,atoms_p,4):
  type[0,i] = 4
for i in range(2,atoms_p,4):
  type[0,i] = 5
for i in range(3,atoms_p,4):
  type[0,i] = 6

### Combine polymer beads into a single coordinate matrix

num = np.array([num]).T
molecule = molecule.T
type = type.T
charge = charge.T

coordinates = np.empty((0,7))

# Add all bead information (num, mol, type, charge, coord) to this matrix
data_p = np.concatenate((num, molecule, type, charge, coords), axis=1) 
# Add the polymer bead info matrix to the total bead info matrix that contains g
coordinates = np.concatenate((coordinates, data_p), axis=0)   

### Polymer connections

bonds_chain = alength-1                           # Calculate number of bonds, angles, and
angles_chain = alength-2                          #   dihedrals in a single polymer chain
dihedrals_chain = alength-3

bonds_p = bonds_chain*chains                       # Total number of polymer bonds, angles, and
angles_p = angles_chain*chains                     #   dihedrals in the assembly
dihedrals_p = dihedrals_chain*chains

bondtypes = 4                                     # Numbers of polymer bond types, angle types,
angletypes = 4                                    #   and dihedral types
dihedraltypes = 4

### Polymer bonds

### Bonds matrix

bond = np.zeros([bonds_p,4])                      # Create an empty matrix for the polymer bond
                                                  #   information
bond[:,0] = np.arange(1,bonds_p+1)                # Number the poly bond IDs (first column)

for k in range(1,chains+1):                       # Loop through the each polymer chain in the
                                                  #   assembly
  bond[range((k-1)*bonds_chain+1-1,               # Bond type 2, A-B
             k*bonds_chain,bondtypes),1] = 2
  bond[range((k-1)*bonds_chain+2-1,               # Bond type 3, B-C
             k*bonds_chain, bondtypes),1] = 3
  bond[range((k-1)*bonds_chain+3-1,               # Bond type 4, C-D
             k*bonds_chain, bondtypes),1] = 4
  bond[range((k-1)*bonds_chain+4-1,               # Bond type 5, D-E
             k*bonds_chain, bondtypes),1] = 5

  bond[range((k-1)*bonds_chain, k*bonds_chain)    # Specify the first bead in each polymer bond
      ,2]=range((k-1)*monomers_chain*4+1, k*monomers_chain*4)
  bond[range((k-1)*bonds_chain, k*bonds_chain)    # Specify the second bead in each poly bond
      ,3]=range((k-1)*monomers_chain*4+1+1, k*monomers_chain*4+1)

bonds = bond.copy()

rownum = 0
for row in bond:
  rownum += 1
  if rownum == 2*(bonds_chain):
    break

###Polymer angles

### Angles matrix

angle = np.zeros([angles_p,5])                    # Create an empty matrix for the poly angle
                                                  #   information
angle[:,0] = np.arange(1,angles_p+1)              # Number the poly angle IDs (first column)

for i in range(1,chains+1):                       # Loop through each polymer chain in the
                                                  #   assembly
  angle[range((i-1)*angles_chain+1-1,             # Angle type 2, A-B-C
              i*angles_chain, angletypes),1] = 2
  angle[range((i-1)*angles_chain+2-1,             # Angle type 3, B-C-D
              i*angles_chain, angletypes),1] = 3
  angle[range((i-1)*angles_chain+3-1,             # Angle type 4, C-D-A
              i*angles_chain, angletypes),1] = 4
  angle[range((i-1)*angles_chain+4-1,             # Angle type 5, D-A-B
              i*angles_chain, angletypes),1] = 5

  angle[range((i-1)*angles_chain+1-1,             # First bead in the polymer angle
              i*angles_chain),2]=range((i-1)
              *monomers_chain*4+1,
              i*monomers_chain*4-1)
  angle[range((i-1)*angles_chain+1-1,             # Second bead in the polymer angle
              i*angles_chain),3]=range((i-1)
              *monomers_chain*4+1+1,
              i*monomers_chain*4)
  angle[range((i-1)*angles_chain+1-1,             # Third bead in the polymer angle
              i*angles_chain),4]=range((i-1)
              *monomers_chain*4+1+2,
              i*monomers_chain*4+1)

angles = angle.copy()

rownum = 0
for row in angle:
  rownum += 1
  if rownum == 2*(angles_chain):
    break

###Polymer dihedrals

### Dihedrals matrix
dihedral = np.zeros([dihedrals_p,6])              # Create an empty matrix for the polymer
                                                  #   dihedral information
dihedral[:,0] = (np.arange(1,dihedrals_p+1))       # Number the polyner dihedral IDs (first
                                                  #   column)

for i in range(1,chains+1):                       # Loop through each polymer chain in the
                                                  #   assembly
  dihedral[range((i-1)*dihedrals_chain+1-1,       # Dihedral type 2, A-B-C-D
                 i*dihedrals_chain,
                 dihedraltypes),1] = 2
  dihedral[range((i-1)*dihedrals_chain+2-1,       # Dihedral type 3, B-C-D-A
                 i*dihedrals_chain,
                 dihedraltypes),1] = 3
  dihedral[range((i-1)*dihedrals_chain+3-1,       # Dihedral type 4, C-D-A-B
                 i*dihedrals_chain,
                 dihedraltypes),1] = 4
  dihedral[range((i-1)*dihedrals_chain+4-1,       # Dihedral type 5, D-A-B-C
                 i*dihedrals_chain,
                 dihedraltypes),1] = 5
  dihedral[range((i-1)*dihedrals_chain+1-1,       # First bead in the polymer dihedral
                 i*dihedrals_chain),2]=range((i-1)
                 *monomers_chain*4+1,
                 i*monomers_chain*4-2)
  dihedral[range((i-1)*dihedrals_chain+1-1,       # Second bead in the polymer dihedral
                 i*dihedrals_chain),3]=range((i-1)
                 *monomers_chain*4+1+1,
                 i*monomers_chain*4-1)
  dihedral[range((i-1)*dihedrals_chain+1-1,       # Third bead in the polymer dihedral
                 i*dihedrals_chain),4]=range((i-1)
                 *monomers_chain*4+1+2,
                 i*monomers_chain*4)
  dihedral[range((i-1)*dihedrals_chain+1-1,       # Fouth bead in the polymer dihededral
                 i*dihedrals_chain),5]=range((i-1)
                 *monomers_chain*4+1+3,
                 i*monomers_chain*4+1)

dihedrals = dihedral.copy()

rownum = 0
for row in dihedral:
  rownum += 1
  if rownum == 2*(dihedrals_chain):
    break

###Write data file
# Create a custom name for the assembly LAMMPS structure file
filename = ('pc_h{height}nm_s{side}nm_{cl}mon.data').format(height=int(round(height/10,1)), side = int(round(side_length/10,1)), cl=clength)

outfile = open(filename, "w")                                     # Open the file for writing
outfile.write("Coarse-grained polycarbonate for LAMMPS\n\n")      # File header

natoms = len(coordinates)                         # Determine the number of coordinates, bonds,
nbonds = len(bonds)                               #   angles, dihedrals, and impropers
nangles = len(angles)
ndihedrals = len(dihedrals)
nimpropers = 0

outfile.write(str(natoms) + " atoms\n")           # Write the above values to the file
outfile.write(str(nbonds) + " bonds\n")
outfile.write(str(nangles) + " angles\n")
outfile.write(str(ndihedrals) + " dihedrals\n")
outfile.write(str(nimpropers) + " impropers\n\n")

natomtypes = 6                                    # Record the number of atom, bond, angle,
nbondtypes = 5                                    #   and dihedral types
nangletypes = 5
ndihedraltypes = 5

outfile.write(str(natomtypes) + " atom types\n")  # Write the above values to the file
outfile.write(str(nbondtypes) + " bond types\n")
outfile.write(str(nangletypes) + " angle types\n")
outfile.write(str(ndihedraltypes) + " dihedral "
              +"types\n\n")

roomz = 10                                       # Extra box space in x, y, and z directions

# Determine the box dimensions with the assembly size and the extra box space
minx = round(xlo, 3); maxx = round(xhi, 3)
miny = round(ylo, 3); maxy = round(yhi, 3)
minz = round(min(coordinates[:,6])-0.5*roomz, 3)
maxz = round(max(coordinates[:,6])+0.5*roomz, 3)

# Write the box information to the file
outfile.write(str(minx) + " " + str(maxx) + " xlo xhi\n")
outfile.write(str(miny) + " " + str(maxy) + " ylo yhi\n")
outfile.write(str(minz) + " " + str(maxz) + " zlo zhi\n\n")

outfile.write("Masses\n\n")                       # Arrange the bead mass info
mass_G = 1.00
mass_G2 = 1.00                                  
masses = np.array([[1, mass_G], [2, mass_G2], [3, massA], [4, massB], [5, massC], [6, massD]])

for i in range(1, len(masses)+1):                 # Write the bead mass info to the file
  prov = masses[i-1,:]
  outfile.write(str(int(prov[0])) + " " + str(prov[1]) + "\n")

### print atoms

outfile.write("\nAtoms\n\n")

for i in range(1,len(coordinates)+1):             # Arrange the coordinate info
  prov = ('{x0} {x1} {x2} {x3} {x4} {x5} {x6}\n').format( \
      x0 = str(int(coordinates[i-1,0])), \
      x1 = str(int(coordinates[i-1,1])), \
      x2 = str(int(coordinates[i-1,2])), \
      x3 = str(round(coordinates[i-1,3],3)), \
      x4 = str(round(coordinates[i-1,4],3)), \
      x5 = str(round(coordinates[i-1,5],3)), \
      x6 = str(round(coordinates[i-1,6],3)))

  outfile.write(str(prov))                        # Write the coordinate info to the file

### print bonds

outfile.write("\nBonds\n\n")

for i in range(1,len(bonds)+1):                   # Arrange the bond info
  prov = ('{x0} {x1} {x2} {x3}\n').format( \
      x0 = str(int(bonds[i-1,0])), \
      x1 = str(int(bonds[i-1,1])), \
      x2 = str(int(bonds[i-1,2])), \
      x3 = str(int(bonds[i-1,3])))

  outfile.write(str(prov))                        # Write the bond info to the file

### print angles

outfile.write("\nAngles\n\n")

for i in range(1,len(angles)+1):                  # Arrange the angle info
  prov = ('{x0} {x1} {x2} {x3} {x4}\n').format( \
      x0 = str(int(angles[i-1,0])), \
      x1 = str(int(angles[i-1,1])), \
      x2 = str(int(angles[i-1,2])), \
      x3 = str(int(angles[i-1,3])), \
      x4 = str(int(angles[i-1,4])))

  outfile.write(str(prov))                        # Write the angle info to the file

### print dihedrals

outfile.write("\nDihedrals\n\n")

for i in range(1,len(dihedrals)+1):               # Arrange the dihedral info
  prov = ('{x0} {x1} {x2} {x3} {x4} {x5}\n').format( \
      x0 = str(int(dihedrals[i-1,0])), \
      x1 = str(int(dihedrals[i-1,1])), \
      x2 = str(int(dihedrals[i-1,2])), \
      x3 = str(int(dihedrals[i-1,3])), \
      x4 = str(int(dihedrals[i-1,4])), \
      x5 = str(int(dihedrals[i-1,5])))

  outfile.write(str(prov))                        # Print the dihedral info to the file

outfile.close()
