#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 12 14:07:45 2018

@author: btreece
"""
from __future__ import print_function
import numpy as np
import igor.igorpy as igor
#import igor
import MDAnalysis

class Radii:
    #########################
    # CLASS INTITIALIZATION #
    #########################
    def __init__(self, Radius_Dictionary = None):
        if not Radius_Dictionary:
            self.Radius_Dictionary = {'H': 1.2, 'C': 1.7, 'N': 1.55, 'O': 1.52,
                                 'Na': 2.27, 'Mg': 1.73, 'P': 1.8, 'S': 1.8,
                                 'Cl': 1.75, 'K': 2.75, 'Ca': 2.31}
        else:
            self.Radius_Dictionary = Radius_Dictionary
            
        self.Mass_Dictionary = {'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999, 'S': 32.06}
        
        # Trivial Interpretations Of Elements
        self.Atom_Type_To_Element_Dictionary = {'H': 'H', 'C': 'C', 'N': 'N',
                                                'O': 'O', 'SOD': 'Na', 'MG': 'Mg',
                                                'P': 'P', 'S': 'S', 'CLA': 'Cl',
                                                'K': 'K', 'CAL': 'Ca'}
        # Common Hydrogen Types
        Charmm_Dict_Hydrogens = {'HA': 'H', 'HA1': 'H', 'HA2': 'H', 'HA3': 'H',
                                 'HAL1': 'H', 'HAL2': 'H', 'HAL3': 'H',
                                 'HB': 'H', 'HB1': 'H', 'HB2': 'H', 'HB3': 'H',
                                 'HC': 'H', 'HCL': 'H',
                                 'HD1': 'H', 'HD2': 'H', 'HD3': 'H',
                                 'HD11': 'H', 'HD12': 'H', 'HD13': 'H',
                                 'HD21': 'H', 'HD22': 'H', 'HD23': 'H',
                                 'HE': 'H', 'HE1': 'H', 'HE2': 'H', 'HE3': 'H',
                                 'HE11': 'H', 'HE12': 'H', 'HE21': 'H', 'HE22': 'H',
                                 'HEL1': 'H', 'HEL2': 'H',
                                 'HF1': 'H', 'HF2': 'H',
                                 'HG': 'H', 'HG1': 'H', 'HG11': 'H', 'HG12': 'H', 'HG13': 'H',
                                 'HG2': 'H', 'HG21': 'H', 'HG22': 'H', 'HG23': 'H',
                                 'HH': 'H', 'HH2': 'H',
                                 'HH11': 'H', 'HH12': 'H', 'HH21': 'H', 'HH22': 'H',
                                 'HL': 'H',
                                 'HN': 'H',
                                 'HR1': 'H', 'HR2': 'H', 'HR3': 'H',
                                 'HS': 'H',
                                 'HT': 'H', 'HT1': 'H', 'HT2': 'H', 'HT3': 'H',
                                 'HZ': 'H', 'HZ1': 'H', 'HZ2': 'H', 'HZ3': 'H'}
        self.Atom_Type_To_Element_Dictionary.update(Charmm_Dict_Hydrogens)
        
        # Common Carbon Types
        Charmm_Dict_Carbons = {'CA': 'C', 'CAP': 'C',
                               'CB': 'C',
                               'CC': 'C', 'CC2': 'C', 'CC1A': 'C', 'CC1B': 'C',
                               'CD': 'C', 'CD1': 'C', 'CD2': 'C',
                               'CE': 'C', 'CE1': 'C', 'CE2': 'C', 'CE3': 'C',
                               'CEL1': 'C', 'CEL2': 'C',
                               'CG': 'C', 'CG1': 'C', 'CG2': 'C',
                               'CH2': 'C',
                               'CL': 'C',
                               'CN': 'C',
                               'COA': 'C',
                               'CP1': 'C', 'CP2': 'C', 'CP3': 'C',
                               'CPH1': 'C', 'CPH2': 'C', 'CPT': 'C',
                               'CS': 'C',
                               'CT': 'C', 'CT1': 'C', 'CT2': 'C', 'CT3': 'C',
                               'CTL1': 'C', 'CTL2': 'C', 'CTL3': 'C',
                               'CY': 'C',
                               'CZ': 'C', 'CZ2': 'C', 'CZ3': 'C',
                               'C3': 'C'}
        self.Atom_Type_To_Element_Dictionary.update(Charmm_Dict_Carbons)
        
        # Common Nitrogen Types
        Charmm_Dict_Nitrogens = {'NC': 'N', 'NC2': 'N',
                                 'ND1': 'N', 'ND2': 'N',
                                 'NE': 'N', 'NE1': 'N', 'NE2': 'N',
                                 'NH1': 'N', 'NH2': 'N', 'NH3': 'N', 'NH3L': 'N',
                                 'NP': 'N', 'NPH': 'N',
                                 'NR1': 'N', 'NR2': 'N', 'NR3': 'N',
                                 'NTL': 'N',
                                 'NY': 'N',
                                 'NZ': 'N'}
        self.Atom_Type_To_Element_Dictionary.update(Charmm_Dict_Nitrogens)

        # Common Oxygen Types
        Charmm_Dict_Oxygens = {'OB': 'O',
                               'OC': 'O', 'OCA': 'O',
                               'OD1': 'O', 'OD2': 'O',
                               'OE1': 'O', 'OE2': 'O',
                               'OG': 'O', 'OG1': 'O',
                               'OH': 'O', 'OH1': 'O',
                               'OM': 'O',
                               'OS': 'O', 'OST': 'O',
                               'OT': 'O', 'OT1': 'O', 'OT2': 'O',
                               'OW': 'O', 'OWT': 'O', 'OWT1': 'O', 'OWT2': 'O', 'OWT3': 'O'}
        self.Atom_Type_To_Element_Dictionary.update(Charmm_Dict_Oxygens)
        
        # Common Sulfur Types
        Charmm_Dict_Sulfurs = {'SD': 'S', 'SG': 'S'}
        self.Atom_Type_To_Element_Dictionary.update(Charmm_Dict_Sulfurs)
        
    ################ THESE NEED TESTED #########################
    def Add_Entry_To_Current_Atom_Type_To_Element_Dictionary(self, Dictionary_To_Be_Added):
        for element in Dictionary_To_Be_Added.itervalues():
            if element not in self.Radius_Dictionary:
                print("Element '" + element + "' not found in radius dictionary. Nothing added.")
                break
        self.Atom_Type_To_Element_Dictionary.update(Dictionary_To_Be_Added)
    
    def Add_Radius_To_Dictionary(self, Dictionary_To_Be_Added):
        self.Radius_Dictionary(Dictionary_To_Be_Added)
        
    def Add_Mass_To_Dictionary(self, Dictionary_To_Be_Added):
        self.Mass_Dictionary(Dictionary_To_Be_Added)
    ############################################################
        
    def What_Are_The_Possible_Atom_Types_Of_An_Element(self, element):
        for atom_type in self.Atom_Type_To_Element_Dictionary.keys():
            if element == self.Atom_Type_To_Element_Dictionary[atom_type]:
                # Python 2.7
#                print(atom_type + ',\t',)
                # Python 3
                print("'" + atom_type + "',\t", end = '')
                    
    ###########################################################################
    # CONSTRUCT LIST OF VAN DER WAALS RADII FOR THE ATOMS IN AN ATOMSELECTION #
    ###########################################################################    
    def Construct_Radii_List(self, atomselection):
        
        ### Check the masses of the atoms against the masses of the corresponding element in dictionary ###
        for name, mass in zip(atomselection.names, atomselection.masses):
            if name in self.Atom_Type_To_Element_Dictionary:
                element = self.Atom_Type_To_Element_Dictionary[name]
                if element not in self.Mass_Dictionary:
                    print("Mass dictionary does not have the element associated with atom type '" + name + "'")
                else:
                    if (round(mass,3) != self.Mass_Dictionary[element]):
                        print("Element '" + element + "' has incorrect mass: {:}".format(round(mass,3)) +" *** The atom type is '" + name + "', check dictionary.")
                        break

        ### Check that all the atomtypes are in the dictionary ###
        for name, mass in zip(atomselection.names, atomselection.masses):
            if name not in self.Atom_Type_To_Element_Dictionary:
                print("Unknown Atom Type '" + name + "' with mass {:}".format(round(mass,3)))
                break
        
        radii = []
        for name in atomselection.names:
            element = self.Atom_Type_To_Element_Dictionary[name]
            radii.append(self.Radius_Dictionary[element])
        return np.array(radii)

##############################
### Protein-Type Densities ###
##############################

class Density_Profile:
    def __init__(self, isneutron = False):
        self.isneutron = isneutron

    def convert_units(self, length_scale_factor, new_units):
        """length_scale_factor is the multiplicative factor that takes the old length to the new one."""
        self.units = new_units
        
        for i in ['zmin', 'zmax', 'zstep', 'mean', 'second_moment', 'norm',
                  'tracking_of_mean', 'radii', 'bilayer_center',
                  'total_density_norm']:
            if hasattr(self, i):
                setattr(self, i, getattr(self, i)*length_scale_factor)
        
        for i in ['density', 'psigma', 'msigma']:
            if hasattr(self, i):
                setattr(self, i, getattr(self, i)/length_scale_factor)
        
        if hasattr(self, 'volume'):
            self.volume = self.volume*(length_scale_factor**3.0)

        if hasattr(self, 'potential_scaling'):
            self.potential_scaling = np.array([self.potential_scaling[0]*length_scale_factor, self.potential_scaling[1]/length_scale_factor**2.0])

        if hasattr(self, 'density_dictionary'):
            for key in self.density_dictionary.keys():
                self.density_dictionary[key] = self.density_dictionary[key]/length_scale_factor
        
##################################
# Profile From Neutron Data File #
##################################
class Protein_From_PXP(Density_Profile):
    def __init__(self, PXP_File = None, data_column_title = 'median_area',
                 msigma_column_title = 'msigma', psigma_column_title = 'psigma',
                 units = 'A'):
        self.isneutron = True
        self.filename = PXP_File
        self.data_column_title = data_column_title
        self.msigma_column_title = msigma_column_title
        self.psigma_column_title = psigma_column_title
        self.units = units
        if PXP_File:
            self.import_pxp_file(PXP_File, data_column_title = data_column_title,
                                 msigma_column_title = msigma_column_title,
                                 psigma_column_title = psigma_column_title,
                                 units = units)
        
    def import_pxp_file(self, filename = None, data_column_title = None,
                        msigma_column_title = None, psigma_column_title = None,
                        units='A', include_confidence = True):
        if ((not filename) and (not self.filename)):
            raise Exception('No filename given to this object or function.')
        elif (filename):
            self.filename = filename

        if (data_column_title):
            self.data_column_title = data_column_title
        if (msigma_column_title):
            self.msigma_column_title = msigma_column_title
        if (psigma_column_title):
            self.psigma_column_title = psigma_column_title
        
        f = igor.load(self.filename)
        
        # Find the z-axis data
        self.zmin = f.zaxis.data[0]
        self.zmax = f.zaxis.data[-1]
        self.zstep = f.zaxis.data[1] - f.zaxis.data[0]
        
        # Locate protein median data
        try:
            self.density = np.array(eval('f.'+self.data_column_title+'.data'))
        except:
            raise Exception('Data for ' + self.data_column_title + ' was not found in '+ self.filename)
            return

        # Locate confidence intervals
        if include_confidence:
            if not (hasattr(f, self.msigma_column_title) and hasattr(f, self.psigma_column_title)):
                raise Exception("Could not find one or both of '" + self.msigma_column_title
                                + "' and '" + self.psigma_column_title + "' in '" + self.filename + "'")
            self.msigma = np.array(eval('f.' + self.msigma_column_title + '.data'))
            self.psigma = np.array(eval('f.' + self.psigma_column_title + '.data'))

        # Check the array lengths for inconsistencies
        if round( (self.zmax - self.zmin)/self.zstep, 8)%1.0 != 0.0:
            raise Exception('zmin and zmax are not separated by an integral multiple of equal zsteps, problem with the file?')
            return
        elif round( ( ( (self.zmax - self.zmin) / self.zstep ) + 1 ), 8) != np.shape(self.density)[0]:
            raise Exception('The number of points in the first dimension of the density array does not match the z-dimensions supplied.\n    This is a weird error for neutron data.')
            return
        
        self.units = units
        print('\n\n\nUnits in the pxp are assumed to be '+self.units+""". Please convert the units if this is incorrect.""")
        print(""" Try 'instance_name'.convert_units(1.0, 'desired_units')\n\n\n""")
        
        # Normalize the density provided
        norm = sum(self.density)*self.zstep
        self.norm = norm
        if norm != 1.0:
            print('\n\n\nThe profile has area = {:}'.format(norm)+""", changing that to be 1.\n The normalization is stored in 'instance_name'.norm\n\n\n""")
            self.density = self.density / norm
            if include_confidence:
                self.msigma = self.msigma / norm
                self.psigma = self.psigma / norm
        
        # Calculate the mean and width of the profile
        self.mean = sum([ self.density[i]*(self.zmin + i*self.zstep) for i in range(len(self.density)) ])*self.zstep
        square_sum = sum([ self.density[i]*(self.zmin + i*self.zstep - self.mean)**2.0 for i in range(len(self.density)) ])*self.zstep
        self.second_moment = np.sqrt(square_sum)

##############################################################
# Profile From GROMACS Configuration File (Used For Biasing) #
##############################################################
class Protein_From_Configuration(Density_Profile):
    def __init__(self, configuration_file = None, units = 'nm'):
        self.isneutron = True
        self.filename = configuration_file
        self.units = units
        if self.filename:
            self.import_configuration_file()
  ################################
    def import_configuration_file(self, filename = None, units = 'nm'):
        if ((not filename) and (not self.filename)):
            raise Exception('No filename given to this object or function.')
        elif (filename):
            self.filename = filename

        file_tmp = open(self.filename,'r')
        lines = file_tmp.readlines()
        file_tmp.close()
        atom_groups = []
        density = []
        weights = []
        radii = []

        for i in lines:
            if i.split()[0] == 'u':
                scaling = [float(j) for j in i.split()[1:]]
            if i.split()[0] == 'n':
                atom_groups.append([int(i.split()[1]), int(i.split()[2])])
            if i.split()[0] == 'p':
                zmin = float(i.split()[1])
                zmax = float(i.split()[2])
                zstep = float(i.split()[3])
            if i.split()[0] == 'd':
                tmp = i.split()[1:]
                for j in tmp:
                    density.append(float(j))
            if i.split()[0] == 'w':
                tmp = i.split()[1:]
                for j in tmp:
                    weights.append(float(j))
            if i.split()[0] == 'r':
                tmp = i.split()[1:]
                for j in tmp:
                    radii.append(float(j))
        
        self.atom_groups = np.array(atom_groups)
        self.potential_scaling = np.array(scaling)
        self.zmin = zmin
        self.zmax = zmax
        self.zstep = zstep
        self.density = np.array(density)
        self.weights = np.array(weights)
        self.radii = np.array(radii)
        
        # Check the lengths of the arrays and endpoints
        if round( (self.zmax - self.zmin)/self.zstep, 8)%1.0 != 0.0:
            raise Exception('zmin and zmax are not separated by an integral multiple of zsteps.')
            return
        elif round( ( ( (self.zmax - self.zmin) / self.zstep ) + 1 ), 8) != np.shape(self.density)[0]:
            raise Exception('The number of points in the first dimension of the density array does not match the z-dimensions supplied.\n    This class assumes zmax is the last entry of the z-array.')
            return

        # Normalize the density provided
        self.norm = sum(self.density)*self.zstep
        if self.norm != 1.0:
            print('\n\n\nThe profile has area = {:}'.format(self.norm)+""", changing that to be 1.\n The normalization is stored in 'instance_name'.norm\n\n\n""")           
            self.density = self.density / self.norm
        
        self.units = units
        print('\n\n\nUnits are assumed to be '+self.units+""". Please change the property 'units' if this is incorrect.\n\n\n""")
        
        # Calculate the mean and width of the profile
        self.mean = sum([ self.density[i]*(self.zmin + i*self.zstep) for i in range(len(self.density)) ])*self.zstep
        sq_sum = sum([ self.density[i]*(self.zmin + i*self.zstep - self.mean)**2.0 for i in range(len(self.density)) ])*self.zstep
        self.second_moment = np.sqrt(sq_sum)
        
        if ( np.shape(self.radii)[0] == 0 ):
            # If no radii for the atoms are supplied, set to 1A
            if self.units == 'nm':
                self.radii = 0.1*np.ones(self.atom_groups[0,1] - self.atom_groups[0,0] + 1)
            elif self.units == 'A':
                self.radii = np.ones(self.atom_groups[0,1] - self.atom_groups[0,0] + 1)
            else:
                print("\n\n\nUnknown units for default radius definition. Using 1.0, please adjust accordingly.\n\n\n")
                self.radii = np.ones(self.atom_groups[0,1] - self.atom_groups[0,0] + 1)
        self.volume = np.sum(self.radii**3.0)
        
######################################
# Profile From Simulation Trajectory #
#######################################
class Protein_From_Simulation(Density_Profile):
    def __init__(self, structure_file = None, trajectory_file = None, atomselection = None, reference_profile_from_configuration = None, units = 'A'):
        self.isneutron = False
        if (structure_file):
            self.structure_file = structure_file
        if (trajectory_file):
            self.trajectory_file = trajectory_file
        if (atomselection):
            self.structure_file = atomselection.universe.filename
            self.trajectory_file = atomselection.universe.trajectory.filename
            self.atomselection = atomselection
        if (reference_profile_from_configuration):
            if not isinstance(reference_profile_from_configuration, Protein_From_Configuration):
                raise Exception("Reference profile is not a 'Protein_From_Configuration' object.")
            if ((not hasattr(reference_profile_from_configuration, 'atom_groups'))
            or (not hasattr(reference_profile_from_configuration, 'density'))
            or (not hasattr(reference_profile_from_configuration, 'radii'))):
                raise Exception("Some attributes were not found in your reference profile, did you import it correctly?")
            self.reference_profile = reference_profile_from_configuration
        self.units = units
        
    def __verify_frames(self, frames):
        if ((not frames) and (not hasattr(self, 'frames'))):
            frames = input("Please provide an iterable containing which frames to look at in the trajectory: ")
        if (frames):
            self.frames = frames
            
    def __verify_zaxis_and_radii(self, reference_profile):
        if (reference_profile):
            if not isinstance(reference_profile, Protein_From_Configuration):
                raise Exception("Reference profile is not a 'Protein_From_Configuration' object.")
            if ((not hasattr(reference_profile, 'atom_groups'))
            or (not hasattr(reference_profile, 'density'))
            or (not hasattr(reference_profile, 'radii'))):
                raise Exception("Some attributes were not found in your reference profile, did you import it correctly?")
            self.reference_profile = reference_profile
        if hasattr(self, 'reference_profile'):
            self.atom_groups = self.reference_profile.atom_groups
            if self.reference_profile.units == 'A':
                length_scale = 1.0
            elif self.reference_profile.units == 'nm':
                length_scale = 10.0
            else:
                raise Exception("Unknown length scale in the reference profile.")
            self.zmin = length_scale*self.reference_profile.zmin
            self.zstep = length_scale*self.reference_profile.zstep
            self.radii = length_scale*self.reference_profile.radii
            self.volume = length_scale**3.0 * self.reference_profile.volume
        else:
            self.zmin = float(input("What is the start of the z-axis for the density (in Angstroms)? "))
            self.zstep = float(input("What is the step size along the z-axis for the density (in Angstroms)? "))
            tmp_1 = int(input("At what index do the atoms of interest start in the '.gro' file? "))
            tmp_2 = int(input("At what index do the atoms of interest end in the '.gro' file? "))
            self.atom_groups = np.array([[tmp_1, tmp_2]])
            self.radii = np.ones(self.atom_groups[0,1] - self.atom_groups[0,0] + 1)
            self.volume = np.sum(self.radii**3.0)
        
        # There are 10000 entries in the GROMACS z-axis array
        self.zmax = self.zmin + 9999.0*self.zstep
        
    def calculate_simulation_density(self):
        
        self.density_trajectory = []
        prot_dens = 0*np.arange(self.zmin, self.zmax + 0.5*self.zstep, self.zstep)
        # Currently, the averaging only happens over one molecule
        nrm = self.volume*(2*np.pi)**0.5
        
        tracking_of_mean = []
        tracking_of_second_moment = []

        for t in self.frames:
            self.atomselection.universe.trajectory[t]
            z_bbox = self.atomselection.universe.dimensions[2]
            prot_dens = 0*prot_dens
            
            positions = self.atomselection.positions
            for i in range(len(positions)):
                z = positions[i,2]
                sigma = self.radii[i]
                # The plus or minus 3 sigma slightly under-represents the area, A = 0.9989759
                if z < self.zmin + 3.0 * sigma:
                    z += z_bbox
                z_low = np.floor( (z - 3.0 * sigma - self.zmin) / self.zstep)
                z_hi  = np.ceil( (z + 3.0 * sigma - self.zmin) / self.zstep)
                for j in range(int(z_low), int(z_hi)):
                    # The integral of each exponential is equal to sigma, compounding to make sigma^3, normalized in nrm by self.volume
                    prot_dens[j] += (sigma**2.0)*np.exp(-0.5*(z - (self.zmin + self.zstep*j))**2.0/((sigma)**2.0)) / nrm
                    
            mean = sum([prot_dens[i]*(self.zmin + i*self.zstep) for i in range(len(prot_dens))])*self.zstep
            second_moment = sum([prot_dens[i]*(self.zmin + i*self.zstep - mean)**2.0 for i in range(len(prot_dens))])*self.zstep
            tracking_of_mean.append(mean)
            tracking_of_second_moment.append(second_moment)
            self.density_trajectory.append(np.array([i for i in prot_dens]))
        
        self.density = sum(self.density_trajectory)/len(self.density_trajectory)
        self.tracking_of_mean = tracking_of_mean
        self.tracking_of_second_moment = tracking_of_second_moment
        self.units = 'A'
        
        # Calculate the mean and width of the profile
        self.mean = sum([ self.density[i]*(self.zmin + i*self.zstep) for i in range(len(self.density)) ])*self.zstep
        second_moment = sum([ self.density[i]*(self.zmin + i*self.zstep - self.mean)**2.0 for i in range(len(self.density)) ])*self.zstep
        self.second_moment = np.sqrt(second_moment)
        
    def import_simulation_density_from_atom_selection(self, atomselection = None, frames = None, reference_profile_from_configuration = None, units = 'A'):
        if units != 'A':
            raise Exception('Sorry, MDAnalysis uses Angstroms. Prepare your z-axis as appropriate and try again.')
            return
        self.units = units
        
        # Figure out what group of atoms is being calculated on.
        if ( (not atomselection) and (not hasattr(self, 'atomselection')) ):
            raise Exception("No atomselection provided in arguments or object.")
        if (atomselection):
            if not isinstance(atomselection, MDAnalysis.core.groups.AtomGroup):
                raise Exception("Was not provided with an appropriate atomselection of type 'MDAnalysis.core.groups.AtomGroup'.")
            self.atomselection = atomselection
            
        # At what times are we looking?
        self.__verify_frames(frames)
            
        # Figure out the z-axis and radii
        self.__verify_zaxis_and_radii(reference_profile_from_configuration)
        
        self.calculate_simulation_density()

    def import_simulation_density_from_trajectory(self, structure_file = None, trajectory_file = None, frames = None, reference_profile_from_configuration = None, units = 'A'):
        if units != 'A':
            raise Exception('Sorry, MDAnalysis uses Angstroms. Please convert post import using class function. Adjust z-axis as appropriate.')
            return
        self.units = units
        
        # Make sure we can find the trajectory
        if ((not structure_file) and (not hasattr(self, 'structure_file'))):
            structure_file = input("Please input the location of the structure file: ")
        if (structure_file):
            self.structure_file = structure_file
            
        if ((not trajectory_file) and (not hasattr(self, 'trajectory_file'))):
            trajectory_file = input("Please input the location of the trajectory file: ")
        if (trajectory_file):
            self.trajectory_file = trajectory_file

        # At what times are we looking?
        self.__verify_frames(frames)
            
        # Figure out the z-axis and radii
        self.__verify_zaxis_and_radii(reference_profile_from_configuration)

        u = MDAnalysis.Universe(self.structure_file, self.trajectory_file)
        if len(self.atom_groups) > 1:
            raise Exception('Sorry, this module does not yet support multiple molecules.')
        # Strangely enough, 'bynum a:b' selects atoms by index in gro file a-1 to b-1
        self.atomselection = u.select_atoms('bynum {:}'.format(self.atom_groups[0][0])+':{:}'.format(self.atom_groups[0][1]))
        
        self.calculate_simulation_density()
        

    def average_density_over_trajectory(self, frame_subset = None):
        if (not hasattr(self, 'density_trajectory')):
            raise Exception("The density trajectory does not exist for this object.")
        elif len(self.density_trajectory) == 0:
            raise Exception("The density trajectory exists but has no entries.")
        
        if not frame_subset:
            print("No subset of the frames specified, using all the frames.")
            frame_subset = range(len(self.frames))
        self.frame_subset = frame_subset
        
        self.density = sum([self.density_trajectory[i] for i in frame_subset]) / len(frame_subset)
        self.mean = sum([ self.density[i]*(self.zmin + i*self.zstep) for i in range(len(self.density)) ])*self.zstep
        sq_sum = sum([ self.density[i]*(self.zmin + i*self.zstep - self.mean)**2.0 for i in range(len(self.density)) ])*self.zstep
        self.second_moment = np.sqrt(sq_sum)
        

##############################
### Bilayer-Type Densities ###
##############################
class Bilayer_Profile(Density_Profile):
#class Bilayer_Profile:
    def __init__(self):
        POPC_RANGES={'headgroups':list(range(0,24)),'tails':list(range(24,130)),'methyls':list(range(130,134))}
        POPG_RANGES={'headgroups':list(range(0,17)),'tails':list(range(17,123)),'methyls':list(range(123,127))}
        DOPC_RANGES={'headgroups':list(range(0,24)),'tails':list(range(24,87))+list(range(91,134)),'methyls':list(range(87,91))+list(range(134,138))}
        DOPS_RANGES={'headgroups':list(range(0,17)),'tails':list(range(17,80))+list(range(84,127)),'methyls':list(range(80,84))+list(range(127,131))}
        CHOL_RANGES={'tails':list(range(0,74))}
        self.lipid_groups_dictionary = {'POPC':POPC_RANGES,'POPG':POPG_RANGES,'DOPC':DOPC_RANGES,'DOPS':DOPS_RANGES,'CHOL':CHOL_RANGES}
    
class Bilayer_From_PXP(Bilayer_Profile):
    def __init__(self, filename = None, protein_column_title='median_area',
                 group_dictionary = {}, units='A',
                 create_group_dictionary_interactively = False):
        Bilayer_Profile.__init__(self)
        self.isneutron = True
        self.filename = filename
        self.protein_column_title = protein_column_title
        self.group_dictionary = group_dictionary
        self.units = units
        if filename:
            self.import_pxp_file(create_group_dictionary_interactively=
                                     create_group_dictionary_interactively)

    def import_pxp_file(self, filename = None, protein_column_title = None,
                        group_dictionary = None, units = None,
                        create_group_dictionary_interactively = False):
        """ group_dictionary is a set of keys correspond to lists;
        each list contains the column headers in the pxp file for that group."""
        ##########        
        if filename:
            self.filename = filename
        elif not self.filename:
            raise Exception("No file name provided or set in attributes of object.")
        
        f = igor.load(self.filename)
        self.zmin = f.zaxis.data[0]
        self.zmax = f.zaxis.data[-1]
        self.zstep = f.zaxis.data[1] - f.zaxis.data[0]
        ##########
        if create_group_dictionary_interactively:
            group_dictionary = {}
            toplevel = True
            while toplevel:
                print("""Interactive Group Definition:\n1. Add A Group To The Dictionary\n2. List Columns In The File\n3. Exit""")
                answer = input("""Please Enter The Integer Of Your Choice - """)
                if answer == "1":
                    secondlevel = True
                    answer = input("""Enter The Group Key - """)
                    group_dictionary.update({answer:[]})
                    while secondlevel:
                        group = input("""Enter Column Title Or EXIT To Leave This Group - """)
                        if group == "EXIT":
                            secondlevel = False
                        elif not hasattr(f, group):
                            print("""\nThat Is Not One Of The Columns, Try Again\n""")
                        else:
                            group_dictionary[answer].append(group)
                elif answer == "2":
                    key_list = f.__dict__.keys()
                    key_string = ''
                    for key in key_list:
                        key_string += key
                        a = 25 - len(key)
                        if a < 1:
                            a = 1
                        for i in range(a):
                            key_string += ' '
                    key_string += '\n'
                    print(key_string)
                elif answer == "3":
                    toplevel = False
                else:
                    print("\nNot One Of The Options, Try Again\n")
        if group_dictionary:
            self.group_dictionary = group_dictionary
        
        if len(self.group_dictionary.keys()) == 0:
            self.group_dictionary = {'headgroups':['headgroup1','headgroup1_2','headgroup1_3','headgroup2','headgroup2_2','headgroup2_3'],
                                'tethers':['bME','tetherg','tether'],
                                'tails':['lipid1','methyl1','lipid2','methyl2'],
                                'substrate':['substrate']}
        if list(self.group_dictionary.keys()).count('tails') != 1:
            raise Exception("""A group named 'tails' (tail distribution) must be included for purposes of centering the bilayer.""")
            return
        ##########
        # Get the protein density
        if protein_column_title:
            self.protein_column_title = protein_column_title
        try:
            protein_density = np.array(eval('f.'+self.protein_column_title+'.data'))
        except:
            raise Exception(self.protein_column_title+' was not found in '+self.filename+'. I do need it for relative normalizing.')
            return
        ##########
        # Check the z-axis variables and the array lengths
        if round( (self.zmax - self.zmin)/self.zstep, 8)%1.0 != 0.0:
            raise Exception('zmin and zmax are not separated by an integral multiple of zsteps, problem with the file?')
            return
        elif round( ( ( (self.zmax - self.zmin) / self.zstep ) + 1 ), 8) != np.shape(protein_density)[0]:
            raise Exception('The number of points in the first dimension of the density array does not match the z-dimensions supplied.\n    This is a weird error for neutron data.')
            return
        ##########
        # Units
        if units:
            self.units = units
        ##########
        # Setting the dictionary of various densities in the file
        temp_dict = {}
        for i in self.group_dictionary.keys():
            temp_dict.update({i:0*protein_density})
            for j in self.group_dictionary[i]:
                try:
                    temp_dict.update( { i:temp_dict[i] + getattr( getattr( f, j), 'data') } )
                except:
                    raise Exception(j+' was not found in '+self.filename+'. Please check your dictionary of component groups for the correct column headings.')
                    return
        self.density_dictionary = temp_dict
        ##########
        # Get the norm of the protein density (this is used for relative normalizing in plots)
        self.protein_norm = sum(protein_density)*self.zstep
        ##########        
        # Calculate the center of the lipidic distribution
        norm = sum(self.density_dictionary['tails'])*self.zstep
        self.bilayer_center = sum([ self.density_dictionary['tails'][i]*(self.zmin + i*self.zstep) for i in range(len(self.density_dictionary['tails'])) ])*self.zstep / norm

class Bilayer_From_Simulation(Bilayer_Profile):
    def __init__(self, structure_file = None, trajectory_file = None, universe = None, frames = None,
                 lipid_resname_dictionary = {'CHOL':'CHL1', 'DOPC':'DOPC', 'DOPS':'DOPS'}):
        Bilayer_Profile.__init__(self)
        self.isneutron = False
        self.structure_file = structure_file
        self.trajectory_file = trajectory_file
        self.universe = universe
        self.frames = frames
        self.resname_dict = lipid_resname_dictionary

    def calculate_simulation_density(self, bilayer_selection_dictionary, frames):
        
        sigma = 1.0
        norm = len(frames)*sigma*(2*np.pi)**0.5
        
        temp_dict = {}
        key_list = list(bilayer_selection_dictionary.keys())
        z_array = np.arange( self.zmin, self.zmax, self.zstep)
        
        for key in key_list:
            temp_dict.update({ key:0*z_array })
        
        for t in frames:
            for i in range(len(key_list)):
                key = key_list[i]
                if i == 0:
                    bilayer_selection_dictionary[key].universe.trajectory[t]
                    
                for r in bilayer_selection_dictionary[key].positions:
                    z = r[2]
                    temp_dict.update({ key: temp_dict[key]+ np.exp(-0.5*(z - z_array)**2.0/((sigma)**2.0)) / norm })
        
        self.density_dictionary = temp_dict
        
        # Calculate the norm of all the densities together (nothing has been normalized accounting for number of atoms)
        norm = 0.0
        for key in key_list:
            norm += sum(temp_dict[key])*self.zstep
        self.total_density_norm = norm
        self.bilayer_center = sum([sum([temp_dict[key][i]*z_array[i] for key in key_list]) for i in range(len(z_array))])*self.zstep / self.total_density_norm
        
        self.bilayer_atomgroup_dictionary = bilayer_selection_dictionary
        self.frames = frames
        self.units = 'A'
        
        
        
    def import_simulation_density_from_universe(self, universe = None, frames = None, lipid_resname_dictionary = None, zmin = None, zmax = None, zstep = None):
        """ lipid_resname_dictionary gives how each component of the bilayer is labeled in the simulation - i.e. {'DOPC':'DOPC, 'CHOL':'CHL1'} """
        ##########
        if universe:
            self.universe = universe
        elif not self.universe:
            raise Exception("No universe supplied to method or as an attribute to object.")
        self.structure_file = self.universe.filename
        self.trajectory_file = self.universe.trajectory.filename
        ##########
        if frames:
            self.frames = frames
        elif not self.frames:
            raise Exception("No list of frames was supplied to method or as an attribute to object.")
        ##########
        if lipid_resname_dictionary:
            self.resname_dict = lipid_resname_dictionary
        ##########
        if zmin and zmax and zstep:
            self.zmin = zmin
            self.zmax = zmax
            self.zstep = zstep
        else:
            z_bbox = self.universe.dimensions[2]
            self.zmin = np.min(self.universe.atoms.positions[:,2]) - 0.1*z_bbox
            self.zmax = np.max(self.universe.atoms.positions[:,2]) + 0.1*z_bbox
            self.zstep = (self.zmax - self.zmin) / 2500.
        ##########
        # Lipid keys are things like 'CHOL' or 'DOPC'
        lipid_keys = self.resname_dict.keys()
        bilayer_selection_dictionary = {}
        
        for group in ['headgroups', 'tails', 'methyls']:
            # Create an 'empty selection' (requires version 0.17 MDAnalysis) to be unioned with other selections when each resname is analyzed
            temp_sel = self.universe.select_atoms('resname NONEXISTENT')
            for key in lipid_keys:
                if key not in self.lipid_groups_dictionary:
                    raise Exception(key+" was not found in the lipid_resname_dictionary.")
                if group in self.lipid_groups_dictionary[key]:
                    # Get the names of the atoms for the resname associated with key
                    NAMES=self.universe.select_atoms("resname "+self.resname_dict[key]).names
                
                    # name_range will have the indices of NAMES that correspond to the group of interest
                    name_range = self.lipid_groups_dictionary[key][group]
                    
                    # Start the string used for atom_selection
                    str_tmp = "name " + NAMES[name_range[0]]
                    for i in name_range[1:]:
                        str_tmp += " or name " + NAMES[i]
                        temp_sel = temp_sel.union(self.universe.select_atoms("resname "+self.resname_dict[key]+" and ("+str_tmp+")") )
            bilayer_selection_dictionary.update({group:temp_sel})
        self.bilayer_atomgroup_dictionary = bilayer_selection_dictionary
        
        self.calculate_simulation_density(self.bilayer_atomgroup_dictionary, self.frames)
    


    def import_simulation_density_from_trajectory(self, structure_file = None, trajectory_file = None, frames = None, lipid_resname_dictionary = None):
        
        ##########
        if structure_file:
            self.structure_file = structure_file
        elif not self.structure_file:
            raise Exception("No structure file supplied to method or as an attribute to object.")
        ##########
        if trajectory_file:
            self.trajectory_file = trajectory_file
        elif not self.trajectory_file:
            raise Exception("No trajectory file supplied to method or as an attribute to object.")
        ##########
        if frames:
            self.frames = frames
        elif not self.frames:
            raise Exception("No list of frames was supplied to method or as an attribute to object.")
        ##########
        if lipid_resname_dictionary:
            self.resname_dict = lipid_resname_dictionary
        ##########
        
        universe = MDAnalysis.Universe(structure_file, trajectory_file)
        self.import_simulation_density_from_universe(universe, frames, lipid_resname_dictionary)


        
#    def convert_units(self, length_scale_factor, new_units):
#        """length_scale_factor is the multiplicative factor that takes the old length to the new one."""
#        self.units = new_units
#        
#        for i in ['zmin', 'zmax', 'zstep', 'bilayer_center']:
#            if hasattr(self, i):
#                setattr(self, i, getattr(self, i)*length_scale_factor)
#        
#        if self.isneutron:
#            if hasattr(self, 'protein_norm'):
#                self.protein_norm = self.protein_norm * length_scale_factor
#        else:
#            if hasattr(self, 'density_dictionary'):
#                for key in self.density_dictionary.keys():
#                    self.density_dictionary[key] = self.density_dictionary[key]/length_scale_factor
