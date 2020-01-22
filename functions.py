#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 10 10:28:12 2018

@author: btreece
"""

#from classes import Density_Profile
from classes import Protein_From_PXP
from classes import Protein_From_Configuration
from classes import Protein_From_Simulation
#from classes import Bilayer_Profile
from classes import Bilayer_From_PXP
from classes import Bilayer_From_Simulation
#from classes import *
import numpy as np

class Profile_Comparison:
    def __init__(self, nframes):
        if not isinstance(nframes, int):
            raise Exception("Need an integral number of frames.")
        self.rmsd_inst = -1.0*np.ones(nframes)
        self.rmsd_avg = -1.0*np.ones(nframes)

def compare_simulation_to_reference(Simulation_Profile, Reference_Profile=None, Subset_Of_Frames_To_Look_At=None):
    if ((not Reference_Profile) and (not hasattr(Simulation_Profile, 'reference_profile'))):
        raise Exception("No reference provided and none found in the simulation profile.")
    if not isinstance(Simulation_Profile, Protein_From_Simulation):
        raise Exception("Simulation not of class 'Protein_From_Simulation'")
    if not isinstance(Reference_Profile, Protein_From_Configuration):
        raise Exception("Only reference profiles from configuration files supported.")
    
    if not Subset_Of_Frames_To_Look_At:
        Subset_Of_Frames_To_Look_At = range(len(Simulation_Profile.frames))
    
    if Reference_Profile:
        Ref = Reference_Profile
    else:
        Ref = Simulation_Profile.reference_profile
    
    # Make the z-axis match
    if ((Ref.zmin != Simulation_Profile.zmin)
    or (Ref.zstep != Simulation_Profile.zstep)
    or (Ref.zmax != Simulation_Profile.zmax)):
        target_density = 0*np.arange(Simulation_Profile.zmin, Simulation_Profile.zmax+1e-6, Simulation_Profile.zstep)
        for i in range(len(target_density)):
            z = Simulation_Profile.zmin + i*Simulation_Profile.zstep
            j = int(np.floor((z - Ref.zmin) / Ref.zstep))
            if ((j < 0) or (j > len(Ref.density) - 2)):
                target_density[i] = 0.0
            else:
                zed = Ref.zmin + j*Ref.zstep
                delta_z = z - zed
                delta_rho = Ref.density[j+1] - Ref.density[j]
                target_density[i] = (delta_rho/Ref.zstep)*(delta_z) + Ref.density[j]
        
    Comparison = Profile_Comparison(len(Subset_Of_Frames_To_Look_At))
    
    avg_dens = 0*Simulation_Profile.density
    for i in range(len(Subset_Of_Frames_To_Look_At)):
        t = Subset_Of_Frames_To_Look_At[i]
        inst_dens = Simulation_Profile.density_trajectory[t]
        avg_dens = (i*avg_dens + inst_dens)/(i+1)
        
        Comparison.rmsd_inst[i] = np.sqrt(np.sum((target_density - inst_dens)**2.0))
        Comparison.rmsd_avg[i] = np.sqrt(np.sum((target_density - avg_dens)**2.0))
    return Comparison

def write_configuration_file(Profile, filename, use_radii = False):
#### Weighted potential not supported yet
#def write_configuration_file(Profile, filename, weight_func = None, existing_weights = False):
    if not isinstance(Profile, Protein_From_PXP):
        raise ValueError(Protein_From_PXP)
#### Weights not supported yet
#    if not callable(weight_func):
#        def weight_func(width, density):
#            if density != 0.0:
#                return np.exp(-1.0*width/density)
#            else:
#                return 1.0
#    if (not hasattr(Profile, 'msigma')) and (not existing_weights):
#        raise Exception('The profile provided does not have confidence intervals, please add them or set them to an array of zeros the size of the density.')
    if Profile.units != 'nm':
        raise Exception('The given units are '+Profile.units+', please convert to nm.')
        return
    if not hasattr(Profile, 'atom_groups'):
        raise Exception('No atom_groups were defined. Please set this attribute.')
    if not hasattr(Profile, 'potential_scaling'):
        raise Exception('No potential_scaling was defined. Please set this attribute.')
    
    # Normalize the density provided
    norm = sum(Profile.density)*Profile.zstep
    if norm != 1.0:
        print('\n\n\nThe profile has area = {:}'.format(norm)+""", changing that to be 1.\n    You're welcome.\n\n\n""")
        Profile.density = Profile.density / norm
        Profile.msigma = Profile.msigma / norm
        Profile.psigma = Profile.psigma / norm

    
    f = open(filename,'w')
    
    tmp = Profile.potential_scaling
    f.write('u\t{:}\t{:}\n'.format(tmp[0],tmp[1]))
    
    for i in Profile.atom_groups:
        f.write('n\t{:}\t{:}\n'.format(i[0],i[1]))
    
    tmp = Profile.zmin
    tmp2 = Profile.zstep
    tmp3 = tmp + (len(Profile.density)-1)*tmp2
    f.write('p\t{:}\t{:}\t{:}\n'.format(tmp,tmp3,tmp2))
    
    f.write('d')
    for i in Profile.density:
        f.write('\t{:}'.format(i))
    f.write('\n')

###### Weights for potential, not supported yet
#    if existing_weights:
#        if not hasattr(Profile, 'weights'):
#            raise Exception('There were no weights supplied with the profile but I was told the weights existed.')
#        f.write('w')
#        for i in Profile.weights:
#            f.write('\t'+`i`)
#        f.write('\n')
#    else:
#        tmp = Profile.density
#        tmp2 = 0.5*(Profile.psigma - Profile.msigma)
#        f.write('w')
#        for i in range(len(tmp)):
#            w = weight_func(tmp2[i], tmp[i])
#            f.write('\t'+`w`)
#        f.write('\n')
    if use_radii:
        f.write('r')
        for i in Profile.radii:
            f.write('\t{:}'.format(i))
            f.write('\n')
    
    f.close()
    
def prepare_neutron_profile_for_configuration(Protein_Profile, Bilayer_Profile_Neutron, Bilayer_Profile_MD):
    if not isinstance(Protein_Profile, Protein_From_PXP):
        raise ValueError(Protein_From_PXP)
    if not isinstance(Bilayer_Profile_MD, Bilayer_From_Simulation):
        raise ValueError(Bilayer_From_Simulation)
    if not isinstance(Bilayer_Profile_Neutron, Bilayer_From_PXP):
        raise ValueError(Bilayer_From_PXP)
        
#    if (not hasattr(Protein_Profile, 'msigma')) and (not hasattr(Protein_Profile, 'weights')):
#        raise Exception('The profile provided does not have confidence intervals, nor does it have weights. Please add one of them or set the weights to an array of ones the size of the density.')

    if Protein_Profile.units == 'A':
        Protein_Profile.convert_units(0.1, 'nm')
        print('Converted your protein profile from Angstroms to nm (not nautical miles)')
    elif Protein_Profile.units != 'nm':
        raise Exception('The given units of the protein profile are '+Protein_Profile.units+', please convert to nm.')
    if Bilayer_Profile_MD.units == 'A':
        Bilayer_Profile_MD.convert_units(0.1, 'nm')
        print('Converted your MD bilayer profile from Angstroms to nm')
    elif Bilayer_Profile_MD.units != 'nm':
        raise Exception('The given units of the MD bilayer profile are '+Bilayer_Profile_MD.units+', please convert to nm.')
    if Bilayer_Profile_Neutron.units == 'A':
        Bilayer_Profile_Neutron.convert_units(0.1, 'nm')
        print('Converted your neutron bilayer profile from Angstroms to nm')
    elif Bilayer_Profile_Neutron.units != 'nm':
        raise Exception('The given units of the neutron bilayer profile are '+Bilayer_Profile_Neutron.units+', please convert to nm.')
    
    # Normalize the density provided
    norm = sum(Protein_Profile.density)*Protein_Profile.zstep
    if norm != 1.0:
        print('\n\n\nThe profile has area = {:}'.format(norm)+""", changing that to be 1.\n    You're welcome.\n\n\n""")
        Protein_Profile.density = Protein_Profile.density / norm
        if hasattr(Protein_Profile, 'msigma'):
            Protein_Profile.msigma = Protein_Profile.msigma / norm
            Protein_Profile.psigma = Protein_Profile.psigma / norm
    
    # Apply the necessary offset
    offset = Bilayer_Profile_MD.bilayer_center - Bilayer_Profile_Neutron.bilayer_center
    Protein_Profile.mean = Protein_Profile.mean + offset
    Protein_Profile.zmin = Protein_Profile.zmin + offset
    Protein_Profile.zmax = Protein_Profile.zmax + offset
    
    # Truncate the density
    dens_temp = np.array([i for i in Protein_Profile.density])
    index_lo = min(np.nonzero(dens_temp)[0]) - 10
    index_hi = max(np.nonzero(dens_temp)[0]) + 11
    if index_lo < 0:
        index_lo = 0
    if index_hi > len(dens_temp):
        index_hi = dens_temp
    
    Protein_Profile.zmin = Protein_Profile.zmin + index_lo*Protein_Profile.zstep
    Protein_Profile.zmax = Protein_Profile.zmax - index_hi*Protein_Profile.zstep
    
    Protein_Profile.density = Protein_Profile.density[index_lo:index_hi]
    if hasattr(Protein_Profile, 'msigma'):
        Protein_Profile.msigma = Protein_Profile.msigma[index_lo:index_hi]
        Protein_Profile.psigma = Protein_Profile.psigma[index_lo:index_hi]
    if hasattr(Protein_Profile, 'weights'):
        Protein_Profile.weights = Protein_Profile.weights[index_lo:index_hi]
    
