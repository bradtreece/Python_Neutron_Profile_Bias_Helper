#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 10 10:28:12 2018

@author: btreece
"""

from classes import Bilayer_Profile
from classes import Protein_Profile
import numpy as np

def write_configuration_file(Profile, filename, weight_func = None, existing_weights = False):
    if not isinstance(Profile, Protein_Profile):
        raise ValueError(Protein_Profile)
    if not Profile.isneutron:
        raise Exception('The profile provided is not a neutron profile, please provide a different profile.')
    if not callable(weight_func):
        def weight_func(width, density):
            if density != 0.0:
                return np.exp(-1.0*width/density)
            else:
                return 1.0
    if (not hasattr(Profile, 'msigma')) and (not existing_weights):
        raise Exception('The profile provided does not have confidence intervals, please add them or set them to an array of zeros the size of the density.')
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
        print('\n\n\nThe profile has area = '+`norm`+""", changing that to be 1.\n    You're welcome.\n\n\n""")
        Profile.density = Profile.density / norm
        Profile.msigma = Profile.msigma / norm
        Profile.psigma = Profile.psigma / norm

    
    f = open(filename,'w')
    
    tmp = Profile.potential_scaling
    f.write('u\t'+`tmp[0]`+'\t'+`tmp[1]`+'\n')
    
    for i in Profile.atom_groups:
        f.write('n\t'+`i[0]`+'\t'+`i[1]`+'\n')
    
    tmp = Profile.zmin
    tmp2 = Profile.zstep
    tmp3 = tmp + (len(Profile.density)-1)*tmp2
    f.write('p\t'+`tmp`+'\t'+`tmp3`+'\t'+`tmp2`+'\n')
    
    f.write('d')
    for i in Profile.density:
        f.write('\t'+`i`)
    f.write('\n')
    
    if existing_weights:
        if not hasattr(Profile, 'weights'):
            raise Exception('There were no weights supplied with the profile but I was told the weights existed.')
        f.write('w')
        for i in Profile.weights:
            f.write('\t'+`i`)
        f.write('\n')
    else:
        tmp = Profile.density
        tmp2 = 0.5*(Profile.psigma - Profile.msigma)
        f.write('w')
        for i in range(len(tmp)):
            w = weight_func(tmp2[i], tmp[i])
            f.write('\t'+`w`)
        f.write('\n')
    
    f.close()
    
def prepare_neutron_profile_for_configuration(protein_profile, Bilayer_Profile_Neutron, Bilayer_Profile_MD):
    if not isinstance(protein_profile, Protein_Profile):
        raise ValueError(Protein_Profile)
    if not isinstance(Bilayer_Profile_MD, Bilayer_Profile):
        raise ValueError(Bilayer_Profile)
    if not isinstance(Bilayer_Profile_Neutron, Bilayer_Profile):
        raise ValueError(Bilayer_Profile)
        
    if not protein_profile.isneutron:
        raise Exception('The protein profile provided is not a neutron profile, please provide a different profile.')
    if not Bilayer_Profile_Neutron.isneutron:
        raise Exception('The neutron bilayer profile provided is not a neutron profile, please provide a different profile.')
    if Bilayer_Profile_MD.isneutron:
        raise Exception('The MD bilayer profile provided is a neutron profile, please provide a different profile.')
        
    if (not hasattr(protein_profile, 'msigma')) and (not hasattr(protein_profile, 'weights')):
        raise Exception('The profile provided does not have confidence intervals, nor does it have weights. Please add one of them or set the weights to an array of ones the size of the density.')

    if protein_profile.units == 'A':
        protein_profile.convert_units(0.1, 'nm')
        print('Converted your protein profile from Angstroms to nm (not nautical miles)')
    elif protein_profile.units != 'nm':
        raise Exception('The given units of the protein profile are '+protein_profile.units+', please convert to nm.')
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
    norm = sum(protein_profile.density)*protein_profile.zstep
    if norm != 1.0:
        print('\n\n\nThe profile has area = '+`norm`+""", changing that to be 1.\n    You're welcome.\n\n\n""")
        protein_profile.density = protein_profile.density / norm
        if hasattr(protein_profile, 'msigma'):
            protein_profile.msigma = protein_profile.msigma / norm
            protein_profile.psigma = protein_profile.psigma / norm
    
    # Apply the necessary offset
    offset = Bilayer_Profile_MD.bilayer_center - Bilayer_Profile_Neutron.bilayer_center
    protein_profile.mean = protein_profile.mean + offset
    protein_profile.zmin = protein_profile.zmin + offset
    protein_profile.zmax = protein_profile.zmax + offset
    
    # Truncate the density
    dens_temp = np.array([i for i in protein_profile.density])
    index_lo = min(np.nonzero(dens_temp)[0]) - 10
    index_hi = max(np.nonzero(dens_temp)[0]) + 11
    if index_lo < 0:
        index_lo = 0
    if index_hi > len(dens_temp):
        index_hi = dens_temp
    
    protein_profile.zmin = protein_profile.zmin + index_lo*protein_profile.zstep
    protein_profile.zmax = protein_profile.zmax - index_hi*protein_profile.zstep
    
    protein_profile.density = protein_profile.density[index_lo:index_hi]
    if hasattr(protein_profile, 'msigma'):
        protein_profile.msigma = protein_profile.msigma[index_lo:index_hi]
        protein_profile.psigma = protein_profile.psigma[index_lo:index_hi]
    if hasattr(protein_profile, 'weights'):
        protein_profile.weights = protein_profile.weights[index_lo:index_hi]
    