#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 15 17:47:18 2018

@author: btreece
"""

from .classes import Density_Profile
from .classes import Protein_From_PXP
from .classes import Protein_From_Configuration
from .classes import Protein_From_Simulation
from .classes import Bilayer_Profile
from .classes import Bilayer_From_PXP
from .classes import Bilayer_From_Simulation
#from classes import *
from matplotlib import pyplot as plt

def Plot_Densities(List_of_Profiles, bilayer_color_dictionary = {}):
    if not isinstance(List_of_Profiles, list):
        raise ValueError(list)
        
    for Profile in List_of_Profiles:
        if not ( isinstance(Profile, Density_Profile) ):
            raise Exception("""An element of 'List_of_Densities' is not a Denisty_Profile""")

    List_of_Bilayer_Profiles = [Profile for Profile in List_of_Profiles if isinstance(Profile, Bilayer_Profile)]
    List_of_Protein_Profiles = [Profile for Profile in List_of_Profiles if not isinstance(Profile, Bilayer_Profile)]
    
    if len(List_of_Bilayer_Profiles) > 0:
        if len(bilayer_color_dictionary.keys()) == 0:
            print("""You gave me no color options for the bilayer, get ready for a black bilayer.""")
        
        if len(List_of_Bilayer_Profiles) > 1:
            print("""\n\nThere are multiple instances of 'Bilayer_Profile', that is odd. Normalizations could get dicey.\n\n""")
    
    list_of_protein_norms = []
    for i in List_of_Protein_Profiles:
        if hasattr(i, 'atom_groups'):
            list_of_protein_norms.append(i.atom_groups[0,1] - i.atom_groups[0,0] + 1)
        else:
            list_of_protein_norms.append(0.0)
    
    tmp = [i for i in list_of_protein_norms if i!=0.0]
    if len(tmp) > 0:
        natoms = tmp[0]
    elif len(list_of_protein_norms) > 0:
        natoms = 1.0
    else:
        natoms = -1.0

    list_of_bilayer_norms = []
    for i in List_of_Bilayer_Profiles:
        if isinstance(i, Bilayer_From_PXP):
            list_of_bilayer_norms.append(i.protein_norm)
        elif isinstance(i, Bilayer_From_Simulation):
#            tmpatoms = sum([i.bilayer_atomgroup_dictionary[key].n_atoms for key in i.bilayer_atomgroup_dictionary.keys()])
            list_of_bilayer_norms.append((1.0*natoms))
        else:
            print("Couldn't construct a norm for one of the bilayers, leaving the scaling alone on this profile.")
            list_of_bilayer_norms.append(1.0)


    fig = plt.figure()
    ax = fig.add_subplot(111)
    plot_data = []
    for n in range(len(List_of_Bilayer_Profiles)):
        Profile = List_of_Bilayer_Profiles[n]
        if natoms == -1.0:
            norm = 1.0
        else:
            norm = list_of_bilayer_norms[n]

        for key in Profile.density_dictionary.keys():
            if key in bilayer_color_dictionary:
                colour = bilayer_color_dictionary[key]
            else:
                colour = 'black'
            plot_data.append(ax.plot([Profile.zmin + i*Profile.zstep for i in range(len(Profile.density_dictionary[key]))], Profile.density_dictionary[key] / norm, color = colour)[0])
    for n in range(len(List_of_Protein_Profiles)):
        Profile = List_of_Protein_Profiles[n]
        
        if hasattr(Profile, 'color'):
            colour = Profile.color
        elif Profile.isneutron:
            colour = 'red'
        else:
            colour = 'black'
        
        plot_data.append(ax.plot([Profile.zmin + i*Profile.zstep for i in range(len(Profile.density))], Profile.density, color = colour))
        if hasattr(Profile, 'msigma'):
            msigma = Profile.msigma
            psigma = Profile.psigma
            plot_data.append(ax.fill_between([Profile.zmin + i*Profile.zstep for i in range(len(msigma))], msigma, psigma, color = colour, alpha = 0.3))
    
    fig.show()

    return ax, plot_data