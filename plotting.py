#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 15 17:47:18 2018

@author: btreece
"""

from classes import Bilayer_Profile
from classes import Protein_Profile
from matplotlib import pyplot as plt

def Plot_Densities(List_of_Profiles, bilayer_color_dictionary = {}):
    if not isinstance(List_of_Profiles, list):
        raise ValueError(list)
        
    for Profile in List_of_Profiles:
        if not ( isinstance(Profile, Bilayer_Profile) or isinstance(Profile, Protein_Profile) ):
            raise Exception("""An element of 'List_of_Densities' is not an instance of 'Protein_Profile' or 'Bilayer_Profile'""")

    protein_norm = 1.0
    List_of_Bilayer_Profiles = [Profile for Profile in List_of_Profiles if isinstance(Profile, Bilayer_Profile)]
    if len(List_of_Bilayer_Profiles) > 0:
        if len(bilayer_color_dictionary.keys()) == 0:
            print("""You gave me no color options for the bilayer, get ready for a black bilayer.""")
        
        if len(List_of_Bilayer_Profiles) > 1:
            print("""\n\nThere are multiple instances of 'Bilayer_Profile', that is odd. Normalizations could get dicey.\n\n""")

        if hasattr(List_of_Bilayer_Profiles[0], 'protein_norm'):
            protein_norm = List_of_Bilayer_Profiles[0].protein_norm
        else:
            protein_norm = 0.0
    
    List_of_Protein_Profiles = [Profile for Profile in List_of_Profiles if isinstance(Profile, Protein_Profile)]
    list_of_protein_norms = []
    if len(List_of_Protein_Profiles) > 0:
        # If there are no Bilayer_Profiles, leave the Protein_Profiles normalized
        if protein_norm == 1.0:
            list_of_protein_norms = [1.0 for i in List_of_Protein_Profiles]
        # If the first Bilayer_Profile is from a pxp, use the normalization factor for the protein in that file. 
        elif protein_norm > 0.0:
            list_of_protein_norms = [protein_norm for i in List_of_Protein_Profiles]
        # If the first Bilayer_Profile is MD, try and use the number of atoms to normalize (MD Bilayers will have area equal to number of atoms in the bilayer)
        elif protein_norm == 0.0:
            if (len(List_of_Protein_Profiles) == 1) and ( not hasattr(List_of_Protein_Profiles[0], 'atom_groups') ):
                raise Exception("""\n\nI have no way to normalize the protein profile with respect to the bilayer.\n\nA workaround is to set the protein profile property 'atom_groups' to the number of atoms in the protein molecule.\n\n""")
            elif not any([hasattr(Profile, 'atom_groups') for Profile in List_of_Protein_Profiles]):
                raise Exception("""\n\nNone of the protein profiles have an associate number of atoms to normalize with respect to the bilayer. Plots would be very unhelpful.\n\n""")
            else:
                for Profile in List_of_Protein_Profiles:
                    if hasattr(Profile, 'atom_groups'):
                        a = Profile.atom_groups
                        list_of_protein_norms.append( a[0][1] - a[0][0] + 1 )
                    else:
                        list_of_protein_norms.append(0.0)
                # When there is no info about the number of atoms, assume there is one number (equivalent to the max of the list) that describes how many atoms are in that molecule.
                for i in range(len(list_of_protein_norms)):
                    if list_of_protein_norms[i] == 0.0:
                        list_of_protein_norms[i] = max(list_of_protein_norms)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    plot_data = []
    for Profile in List_of_Bilayer_Profiles:
        for key in Profile.density_dictionary.keys():
            if bilayer_color_dictionary.has_key(key):
                colour = bilayer_color_dictionary[key]
            else:
                colour = 'black'
            plot_data.append(ax.plot([Profile.zmin + i*Profile.zstep for i in range(len(Profile.density_dictionary[key]))], Profile.density_dictionary[key], color = colour)[0])
    for i in range(len(List_of_Protein_Profiles)):
        Profile = List_of_Protein_Profiles[i]
        norm = list_of_protein_norms[i]

        if hasattr(Profile, 'color'):
            colour = Profile.color
        elif Profile.isneutron:
            colour = 'red'
        else:
            colour = 'black'
        
        density = Profile.density * norm
        plot_data.append(ax.plot([Profile.zmin + i*Profile.zstep for i in range(len(density))], density, color = colour))
        if hasattr(Profile, 'msigma'):
            msigma = Profile.msigma * norm
            psigma = Profile.psigma * norm
            plot_data.append(ax.fill_between([Profile.zmin + i*Profile.zstep for i in range(len(msigma))], msigma, psigma, color = colour, alpha = 0.3))
    
    fig.show()

    return ax, plot_data