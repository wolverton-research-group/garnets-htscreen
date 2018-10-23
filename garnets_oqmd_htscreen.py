
# The code below can be used to generate OQMD derived thermodynamic quantities 
# for the Li3-garnets predicted by Aykol, Kim, Hegde and Wolverton (2018).
# Please note that this code requires the qmpy library, which can be
# accessed at https://github.com/wolverton-research-group/qmpy and also the OQMD
# installed locally, per the instructions provided therein. Present authors used the
# version 1.1 of the database as available for download at oqmd.org/download.


import matplotlib.pyplot as plt
import numpy as np
from qmpy import *

garnets = Formation.objects.filter(entry__path__contains='garnet')

print 'Garnet, formation_energy stability band_gap V_max V_min stable_phases'

for garnet in garnets:   
    
    # Calculating stability with respect to other materials in OQMD within the 
    # garnet's chemical space
    space = PhaseSpace('-'.join(garnet.entry.space))
    phase_O = Phase(composition='O',energy=0)
    phase_O.energy = -0.31695794 #298 K -TS of oxygen gas at 1 atm
    space.phase_dict['O'] = phase_O
    energy, phases = space.gclp(garnet.entry.name)
    energy = energy/20.0
    space.clear_all()
    stability = garnet.delta_e-energy
    if stability > 0.05:
        continue
    
    # COMPUTATION OF BULK VOLTAGE PROFILES
    # Taking a slice of the convex hull from the garnet towards Li.
    
    base_comp = str(garnet.entry.name[3:])
    slice_comp = base_comp + "-Li15"+base_comp
    ps = PhaseSpace(slice_comp)
    phase_O = Phase(composition='O',energy=0)
    phase_O.energy = -0.31695794
    ps.phase_dict['O'] = phase_O
    ps.get_hull()
    
    # List of phases on the convex hull slice
    eq_list = list(ps.hull)
    
    # Voltages
    V_list = []
   
    
    # We keep track of oxygen evolution as well.
    O_x, O_y = [], []
    
    # Scanning the sets of equilibrium phases through the slice
    for eq in eq_list:   
        energy1 = eq.phases[0].energy/eq.phases[0].unit_comp['O']*12.0
        energy2 = eq.phases[1].energy/eq.phases[1].unit_comp['O']*12.0
        if 'Li' in eq.phases[0].unit_comp:
            x_Li1 = eq.phases[0].unit_comp['Li']/eq.phases[0].unit_comp['O']*12.0
        else:
            x_Li1 = 0
        if 'Li' in eq.phases[1].unit_comp:
            x_Li2 = eq.phases[1].unit_comp['Li']/eq.phases[1].unit_comp['O']*12.0
        else:
            x_Li2 = 0
        V = -(energy2-energy1)/(x_Li2-x_Li1)
        V_list.append([x_Li1,V])
        V_list.append([x_Li2,V])

        # Check if oxygen is one of the phases (to record evolution). 
        for phase in eq.phases[0].phase_dict:
            if phase.name == 'O': 
                O_x.append(x_Li1)
                O_y.append(V)
        for phase in eq.phases[1].phase_dict:
            if phase.name == 'O':
                O_x.append(x_Li2)
                O_y.append(V)
            
    V_list.sort(key=lambda e: e[1], reverse=True)
    V_list.sort(key=lambda e: e[0])
    
    window_V = []
    
    # For stable materials, get the voltage window
    # For metastable materials, get the average voltage at x=3.0
    if np.allclose(stability, 0.0, atol=0.0001):
        for i in range(len(V_list)):
            if np.allclose(V_list[i][0], 3.0, atol=0.0001):
                window_V.append(V_list[i][1])
    else:
        for i in range(len(V_list)):
            if np.allclose(V_list[i][0], 3.0, atol=0.0001):
                window_V.append((V_list[i][1]+V_list[i+1][1])/2.0)
                break
            elif V_list[i][0]>3.0:
                window_V.append(V_list[i][1])
                break
                
    
    x_data, y_data = [], []
    
    arrow_y = 3.0 
    for i in xrange(len(V_list)):
        x_data.append(V_list[i][0])
        y_data.append(V_list[i][1])
        if V_list[i][0]==3.0:
            arrow_y=V_list[i][1]
            
    # Plotting the voltage profile
    plt.plot(x_data,y_data, linewidth=2.0, color='red')
    plt.xlabel(r'x in Li$\mathsf{_x}$'+base_comp,fontsize=14)
    plt.ylabel(r'Voltage vs Li/Li+ (V)',fontsize=14)
    plt.title('Lix'+base_comp,fontsize=14,fontweight='bold')
    plt.axis([0,16,0,5.5])
    
    # Annotating the garnet composition on the plot
    string = ' '.join([garnet.entry.latex, 
                      r'$E_g =$ ', str(garnet.calculation.band_gap), 
                      r' eV'])
    plt.annotate(string, xy=(3.25, arrow_y+0.1), xytext=(3.6,arrow_y+0.6), fontsize=14,
                  arrowprops=dict(arrowstyle="-|>", relpos=(0.1, 0.), fc='b') )
    
    # Find the lowest ooxygen release potential
    O_volts = 1000.0
    i_volts = -1
    for i in xrange(len(O_y)):
        if O_y[i] < O_volts:
            i_volts = i
            O_volts = O_y[i]
            
    # Annotate where O is released.
    if i_volts != -1:
        plt.annotate(r'$\mathrm{O_2-released}$', xy=(O_x[i_volts]+0.25, =O_y[i_volts]+0.1), fontsize=14, 
                     xytext=(O_x[i_volts]+0.6,O_y[i_volts]+0.6),
                     arrowprops=dict(arrowstyle='-|>',fc='b'))

    plt.tight_layout()
    plt.savefig('Li3'+base_comp+'.png')
    plt.close()
    
    print garnet.entry.name + " {0:.4g}".format(garnet.delta_e) \
                + " {0:.4g}".format(stability) + " " + str(garnet.calculation.band_gap) + " " + \
                " {0:.4g}".format(max(window_V)), " {0:.4g}".format(min(window_V)), str(phases)


    ps.clear_all()

# It is possible to access the structure (POSCAR) or 
# other input files (except POTCARs) as demonstrated below:

for i in garnets:
    if i.entry.name=='Li3Nd3Te2O12':
        print i.calculation.POSCAR
        print i.calculation.INCAR