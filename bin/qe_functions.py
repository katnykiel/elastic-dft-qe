def read_key():
    """
    Read in new Materials Project API key
    """
    import os, stat
    from IPython.display import clear_output

    # Read in new Materials Project API key, if one exists
    try:
        with open(os.path.expanduser('~/.mpkey.txt'), 'r') as f:
            key = f.readlines()[0]
            return key
    except:
        key = ""

    # Check if API key already exists, skip try-except
    if not key:
        # Prompt user for API key
        try:
            user = str(input())
            clear_output()
            if not user.isalnum():
                raise TypeError('Wrong Key')
            if user == None:
                raise TypeError('Empty')
            with open(os.path.expanduser('~/.mpkey.txt'), 'w') as keyfile:
                keyfile.write(user)
            os.chmod(os.path.expanduser('~/.mpkey.txt'), stat.S_IREAD | stat.S_IWRITE)
            del user

            with open(os.path.expanduser('~/.mpkey.txt'),'r') as f:
                key = f.readlines()[0]
                return key
            print("Success")
        except:
            print("Something seems wrong with your key")
            
def view_struct(struct):
    """
    Plot interactive 3D crystal structure 
    
    input: 
        struct (pymatgen structure object)
    output: 
        n/a
    """
    
    # Import libraries
    import plotly.express as px
    import plotly.graph_objects as go
    import pandas as pd
    import numpy as np
    
    # Convert list of sites to pandas DataFrame object (there's gotta be a better way?)
    x_list = [site.coords[0] for site in struct.sites]
    y_list = [site.coords[1] for site in struct.sites]
    z_list = [site.coords[2] for site in struct.sites]
    e_list = [site.species.elements[0] for site in struct.sites]
    
    scaling_factor = 20
    r_list = [scaling_factor*site.species.elements[0].atomic_radius_calculated 
              for site in struct.sites]
    
    site_df = pd.DataFrame({'x':x_list,'y':y_list,'z':z_list,'element':e_list,'r':r_list})
    
    # Draw spheres at each site
    fig = px.scatter_3d(site_df,x='x',y='y',z='z',size='r',size_max=100,
                        color='element',opacity=1,
                        color_discrete_sequence=px.colors.qualitative.Bold)

    # Convert lattice parameters to unit cell box
    lattice = struct.lattice.matrix
    corners = np.array([[0,1,1,0,0,0,0,1,1,0,0,1,1,1,1,0,0],
                        [0,0,1,1,0,0,0,0,1,1,0,0,0,1,1,1,1],
                        [0,0,0,0,0,1,1,1,1,1,1,1,0,0,1,1,0]]).T
    cell = np.matmul(corners,lattice).transpose()
    box = pd.DataFrame({'x':cell[0],'y':cell[1],'z':cell[2]})
    
    # Draw box
    box = px.line_3d(box,x='x', y='y', z='z')
    box.data[0]['line']['color']='black'
    
    fig.add_traces(list(box.select_traces()))
    fig.update_layout(scene = dict(
                      xaxis = dict(
                        nticks=0,showbackground=False,showticklabels=False,visible=False),
                      yaxis = dict(
                        nticks=0,showbackground=False,showticklabels=False,visible=False),
                      zaxis = dict(
                        nticks=0,showbackground=False,showticklabels=False,visible=False),),
                      width=600,
                      font_size=20,
                      title_text=f"{struct.formula}: {struct.get_space_group_info()}",
                    
                  )
    

    
    fig.show()
    
def get_qe_outputs_ionic_relax(file,lattice_init):
    """
    Extract outputs (energies, forces, structures) from qe .stdout files
    
    inputs:
        file: path to the file we want to extract outputs from
    outputs:
        dict: dictionary of the extracted outputs
    """
    
    # TODO: this is very VERY hardcoded. you can do better, fix this (past kat to future kat)
    import numpy as np
    
    output = open(file, "r")
    lines = output.readlines()
    lattice = lattice_init
    iE = [] # energy at each ionic step, Ry
    eE = [[]] # energy at each electronic step, Ry
    P = [] # pressure, kbar
    F = [] # total force, Ry/au
    stresses = [] # stress tensor, kbar
    structures = [] # pymatgen structure objects, angstrom

    from pymatgen.core import Lattice, Structure

    # Check for certain tags on lines, add variables to lists
    for i,line in enumerate(lines):
        if 'total energy' in line and '!' not in line and 'The' not in line:
            eE[-1].append(float(line.split()[3]))
        elif '!' in line:
            eE.append([])
            iE.append(float(line.split()[4]))
        elif 'P=' in line:
            P.append(float(line.split()[5]))
            stresses.append(np.array([lines[i+1].split()[3:6],lines[i+2].split()[3:6],lines[i+3].split()[3:6]]).astype(float))
        elif "Total force" in line:
            F.append(float(line.split()[3]))
        # TODO: come back and fix this, make it more robust 
        # figure out why QE only sometimes gives cell outputs
        elif 'ATOMIC_POSITIONS' in line:
            try:
                sites = []
                atoms = []
                j=1
                while ("End" not in lines[i+j].strip()) and (lines[i+j].strip()!=""):
                    sites.append(np.array(lines[i+j].split()[1:]).astype(float))
                    atoms.append(lines[i+j].split()[0])
                    j=j+1
                #lattice_obj = Lattice(lattice)
                test_struct = Structure(lattice,atoms,sites)
                #print(test_struct)
                structures.append(test_struct)
            except:
                pass
    eE = eE[:-1]

    # return output dictionary
    return {'ionic_energies':iE,'electronic_energies':eE,'pressure':P,'forces':F,'stresses':stresses,'structures':structures}

    
def get_qe_outputs_vc_relax(file):
    """
    Extract outputs (energies, forces, structures) from qe .stdout files
    
    inputs:
        file: path to the file we want to extract outputs from
    outputs:
        dict: dictionary of the extracted outputs
    """
    
    # TODO: this is very VERY hardcoded. you can do better, fix this (past kat to future kat)
    import numpy as np
    
    output = open(file, "r")
    lines = output.readlines()
    iE = [] # energy at each ionic step, Ry
    eE = [[]] # energy at each electronic step, Ry
    P = [] # pressure, kbar
    F = [] # total force, Ry/au
    stresses = [] # stress tensor, kbar
    structures = [] # pymatgen structure objects, angstrom

    from pymatgen.core import Lattice, Structure

    # Check for certain tags on lines, add variables to lists
    for i,line in enumerate(lines):
        if 'total energy' in line and '!' not in line and 'The' not in line:
            eE[-1].append(float(line.split()[3]))
        elif '!' in line:
            eE.append([])
            iE.append(float(line.split()[4]))
        elif 'P=' in line:
            P.append(float(line.split()[5]))
            stresses.append(np.array([lines[i+1].split()[3:6],lines[i+2].split()[3:6],lines[i+3].split()[3:6]]).astype(float))
        elif "Total force" in line:
            F.append(float(line.split()[3]))
        # TODO: come back and fix this, make it more robust 
        # figure out why QE only sometimes gives cell outputs
        elif 'CELL_PARAMETERS' in line:
            try:
                if 'alat' in line:
                    scale = float(line.split()[-1].split(')')[0])*0.529177
                else:
                    scale = 1.0
                lattice = scale*np.array([lines[i+1].split(),lines[i+2].split(),lines[i+3].split()]).astype(float)
                sites = []
                atoms = []
                j=6
                while ("End" not in lines[i+j].strip()) and (lines[i+j].strip()!=""):
                    sites.append(np.array(lines[i+j].split()[1:]).astype(float))
                    atoms.append(lines[i+j].split()[0])
                    j=j+1
                lattice_obj = Lattice(lattice)
                print(lattice, atoms, sites)
                test_struct = Structure(lattice,atoms,sites)
                structures.append(test_struct)
            except:
                pass
    eE = eE[:-1]

    # return output dictionary
    return {'ionic_energies':iE,'electronic_energies':eE,'pressure':P,'forces':F,'stresses':stresses,'structures':structures}

def temp_en_plot(step_dict):
    import matplotlib.pyplot as plt
    import numpy as np
    
    i_energies = step_dict['ionic_energies']
    e_energies_array = step_dict['electronic_energies']
    e_count = [len(e) for e in e_energies_array]
    i_steps = [sum(e_count[0:n+1]) for n in range(len(e_count))]
    e_energies = [item for sublist in e_energies_array for item in sublist]
    e_steps = np.linspace(1,len(e_energies),len(e_energies))

    plt.plot(e_steps, e_energies, label = 'electronic')
    plt.plot(i_steps, i_energies, label = 'ionic')
    
    scaling_factor = 1.005
    plt.xlabel('electronic steps',fontsize=16)
    plt.ylabel('energy (Ry)',fontsize=16)
    plt.ylim(min(i_energies)*scaling_factor,max(i_energies)/scaling_factor)
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig('energy.png',dpi=333)
    
def temp_press_plot(step_dict):
    import matplotlib.pyplot as plt
    import numpy as np
    
    pressures = step_dict['pressure']
    steps = np.linspace(0,len(pressures))

    plt.plot(steps, pressures, label = 'pressure')

    
    scaling_factor = 1.005
    plt.xlabel('ionic steps',fontsize=16)
    plt.ylabel('pressure (kbar)',fontsize=16)
    plt.ylim(min(pressures)*scaling_factor,max(pressures)/scaling_factor)
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig('pressure.png',dpi=333)

def get_en_convergence_plots(step_dict, sim_name = ""):
    """
    Plot both ionic and electronic energy at each SCF step
    
    inputs:
        step_dict: dictionary output from get_qe_outputs()
        sim_name: optional name to add to title on plot
    outputs:
        n/a
    """
    
    import plotly.graph_objects as go
    import numpy as np
    import os
    # Extract the energies, lining up electronic and ionic steps
    i_energies = step_dict['ionic_energies']
    e_energies_array = step_dict['electronic_energies']
    e_count = [len(e) for e in e_energies_array]
    i_steps = [sum(e_count[0:n+1]) for n in range(len(e_count))]
    e_energies = [item for sublist in e_energies_array for item in sublist]
    e_steps = np.linspace(1,len(e_energies),len(e_energies))

    template='simple_white'
    # Create and save a plotly figure with the energy at each ionic and electronic step
    fig_energies = go.Figure()
    fig_energies.add_trace(go.Scatter(x = e_steps, y = e_energies, name = 'electronic'))
    fig_energies.add_trace(go.Scatter(x = i_steps, y = i_energies, name = 'ionic'))

    scaling_factor = 1.005
    fig_energies.update_layout(
        title = f'{sim_name} energy convergence',
        xaxis_title = 'electronic steps',
        yaxis_title = 'energy (Ry)',
        yaxis_range = [min(i_energies)*scaling_factor,max(i_energies)/scaling_factor],
        template=template
    )
    
    fig_energies.show()
    fig_energies.write_image("energy.png")
    
    return True

def get_press_convergence_plots(step_dict, sim_name = ""):
    """
    Plot external pressure at each SCF step
    
    inputs:
        step_dict: dictionary output from get_qe_outputs()
        sim_name: optional name to add to title on plot
    outputs:
        n/a
    """
    
    import plotly.graph_objects as go
    import numpy as np
    import os
    # Extract the pressures
    pressures = step_dict['pressure']
    steps = np.linspace(0,len(pressures))
    template='simple_white'
    
    fig_pressures = go.Figure()
    fig_pressures.add_trace(go.Scatter(x = steps, y = pressures, name = 'external pressure'))

    scaling_factor = 1.005
    fig_pressures.update_layout(
        title = f'{sim_name} pressure convergence',
        xaxis_title = 'ionic steps',
        yaxis_title = 'external pressure (kB)',
        yaxis_range = [min(pressures)*scaling_factor,max(pressures)/scaling_factor],
        template=template
    )
    
    fig_pressures.show()
    fig_pressures.write_image("press.png")
    
    return True 

            
def main():
    # Main loop 
    pass

if __name__ == "__main__":
    main()