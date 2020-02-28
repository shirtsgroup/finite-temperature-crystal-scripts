import sys
import os
import subprocess
import numpy as np
import intertools as it
import mdtraj as md
import yaml

def write_tinker_xyz(template, save_path):
    # Setting the first line
    s = template[0, 0] + ' \n'

    # Adding in the lattice parameters
    for i in range(6):
        s += '{:14.8f}'.format(template[1, i])
    s += ' \n'

    # Adding in the coordinate
    for i in range(len(template[2:, 0])):
        s += '{:6d}'.format(template[2 + i, 0]) + '{:4d}'.format(template[2 + i, 1])
        s += '{:14.8f}'.format(template[2 + i, 2]) + '{:14.8f}'.format(template[2 + i, 3]) \
             + '{:14.8f}'.format(template[2 + i, 4])
        for j in template[2 + i, 5:]:
            if j != '':
                s += '{:6d}'.format(j)
        if not i == len(template[2:, 0]) - 1:
            s += ' \n'

    # Save the coordinate file
    with open(save_path, 'w') as file_out:
        file_out.write(s)


def setup_boltz_QHA(inputs):
    # Check to see if boltz_QHA should even be run
    if not inputs['QHA_in']['boltz_QHA']:
        return

    # Load in the template .xyz file
    with open('%s' % inputs['QHA_in']['Tinker_template'], 'r') as l:  # Opening template file
        template = [lines.split() for lines in l]
        template = np.array(list(it.zip_longest(*template, fillvalue=''))).T

    # Cycle through the MD temperatuers and determine, which ones will be used for boltzmann weighting of QHA
    for i, poly in enumerate(inputs['gen_in']['polymorph_num'].split()):
        if not os.path.isdir(poly + 'boltz_QHA'):
            subprocess.call(['mkdir', poly + 'boltz_QHA'])

        for j, T_MD in enumerate(np.array(inputs['temp_in']['temperatures'].split()).astype(float)):
            if T_MD in np.array(inputs['QHA_in']['boltz_temperatures'].split()).astype(float):
                dirpath = poly + '/temperature/' + str(j) + '/'
                if os.path.isfile(dirpath + 'PROD.trr') and not os.path.isdir(poly + 'boltz_QHA/' + str(j)):
                    # Making a directory to store the files in for this temperature
                    subprocess.call(['mkdir', poly + 'boltz_QHA/' + str(j)])

                    # Centering the files
                    c = subprocess.Popen(['echo', '0', '0'], stdout=subprocess.PIPE)
                    output = subprocess.check_output(['gmx', 'trjconv', '-f', dirpath + 'PROD.trr', '-s',
                                                      dirpath + 'PROD_0.tpr', '-o', 'hold.trr', '-center', '-pbc',
                                                      'whole'], stdin=c.stdout)
                    c.wait()

                    # Loading in the trajectory
                    trajectory = md.load('hold.trr', top=dirpath + 'pre_EQ.gro')
                    xyz = trajectory.xyz[::-1][:30] * 10.
                    lattice_vectors = trajectory.unitcell_lengths[::-1][:inputs['QHA_in']['max_frames']] * 10.
                    lattice_angles = trajectory.unitcell_angles[::-1][:inputs['QHA_in']['max_frames']]

                    # Lattice dynamics inputs
                    LD_inputs = {'temperature': T_MD,
                                 'method': 'HA',
                                 'number_of_molecules': inputs['gen_in']['number_of_molecules']}

                    # Run through the frames and create directories
                    for k in range(inputs['QHA_in']['max_frames']):
                        # Create the directory to store the frame in
                        frame_dirpath = poly + 'boltz_QHA/' + str(j) + '/' + str(k) + '/'
                        subprocess.call(['mkdir', frame_dirpath])

                        # Create the .pdb file
                        template[1, :3] = lattice_vectors[k].astype('U32')
                        template[1, 3:6] = lattice_angles[k].astype('U32')
                        template[2:, 2:5] = xyz[k].astype('U32')
                        write_tinker_xyz(template, frame_dirpath + 'molecule.xyz')

                        # Place the other files in the directory
                        for l in inputs['QHA_in']['runfiles']:
                            subprocess.call(['cp', l, frame_dirpath])
                            with open(frame_dirpath + 'input.yaml', 'w') as yaml_file:
                                yaml.dump(LD_inputs, yaml_file, default_flow_style=False)

                    # Removing the extra .trr file
                    subprocess.call(['rm', 'hold.trr'])

