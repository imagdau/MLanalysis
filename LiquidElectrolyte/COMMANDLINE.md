# Command line options

## Corp trajectory

`ase convert -n 0:10000:10 traj_2.1.xyz traj_2.1_0-10000-10.xyz`

## Computing SOAP

`quip atoms_filename="traj_2.1_0-1000.xyz" descriptor_str="soap cutoff=6.0 cutoff_transition_width=1.0  n_max=8 l_max=4 atom_sigma=0.5  n_Z=1 Z={3} n_species=6 species_Z={1 3 6 8 9 15}" output_file="in.out"`
