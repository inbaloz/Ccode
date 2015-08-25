import os
import sys

def convert_coords_to_traj(path):
	atoms = open(path,'rb').read().splitlines()[1:]
	new_traj_text = "%d\n0.0000000000000000\n" % len(atoms)
	new_traj_text += "\n".join("C " + " ".join(atom.split()[2:]) for atom in atoms)
	open(path+".xyz", 'wb').write(new_traj_text)

def main():
	convert_coords_to_traj(sys.argv[1])

if __name__ == '__main__':
	main()
