"""
for each Configuration in cinfigurations: (181)
	1. copy configuration into Coords.xyz
	2. create run.pbs that will copy Single_point_Energy.dat into
   	   SPE(index).dat
	3. qsub run.pbs
"""
import os
import shutil
import subprocess
import time
from glob import glob

PBS_SEND_COMMAND = "qsub"
DEFAULT_COORDS_FILENAME = "Coords.xyz"
DEFAULT_PBS_FILENAME = "run.pbs"
DEFAULT_PBS_CODE = """#PBS -N Zig_CNT.243.0.241.0.Ters.KC.Opt
#PBS -l nodes=1:ppn=8,walltime=2400:00:00
#PBS -l mem=7000mb
#PBS -q cores-8
#PBS -S /bin/tcsh
#--------------------------------
echo "Node Name: $HOSTNAME"
setenv SCRDIR /scratch/${LOGNAME}_%d/$PBS_JOBID
/bin/mkdir -p $SCRDIR
cd $SCRDIR
cp $PBS_O_WORKDIR/27.5.14_MD.x $SCRDIR
cp $PBS_O_WORKDIR/CH.airebo $SCRDIR
cp $PBS_O_WORKDIR/run.par $SCRDIR
cp $PBS_O_WORKDIR/Coords.xyz $SCRDIR
./27.5.14_MD.x
cp -f $SCRDIR/* $PBS_O_WORKDIR
cp -f $SCRDIR/Single_point_Energy.dat $PBS_O_WORKDIR/%s
rm -rf $SCRDIR
"""


def send_config_to_pbs(config_file, index, code_dir):
	# copy configuration into Coords.xyz
	copy_config_to_coords(config_file, code_dir)
	# create related run.pbs file
	create_run_pbs(index, code_dir)
	# run the pbs file
	run_pbs(code_dir)

def copy_config_to_coords(config_file, code_dir):
	dst_path = os.path.join(code_dir, DEFAULT_COORDS_FILENAME)
	shutil.copy(config_file, dst_path)

def create_run_pbs(index, code_dir):
	output_file = "SPE.%d.dat" % index
	pbs_path = os.path.join(code_dir, DEFAULT_PBS_FILENAME)
	pbs_code = DEFAULT_PBS_CODE % (index, output_file)
	open(pbs_path, 'wb').write(pbs_code)

def run_pbs(code_dir):
	os.chdir(code_dir)
	pbs_file = os.path.join(code_dir, DEFAULT_PBS_FILENAME)
	process = subprocess.Popen([PBS_SEND_COMMAND, pbs_file],
		stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = process.communicate()
	print out
	print err
	return not err

def main():
	code_path = "/home/inbaloz/Projects/NT_on_surfaces/interlayer_potential_Itai"
	all_configs = glob("/home/inbaloz/Projects/configs_RI/sliding_paper*atoms*")
	all_configs.sort(key=lambda x: int(x.split(" ")[-1]))
	print all_configs
	for idx, config_file in enumerate(all_configs):
		if 100 < idx <= 120:
			time.sleep(5)
			send_config_to_pbs(config_file, idx, code_path)

if __name__ == '__main__':
	main()

