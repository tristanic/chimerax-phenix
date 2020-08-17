# @Author: Tristan Croll <tic20>
# @Date:   17-Aug-2020
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 17-Aug-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

def find_cif_file(code, working_dir):
    import os, shutil
    local_filename = "{}_restraints.cif".format(code)
    lf = os.path.join(working_dir, local_filename)
    if os.path.isfile(lf):
        shutil.copy(lf, os.curdir)
        return local_filename
    from elbow.utilities import geostd_utils, mmtbx_utils
    filename = None
    filename = geostd_utils.get_geostd_cif_file(code)
    if filename is None:
        filename = mmtbx_utils.get_monomer_cif_file(code)
    if filename is not None:
        shutil.copyfile(filename, os.path.join(os.curdir, local_filename))
        shutil.copyfile(filename, lf)
        return local_filename
    run_elbow(code)
    shutil.copy(os.path.join(os.curdir, local_filename), lf)
    return local_filename

def run_elbow(code):
    import subprocess, shutil
    cmd_args = ['phenix.elbow', '--chemical_component', code]
    pipes = subprocess.Popen(cmd_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    std_out, std_err = pipes.communicate()

    if pipes.returncode != 0:
        err_str = ("Tried to generate restraints for ligand {} with phenix.elbow, "
            "but it failed with the below error message. Please generate restraints "
            "externally, and save in your working directory as {}. \n\n"
            "Error message:\n{}").format(
                code, code+'_restraints.cif', std_err.strip().decode('utf-8')
            )
        raise RuntimeError(err_str)
    shutil.move(code+'.cif', code+'_restraints.cif')

if __name__=='__main__':
    import sys
    print(find_cif_file(*sys.argv[1:]))
