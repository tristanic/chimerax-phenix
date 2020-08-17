# @Author: Tristan Croll <tic20>
# @Date:   17-Aug-2020
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 17-Aug-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

def run_phenix_process(session, command, args,
        required_phenix_version = None):
    if required_phenix_version is None:
        from . import DEFAULT_MIN_PHENIX_VERSION
        required_phenix_version = DEFAULT_MIN_PHENIX_VERSION
    from . import check_for_needed_phenix_version
    check_for_needed_phenix_version(required_phenix_version)
    with SafeTempDir() as _:
        import sys, os, subprocess
        call_args = [command, *args]
        pipes = subprocess.Popen(call_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        std_out, std_err = pipes.communicate()
        cmd_str = command + ' ' + ' '.join(args)
        session.logger.info('Running Phenix command: \n{}'.format(cmd_str))
        # session.logger.info(std_out.strip().decode('utf-8'))
        if pipes.returncode !=0:
            from chimerax.core.errors import UserError
            err_str = ('Attempt to run the Phenix command {} failed with the following '
                'error: {}'.format(
                    cmd_str, std_err.strip().decode('utf-8')
                ))
            raise UserError(err_str)
        result = std_out.strip().decode('utf-8').split('\n')
        return result

class SafeTempDir:
    '''
    Create and change to a temporary directory. When called with

    with SafeTempDir():
        do_something()

    ... it is guaranteed to change back to the original working directory and
    delete the temporary directory.
    '''
    def __init__(self):
        import os
        self.cwd = os.path.abspath(os.curdir)

    def __enter__(self):
        import tempfile, os
        self._tdmgr = tempfile.TemporaryDirectory()
        td = self.temp_dir = self._tdmgr.__enter__()
        os.chdir(td)

    def __exit__(self, exc_type, exc_value, exc_traceback):
        import os
        os.chdir(self.cwd)
        self._tdmgr.__exit__(exc_type, exc_value, exc_traceback)
