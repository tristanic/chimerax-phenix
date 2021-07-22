from http.server import BaseHTTPRequestHandler
import json

def cleanup_temp_files():
    import os, shutil
    for filename in os.listdir(os.curdir):
        if os.path.isfile(filename) or os.path.islink(filename):
            os.unlink(filename)
        elif os.path.isdir(filename):
            shutil.rmtree(filename)

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

class ServerMethods:

    @staticmethod
    def phenix_version_info():
        from phenix import phenix_info
        version, tag = phenix_info.version_and_release_tag()
        return {
            'version':  version,
            'tag':      tag
        }

    @staticmethod
    def find_cif_file(code, working_dir):
        import os, shutil
        local_filename = "{}_restraints.cif".format(code)
        already_present=False
        lf = os.path.join(working_dir, local_filename)
        if os.path.isfile(lf):
            shutil.copy(lf, os.curdir)
            already_present=True
        if not already_present:
            from elbow.utilities import geostd_utils, mmtbx_utils
            filename = None
            filename = geostd_utils.get_geostd_cif_file(code)
            if filename is None:
                filename = mmtbx_utils.get_monomer_cif_file(code)
            if filename is not None:
                # shutil.copyfile(filename, os.path.join(os.curdir, local_filename))
                shutil.copyfile(filename, lf)
            else:
                run_elbow(code)
                shutil.copy(os.path.join(os.curdir, local_filename), lf)
            cleanup_temp_files()
        return {
            'filename': lf,
        }

    @staticmethod
    def working_dir():
        import os
        return {'dir': os.path.abspath(os.curdir)}

    @staticmethod
    def fit_ligand(working_dir, map_file, model_file, ligand_restraint_file, ligand_coord_file,
            resolution=3.0, is_xray_map=False, nproc=4, thoroughness='medium'):
        from iotbx.data_manager import DataManager
        from iotbx.map_model_manager import map_model_manager
        from phenix.model_building import local_model_building
        import os

        dm = DataManager()
        map_manager = dm.get_real_map(map_file)
        model = dm.get_model(model_file)
        mam = map_model_manager(map_manager=map_manager, wrapping=False, model=model)

        mam.remove_model_outside_map(boundary=1.5)
        ligand_model = dm.get_model(ligand_coord_file)
        dm.process_restraint_file(ligand_restraint_file)
        ligand_restraints=dm.get_restraint(ligand_restraint_file)

        rlmb = local_model_building(
            map_model_manager=mam,
           nproc=nproc,
        )

        if is_xray_map:
            scattering_table = 'n_gaussian'
        else:
            scattering_table = 'electron'

        rlmb.set_defaults(
            scattering_table=scattering_table,
            thoroughness=thoroughness
        )

        fitted_ligand_model=rlmb.fit_ligand(
            ligand_model=ligand_model,
            restraints_object=ligand_restraints,
            good_enough_score=0.75,
        )

        fitted_ligand_filename = os.path.join(working_dir, 'fitted_ligand.pdb')
        dm.write_model_file(fitted_ligand_model, fitted_ligand_filename, overwrite=True)
        cleanup_temp_files()
        return {
            'ligand_file': fitted_ligand_filename,
        }

default_server_methods = {
    'version_info': ServerMethods.phenix_version_info,
    'find_cif_file': ServerMethods.find_cif_file,
    'working_dir': ServerMethods.working_dir,
    'fit_ligand': ServerMethods.fit_ligand,

}

class PhenixRESTServer:
    '''
    Listen for HTTP/REST requests, and return machine-friendly
    JSON descriptors of results.
    '''

    def __init__(self, *args, **kw):
        self.httpd = None
        self._server_methods = {}
        self.standard_functions = default_server_methods

        for fname, func in self.standard_functions.items():
            self.register_server_method(fname, func)
        self._server_methods['batch'] = self.batch_run



    def batch_run(self, batch_commands):
        '''
        Run a series of commands in a single call, to avoid communication delays.

        Args:

            batch_commands: a list, where each entry in the list is:
                ['func_name', [args], {kwargs}]
        '''
        ret_dict = {}
        for i, cmd in enumerate(batch_commands):
            fname, args, kwargs = cmd
            f_dict = ret_dict[str(i) + ': ' + fname] = {}
            func = self.standard_functions.get(fname, None)
            if func is None:
                err_msg = 'Unrecognised function name: {}'.format(fname)
                ret_dict['error'] = err_msg
                break
            f_dict.update(func(*args, **kwargs))
        return ret_dict

    @property
    def server_address(self):
        if self.httpd is not None:
            return self.httpd.server_address
        return None

    @property
    def port(self):
        return self.server_address[1]

    def run(self, port):
        import os
        from http.server import HTTPServer
        import sys
        if port is None:
            # Defaults to any available port
            port = 0

        httpd = self.httpd = HTTPServer(('localhost', port), RESTHandler)
        # httpd.chimerax_session = self.session
        httpd.manager = self
        msg = "REST server started on host {} port {}".format(
            *httpd.server_address
        )
        print(msg)
        httpd.serve_forever()

    def terminate(self):
        if self.httpd is not None:
            self.httpd.shutdown()
            self.httpd = None
        super().terminate()

    def register_server_method(self, func_name, func):
        '''
        Register a method to be available as a remote command. The method should
        take only JSON-serialisable types (e.g. (lists of) strings or
        numbers).

        The return type should be a JSON-serialisable dict.
        '''
        if func_name not in self._server_methods:
            self._server_methods[func_name] = func

    @property
    def server_methods(self):
        return self._server_methods

    def list_server_methods(self):
        import inspect
        from collections import defaultdict
        ret_dict = defaultdict(lambda: dict())
        for func_name, func in self.server_methods.items():
            func_dict = ret_dict[func_name]
            func_dict['docstring'] = inspect.getdoc(func)
            arg_list = func_dict['args'] = list()
            kwarg_dict = func_dict['kwargs'] = dict()
            args, varargs, kwargs, defaults = inspect.getargspec(func)
            if defaults is None:
                num_defaults = 0
            else:
                num_defaults = len(defaults)
            if num_defaults:
                for i, default_val in enumerate(reversed(defaults)):
                    arg_name = args[-i-1]
                    arg_props = kwarg_dict[arg_name] = dict()
                    if isinstance(default_val, str):
                        default_val='"{}"'.format(default_val)
                    arg_props['default'] = default_val
                    arg_props['type'] = 'unspecified'
            for arg in args[:len(args)-num_defaults]:
                if arg != 'self':
                    arg_list.append((arg, 'unspecified'))

        return ret_dict



class RESTHandler(BaseHTTPRequestHandler):
    '''Process one REST request.'''

    def _set_headers(self, status_code=200):
        self.send_response(status_code)
        self.send_header('Content-type', 'application/json')
        self.end_headers()

    def do_HEAD(self):
        self._set_headers()

    def do_GET(self):
        '''
        Return a JSON dict listing all available methods
        '''
        self._list_methods()

    def _list_methods(self):
        mgr = self.server.manager
        self._set_headers()
        msg = json.dumps(mgr.list_server_methods())
        self.wfile.write(msg.encode('utf-8'))

    def do_POST(self):
        from cgi import parse_header
        ctype, pdict = parse_header(self.headers.get('content-type'))

        # refuse to receive non-json requests
        if ctype != 'application/json':
            self.send_response(400)
            self.end_headers()
            return

        l = int(self.headers.get('content-length'))
        request = json.loads(self.rfile.read(l).decode('utf-8'))
        return_dict = {}
        try:
            return_dict = self._run_post_job(request)
        except Exception as e:
            import traceback
            err_dict = {'error': str(e),
                'traceback': traceback.format_exc()
                }
            self._set_headers(400)
            self.wfile.write(json.dumps(err_dict).encode('utf-8'))
            return
        self._set_headers()
        self.wfile.write(json.dumps(return_dict).encode('utf-8'))

    def _run_post_job(self, request_dict):
        try:
            func_name = request_dict['cmd']
        except KeyError:
            err_dict = {'error': 'You must provide a command name with the key "cmd"!'}
            return err_dict
        mgr = self.server.manager
        f = mgr.server_methods.get(func_name, None)
        if f is None:
            err_dict = {'error': 'No registered server method with the name {}'.format(func_name)}
            return err_dict
        args = request_dict.get('args', [])
        # Required because json.loads() turns strings to UTF-8 unicode, but
        # various Phenix methods can only handle ASCII. Should be able to dump
        # this once Phenix migrates to Python 3.
        for i, arg in enumerate(args):
            if isinstance(arg, unicode):
                args[i] = arg.decode('utf-8').encode('ascii')
        kwargs = request_dict.get('kwargs', {})
        for kwarg, val in kwargs.items():
            if isinstance(val, unicode):
                kwargs[kwarg] = val.decode('utf-8').encode('ascii')
        return f(*args, **kwargs)

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
        td = self._temp_dir = tempfile.mkdtemp()
        os.chdir(td)

    def __exit__(self, exc_type, exc_value, exc_traceback):
        import os
        os.chdir(self.cwd)
        import shutil
        shutil.rmtree(self._temp_dir)



if __name__ == '__main__':
    server = PhenixRESTServer()
    import sys

    with SafeTempDir():
        server.run(int(sys.argv[1]))
