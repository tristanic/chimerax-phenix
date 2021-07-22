from chimerax.core.tasks import Task

_MIN_PORT=49152
_MAX_PORT=65535

class PhenixServer(Task):

    def __init__(self, session, phenix_base_path, retries=5):
        super().__init__(session)
        self._port = None
        self._retries = retries
        import os
        self._phenix_executable = os.path.join(phenix_base_path, 'phenix.python')
    
    @property
    def port(self):
        return self._port
    
    def run(self):
        for i in range(self._retries):
            try:
                self._start_server()
                break
            except RuntimeError as e:
                self.session.logger.warning(str(e))
                continue
        else:
            self.session.logger.warning("Something went wrong starting the Phenix server. Please report this error.")
    
    def _start_server(self):
        import os, subprocess
        curdir = os.path.abspath(os.path.dirname(__file__))
        from random import randint
        port = self._port = randint(_MIN_PORT,_MAX_PORT)
        cmd_args = [self._phenix_executable, os.path.join(curdir, "phenix_side", "server.py"), str(port)]
        pipes = self._pipes = subprocess.Popen(cmd_args)#, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #std_out, std_err = pipes.communicate()
        # print(std_out.strip())
        # if pipes.returncode != 0:
        #     err_message = std_err.strip()
        #     raise RuntimeError(err_message)

    
