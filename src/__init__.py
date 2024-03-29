# vim: set expandtab shiftwidth=4 softtabstop=4:

from chimerax.core.toolshed import BundleAPI

DEFAULT_MIN_PHENIX_VERSION=3964
DEV_PHENIX_VERSION=int(1e9)

# global _installed_phenix_version
_installed_phenix_version = -1
# try:
#     _installed_phenix_version = check_for_phenix_version()
# except:
#     pass

def _choose_phenix_directory(session):
    satisfied = False
    from Qt.QtWidgets import QFileDialog, QMessageBox
    parent = session.ui.main_window
    import subprocess, os
    while not satisfied:
        result = QFileDialog.getExistingDirectory(parent, 'Please provide the directory containing the Phenix executables.', options=QFileDialog.Options())
        if not result:
            break
        try:
            subprocess.call([os.path.join(result,'phenix.version')])
            satisfied = True
        except FileNotFoundError:
            choice = QMessageBox.warning(parent, 'This directory does not appear to contain Phenix executables. Would you like to try again?',
                QMessageBox.Ok|QMessageBox.Cancel)
        except:
            raise
    if not satisfied:
        from chimerax.core.errors import UserError
        raise UserError('Could not find Phenix installation. Operation cancelled')
    return result
            

        
    

def check_for_phenix_version(session):
    global _installed_phenix_version
    if _installed_phenix_version != -1:
        return _installed_phenix_version
    from chimerax.core.errors import UserError
    import sys, os, subprocess
    from .phenix_bridge import settings
    phenix_path = settings.phenix_base_path
    if phenix_path is None:
        phenix_path = _choose_phenix_directory(session)
        settings.phenix_base_path = phenix_path
    cmd_args = [os.path.join(phenix_path,"phenix.version")]
    try:
        pipes = subprocess.Popen(cmd_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except FileNotFoundError:
        session.logger.warning('Phenix installation appears to have moved or been deleted. Please provide a new path.')
        settings.phenix_base_path = None
        return check_for_phenix_version(session)
    std_out, std_err = pipes.communicate()
    lines = [line.strip() for line in std_out.decode('utf-8').split('\n')]
    version = _parse_version(lines)
    _installed_phenix_version = version
    return version

def _parse_version(lines):
    version = -1
    for line in lines:
        if line.startswith('Release'):
            version = line.split(':')[-1].strip()
            if version.isnumeric():
                version = int(version)
            else:
                version = DEV_PHENIX_VERSION
    return version

def check_for_phenix(session):
    version = check_for_phenix_version(session)
    return (version >= DEFAULT_MIN_PHENIX_VERSION, version)

def check_for_needed_phenix_version(version):
    global _installed_phenix_version
    from chimerax.core.errors import UserError
    installed_version = _installed_phenix_version
    if installed_version < version:
        err_str = ('This method requires Phenix version {}, but you have version '
            '{} installed. Please update your Phenix installation.').format(
                version, installed_version
            )
        raise UserError(err_str)



class _PhenixPluginAPI(BundleAPI):

    api_version = 1     # register_command called with CommandInfo instance
                        # instead of string

    # Override method for starting tool
    # @staticmethod
    # def start_tool(session, bi, ti, **kw):
    #     from .tool import SampleTool
    #     return SampleTool(session, ti.name, **kw)

    # Override method for registering commands
    @staticmethod
    def register_command(bi, ci, logger):
        # We expect that there is a function in "cmd"
        # corresponding to every registered command
        # in "setup.py.in" and that they are named
        # identically (except with '_' replacing spaces)
        from . import cmd
        from chimerax.core.commands import register
        command_name = ci.name
        base_name = command_name.replace(" ", "_")
        func = getattr(cmd, base_name)
        desc = getattr(cmd, base_name + "_desc")
        if desc.synopsis is None:
            desc.synopsis = ci.synopsis
        register(command_name, desc, func)

    # Implement provider method for opening file
    # @staticmethod
    # def run_provider(session, name, mgr):
    #     # 'run_provider' is called by a manager to invoke the
    #     # functionality of the provider.  Since the "data formats"
    #     # manager never calls run_provider (all the info it needs
    #     # is in the Provider tag), we know that only the "open
    #     # command" manager will call this function, and customize
    #     # it accordingly.
    #     #
    #     # The 'name' arg will be the same as the 'name' attribute
    #     # of your Provider tag, and mgr will be the corresponding
    #     # Manager instance
    #     #
    #     # For the "open command" manager, this method must return
    #     # a chimerax.open_command.OpenerInfo subclass instance.
    #     #
    #     # If your bundle also saved or fetched files, then it
    #     # need to distinguish what to do based on the 'name'
    #     # argument or by testing the 'mgr' argument against
    #     # session.open_command or session.fetch_command.  See
    #     # the developer tutorial here:
    #     #   http://www.cgl.ucsf.edu/chimerax/docs/devel/tutorials/introduction.html#writing-bundles-in-seven-easy-steps
    #     # for more info
    #     from chimerax.open_command import OpenerInfo
    #     class XyzOpenerInfo(OpenerInfo):
    #         def open(self, session, data, file_name, **kw):
    #             # The 'open' method is called to open a file,
    #             # and must return a (list of models created,
    #             # status message) tuple.
    #             from .io import open_xyz
    #             return open_xyz(session, data)
    #     return XyzOpenerInfo()

    # Override method for initialization function called each time
    # ChimeraX starts.  Only invoked if the custom initialization
    # flag is set in bundle_info.xml.
    @staticmethod
    def initialize(session, bi):
        # bundle-specific initialization (causes import)
        from . import phenix_bridge
        phenix_bridge.settings = phenix_bridge._PhenixSettings(session, 'phenix plugin')

    # Override method for finalization function.
    # Only invoked if the custom initialization
    # flag is set in bundle_info.xml.
    @staticmethod
    def finish(session, bi):
        # deinitialize bundle in session (causes import)
        raise NotImplementedError

    # Override method to support saving tools in sessions
    # @staticmethod
    # def get_class(class_name):
    #     # 'get_class' is called by session code to get class from bundle that
    #     # was saved in a session
    #     raise NotImplementedError
    #     # "class_name" should be the name of one of the tools
    #     # in this bundle, so code might look something like:
    #     if class_name == 'ToolUI':
    #         from . import tool
    #         return tool.ToolUI
    #     return None


bundle_api = _PhenixPluginAPI()
