# vim: set expandtab ts=4 sw=4:

from chimerax.core.commands import CmdDesc, StringArg, IntArg, BoolArg
from chimerax.atomic import AtomicStructure, AtomicStructuresArg

import os
_base_dir = os.path.dirname(os.path.abspath(__file__))

global restclient
restclient = None

def phenix_connect(session, address='localhost', port=13010, timeout=300,
        reconnect=False, warn_if_already_connected=True):
    global restclient
    if not reconnect:
        if restclient is not None:
            if warn_if_already_connected:
                session.logger.warning('Phenix REST client is already connected!')
            return
    from .rest_server import PhenixRESTClient
    rc = PhenixRESTClient(address, port, timeout=timeout)
    try:
        rc.connect()
    except ConnectionRefusedError:
        import os
        base_dir = os.path.dirname(os.path.abspath(__file__))
        err_str = ('It looks like the Phenix server is not running. To get it '
            'started, run the following command in a separate terminal window '
            'with the Phenix environment initialised: \n\n'
            'phenix.python {}'.format(
                os.path.join(base_dir, 'rest_server','phenix_side', 'server.py')
            ))
        from chimerax.core.errors import UserError
        raise UserError(err_str)
    session.logger.info('Connected to Phenix server on {}:{}'.format(address, port))
    restclient = rc
phenix_connect_desc=CmdDesc(
    synopsis="Initialise the connection to Phenix command server",
    keyword=[
        ('address', StringArg),
        ('port', IntArg),
        ('timeout', IntArg),
        ('reconnect', BoolArg)
    ]
)

def phenix_version(session):
    phenix_connect(session, warn_if_already_connected=False)
    result = restclient.version_info()
    session.logger.info('Phenix Version: {}\nRelease: {}'.format(result['version'], result['tag']))
phenix_version_desc=CmdDesc(
    synopsis="Report the Phenix version to the log."
)


def phenix_fit_ligand(session, model, volume, ligand_id, position=None, radius=None,
        markers=None, mask_radius=8):
    m = model
    v = volume

def phenix_get_restraints(session, ligand_id):
    phenix_connect(session, warn_if_already_connected=False)
    import os
    result = restclient.find_cif_file(ligand_id, os.path.abspath(os.curdir))
    if not result['already_present']:
        session.logger.info('Restraints written to {}'.format(result['filename']))
    else:
        session.logger.info('Restraints for {} already present in working directory'.format(ligand_id))

phenix_get_restraints_desc = CmdDesc(
    required=[('ligand_id', StringArg)]
)
