# vim: set expandtab ts=4 sw=4:

from chimerax.core.commands import (
    CmdDesc, StringArg, IntArg, BoolArg, FloatArg, Float3Arg,
)
from chimerax.atomic import (
    AtomicStructure, AtomicStructuresArg,
    AtomsArg, ResiduesArg,
)
from chimerax.map import MapArg
from chimerax.core.errors import UserError

DEFAULT_PORT=15010
from . import DEFAULT_MIN_PHENIX_VERSION

import os
_base_dir = os.path.dirname(os.path.abspath(__file__))

global restclient
restclient = None

global _phenix_version
_phenix_version = 0

def phenix_connect(session, address='localhost', port=DEFAULT_PORT, timeout=300,
        reconnect=False, warn_if_already_connected=True):
    global restclient, _phenix_version
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
        raise UserError(err_str)
    session.logger.info('Connected to Phenix server on {}:{}'.format(address, port))
    restclient = rc
    version_info = rc.version_info()
    release = version_info['tag']
    if release.isnumeric():
        _phenix_version = int(release)
    else:
        from .. import DEV_PHENIX_VERSION
        _phenix_version = DEV_PHENIX_VERSION

def check_for_needed_phenix_version(session, version):
    phenix_connect(session, warn_if_already_connected=False)
    from chimerax.core.errors import UserError
    if _phenix_version < version:
        err_str = ('This method requires Phenix version {}, but you have version '
            '{} installed. Please update your Phenix installation.').format(
                version, installed_version
            )
        raise UserError(err_str)


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


def phenix_fit_ligand(session, model, ligand_id, chain_id, volume, ligand_residue=None,
        position=None, radius=20, resolution=3.0):
    import numpy, os
    from chimerax.save_command.cmd import provider_save
    from .building import new_residue_from_template
    phenix_connect(session, warn_if_already_connected=False)
    curdir = os.path.abspath(os.curdir)
    server_dir = restclient.working_dir()['dir']
    check_for_needed_phenix_version(session, DEFAULT_MIN_PHENIX_VERSION)
    if len(model) != 1:
        raise UserError('"model" must specify a single atomic structure!')
    model = model[0]
    cleanup_ligand_model = False
    if ligand_residue is None:
        cleanup_ligand_model = True
        from chimerax.atomic import AtomicStructure, mmcif
        try:
            tmpl = mmcif.find_template_residue(session, ligand_id)
        except ValueError:
            raise UserError('If ligandResidue is not provided, ligandId must be in the chemical components dictionary.')
        ligand_model = AtomicStructure(session)
        ligand_residue = new_residue_from_template(ligand_model, tmpl, 'A', [0,0,0])
    else:
        if not len(ligand_residue):
            raise UserError('ligandResidue argument specifies an empty selection!')
        elif len(ligand_residue) != 1:
            raise UserError('ligandResidue argument must specify a single residue')
        ligand_residue = ligand_residue[0]
        if ligand_residue.name != ligand_id:
            raise UserError('Name of ligandResidue must match ligandId!')
        ligand_model = ligand_residue.structure

    ligand_file = os.path.join(server_dir, 'ligand.pdb')
    provider_save(session, ligand_file, models=[ligand_model])
    if cleanup_ligand_model:
        ligand_model.delete()


    phenix_get_restraints(session, ligand_id, log_if_present=False)
    ligand_restraint_file = os.path.join(curdir, '{}_restraints.cif'.format(ligand_id))
    m = model
    v = volume
    if position is None:
        position = session.view.center_of_rotation
    from chimerax.map import data as vd, volume_from_grid_data
    subgrid = vd.zone_masked_grid_data(v.data, [position], radius, minimal_bounds=True)
    temp_volume = volume_from_grid_data(subgrid, session, open_model=False)
    offset = v.scene_position.inverse()*numpy.array(temp_volume.data.origin)
    temp_volume.data.set_origin([0,0,0])
    from chimerax.geometry import Place
    temp_volume.position = Place()
    from chimerax.geometry import find_close_points
    i1, i2 = find_close_points([position], model.atoms.scene_coords, radius)
    final_atoms = model.atoms[i2].unique_residues.atoms
    from chimerax.std_commands.split import molecule_from_atoms
    temp_model = molecule_from_atoms(model, final_atoms)
    from chimerax.atomic import Residue
    temp_model.residues.ss_types=Residue.SS_COIL
    temp_model.atoms.coords -= offset
    temp_model.id=(1,)
    map_file = os.path.join(server_dir, 'map_to_fit.mrc')
    provider_save(session, map_file, models=[temp_volume])
    model_file = os.path.join(server_dir, 'surrounding_model.cif')
    provider_save(session, model_file, models=[temp_model])
    temp_volume.delete()
    temp_model.delete()
    result = restclient.fit_ligand(curdir, map_file, model_file, ligand_restraint_file, ligand_file,
        resolution=resolution)
    from chimerax.open_command.cmd import provider_open
    fitted_ligand_model = provider_open(session, [result['ligand_file']], _add_models=False)[0]
    fitted_ligand_model.atoms.coords += offset
    new_residue_from_template(model, fitted_ligand_model.residues[0], chain_id)
    fitted_ligand_model.delete()

phenix_fit_ligand_desc=CmdDesc(
    synopsis='Fit a ligand into a density blob',
    required=[
        ('model', AtomicStructuresArg),
        ('ligand_id', StringArg),
        ('chain_id', StringArg),
        ('volume', MapArg),
    ],
    keyword=[
        ('ligand_residue', ResiduesArg),
        ('position', Float3Arg),
        ('radius', FloatArg),
        ('resolution', FloatArg),
    ]
)


def phenix_get_restraints(session, ligand_id, log_if_present=True):
    phenix_connect(session, warn_if_already_connected=False)
    import os
    curdir = os.path.abspath(os.curdir)
    if os.path.isfile(os.path.join(curdir, '{}_restraints.cif'.format(ligand_id))):
        if log_if_present:
            session.logger.info('Restraints for {} already present in working directory'.format(ligand_id))
            return
    session.logger.status('Fetching CIF restraints for {}. This may take a while if they need to be generated.'.format(ligand_id))
    result = restclient.find_cif_file(ligand_id, curdir)
    session.logger.info('Restraints written to {}'.format(result['filename']))
    session.logger.status('')

phenix_get_restraints_desc = CmdDesc(
    required=[('ligand_id', StringArg)]
)
