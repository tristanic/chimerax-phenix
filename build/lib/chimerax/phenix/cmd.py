# vim: set expandtab ts=4 sw=4:

from chimerax.core.commands import CmdDesc, StringArg
from chimerax.atomic import AtomicStructure, AtomicStructuresArg

import os
_base_dir = os.path.dirname(os.path.abspath(__file__))

def phenix_fit_ligand(session, model, volume, ligand_id, position=None, radius=None,
        markers=None, mask_radius=8):
    m = model
    v = volume

def phenix_get_restraints(session, ligand_id):
    import os
    from .phenix_bridge import run_phenix_process
    command = 'phenix.python'
    args = [
        os.path.join(_base_dir, 'phenix_scripts', 'get_cif_restraints.py'),
        ligand_id,
        os.path.abspath(os.curdir)
    ]
    result = run_phenix_process(session, command, args)
    session.logger.info('Restraints written to {}'.format(result[-1]))
phenix_get_restraints_desc = CmdDesc(
    required=[('ligand_id', StringArg)]
)


def sample_count(session, structures=None, api="python"):
    if structures is None:
        structures = session.models.list(type=AtomicStructure)
    module_def = api_modules.get(api.lower(), None)
    if module_def is None:
        from chimerax.core.errors import UserError
        err_string = ('Unrecognised API! Allowed choices are: {}').format(
            ', '.join(api_modules.keys())
        )
    c_module, attribute = module_def
    for m in structures:
        if attribute is None:
            arg = m
        else:
            arg = getattr(m, attribute)
        atoms, bonds = c_module.counts(arg)
        session.logger.info("%s: %d atoms, %d bonds" % (m, atoms, bonds))
sample_count_desc = CmdDesc(
    optional=[("structures", AtomicStructuresArg),],
    keyword=[("api", StringArg),]
)
