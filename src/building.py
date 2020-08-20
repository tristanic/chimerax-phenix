# @Author: Tristan Croll <tic20>
# @Date:   20-Aug-2020
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 20-Aug-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

def suggest_new_residue_number_for_ligand(model, chain_id):
    '''
    Suggest a suitable residue number for a new ligand, based on what is already
    in the chain.
    '''
    from chimerax.atomic import Residue
    residues = model.residues[model.residues.chain_ids == chain_id]
    ligand_residues = residues[residues.polymer_types==Residue.PT_NONE]
    if ligand_residues:
        return max(ligand_residues.numbers)+1

    last_polymeric_residue_num = max(residues.numbers)
    # Start ligands on a nice round multiple of 1000
    ligand_num = round(last_polymeric_residue_num+1000, -3)
    if ligand_num > 9999:
        raise TypeError('The PDB format does not support residue numbers greater '
            'than 9999. Consider adding your ligand to a different chain.')
    return ligand_num


def new_residue_from_template(model, template, chain_id, center=None,
        residue_number=None, insert_code=' ', b_factor=50, precedes=None):
    '''
    Create a new residue based on a template, and add it to the model.
    '''
    if residue_number is None:
        if chain_id in model.residues.chain_ids:
            residue_number = suggest_new_residue_number_for_ligand(model, chain_id)
        else:
            residue_number = 0
    import numpy
    from chimerax.atomic import Atom
    t_coords = numpy.array([a.coord for a in template.atoms])
    if center is not None:
        t_center = t_coords.mean(axis=0)
        t_coords += numpy.array(center) - t_center
    tatom_to_atom = {}
    r = model.new_residue(template.name, chain_id, residue_number,
        insert=insert_code, precedes=precedes)
    for i, ta in enumerate(template.atoms):
        a = tatom_to_atom[ta] = model.new_atom(ta.name, ta.element)
        a.scene_coord = t_coords[i]
        a.bfactor = b_factor
        r.add_atom(a)
        for tn in ta.neighbors:
            n = tatom_to_atom.get(tn, None)
            if n is not None:
                model.new_bond(a, n)
    set_new_atom_style(model.session, r.atoms)
    return r

def set_new_atom_style(session, atoms):
    from chimerax.atomic import selected_atoms, Atom
    from chimerax.core.commands import run
    atoms.draw_modes = Atom.STICK_STYLE
    residues = atoms.unique_residues
    residues.ribbon_displays=True
    residues.ribbon_hide_backbones=False
    current_sel = selected_atoms(session)
    session.selection.clear()
    atoms.selected = True
    run(session, "color sel bychain; color sel byhetero", log=False)
    session.selection.clear()
    current_sel.selected = True
