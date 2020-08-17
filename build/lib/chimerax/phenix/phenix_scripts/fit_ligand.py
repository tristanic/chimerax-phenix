
from iotbx.data_manager import DataManager
from iotbx.map_model_manager import map_model_manager
from phenix.model_building import local_model_building
#from phenix.model_building.high_level_tools import run_local_model_building



# Set up a data_manager to read and write files
dm=DataManager()

# Set up a cut-out map containing density for a ligand and part of a model
#  Normally this will be your own map and existing model
map_manager=dm.get_real_map('map_to_fit.ccp4')
model=dm.get_model('model_noligand.cif')
mmm=map_model_manager(map_manager=map_manager, wrapping=False,model=model)
# mam.box_all_maps_around_model_and_shift_origin()
# mam.write_map('boxed_map.ccp4')
#
# # Read in a cut_out map and a model and find ligand in unmodelled density
#
# # Read in the map and model
# map_manager=dm.get_real_map('boxed_map.ccp4')
# model=dm.get_model('model_noligand.pdb')
#
# # make a map-model manager and remove the model outside the map
# mmm=map_model_manager(map_manager=map_manager, model=model)
mmm.remove_model_outside_map(boundary=1.5)

# Get a ligand and restraints (cif) for the ligand (if necessary)
atp_model=dm.get_model('atp.pdb')
dm.process_restraint_file('atp.cif')
atp_restraints=dm.get_restraint('atp.cif')

# Ready with small map and model that is inside the map region and a ligand to fit

# Set up local model building
rlmb=local_model_building(
   map_model_manager=mmm, # map_model manager
   resolution=3,        # d_min
   is_xray_map=False,   # for cryo-em map
   nproc=4,   # more processors mean more tries also
  )

# set any defaults
rlmb.set_defaults(
   scattering_table='electron',  # for a cryo-em map
   thoroughness='medium',  # quick/medium/thorough/extra_thorough
   )

# run ligand fitting (NOTE: you can run again and get a different answer each time)
fitted_ligand_model=rlmb.fit_ligand(
   ligand_model=atp_model,           #ligand model object
   restraints_object=atp_restraints, # optional restraints for ligand
   good_enough_score=0.75,           # stop looking if this is achieved
  )

# write out the ligand
dm.write_model_file(fitted_ligand_model,'fitted_ligand.pdb',overwrite=True)
