
import ansys.mapdl.core as core
from ansys.mapdl.core import launch_mapdl
from materials import material 
from materials.material import Material
from solver.mapdl import _material


from pathlib import Path
import time
from dataclasses import dataclass
from enum import Enum

class ComponentType(Enum):
    NODE = 1
    ELEMENT = 2

'''@dataclass
class Component:

    type : ComponentType
    element_type : str
    material : MaterialModel
    #section'''


class Model:

    def __init__(self, solver='MAPDL', name = '', cpus = 4, interactive = True, work_folder = ''):
        
        self.solver = solver
        self.modelname = name
        self.cpus = cpus
        self.interactive = interactive
        self.work_folder = work_folder
        _creation_time = time.time()

        self._model_definitions = {} # {'materials' : {'material_name' : [material_obj, mat_id],}}
                                    # {'components' : {'comp_name' : [comp_obj, comp_id],}}
                                    # ...
        # Mesh
        # Components
        # Sections?
        # Element Types
        # Materials
        # Contacts
        # Kinematic Constraints & Couplings
        # Analysis & Solver Parameters
        # Load Steps
        # Loads & BCs
        # Outputs
        # Restart Parameters

        if work_folder == '':
            _current_path = Path.cwd()
            self.working_directory = _current_path / f'{name}_{_creation_time}'
            self.working_directory.mkdir()
        
        if interactive:
            self._start()

    def _start(self):

        self.mapdl = launch_mapdl(run_location=self.working_directory,
                                  jobname=self.modelname,
                                  nproc=self.cpus,
                                  override=True)
        
        # First, reset the MAPDL database.
        self.mapdl.clear()

        return self.mapdl
    
    def stop(self):
        self.mapdl.exit()

    def add_mesh(self, mesh_file_path : str):
        
        self.mapdl.input(mesh_file_path)
        # TODO : create component objects after import

    def components(self):
        pass

    def add_material(self, mat : str | Material):

        if "materials" in self._model_definitions:
            id = len(self._model_definitions["materials"]) + 1
        else: 
            self._model_definitions["materials"] = {}
            id = 1

        if isinstance(mat, str):
            mat = material.load(mat_name=mat)
            self._model_definitions["materials"][mat.name] = [mat, id]
        elif isinstance(mat, Material):
            self._model_definitions["materials"][mat.name] = [mat, id]
        else: return

    def add_contact(self, master='', slave=''):
        pass

    def add_load_step(self):
        pass

    def add_bc(self):
        pass

    def add_load(self):
        pass

    def plot(self):
        # plot mesh & BC
        self.mapdl.eplot()

    def info(self):
        # report model definitions
        pass

    def solve(self):
        pass
    
def _create_sets_from_components(mapdl : core.Mapdl):
    pass

def import_mesh_file(mapdl : core.Mapdl, path):

    mapdl.input(path)

def show_mesh(mapdl : core.Mapdl):

    mapdl.eplot()



####
## TEST
####

def main():
    
    model = Model(name="test_model")
    model.add_material(mat="My Steel")
    model.add_material(mat="My Steel-2")

    for mat_data in model._model_definitions["materials"].values():
        _material(mapdl=model.mapdl, mat=mat_data[0], mat_id=mat_data[1])
    


if __name__ == '__main__':
    main()






    