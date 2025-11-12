
import ansys.mapdl.core as core
from ansys.mapdl.core import launch_mapdl
from dataclasses import dataclass
from pathlib import Path
from enum import Enum
import time

from .materials import material 
from .materials.material import Material
from .constraints import constraint


class ComponentType(Enum):
    NODE = 1
    ELEMENT = 2

@dataclass
class Component:
    name : str
    type : str
    items : list

@dataclass
class Nset(Component):
    type = "Node"

@dataclass
class Elset(Component):
    type = "Element"


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
        # Sets
        # Sections?
        # Element Types
        # Materials
        # Contacts
        # Constraints
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

        # TODO : move it to the solver module
        self.mapdl.exit()

    def add_mesh(self, mesh_file_path : str):
        
        # TODO : move it to the solver module
        self.mapdl.input(mesh_file_path)

        # create component objects after import
        for name, type in self.mapdl.components.items():
            if type == "ELEM":
                comp = Elset(name, type, self.mapdl.components[name].items)
                self.add_set(comp)
            elif type == "NODE":
                comp = Nset(name, type, self.mapdl.components[name].items)
                self.add_set(comp)

    def add_set(self, comp : Component):
        
        if "sets" in self._model_definitions:
            id = len(self._model_definitions["sets"]) + 1
        else: 
            self._model_definitions["sets"] = {}
            id = 1

        self._model_definitions["sets"][comp.name] = (comp, id)

    def add_material(self, mat : str | Material):

        if "materials" in self._model_definitions:
            id = len(self._model_definitions["materials"]) + 1
        else: 
            self._model_definitions["materials"] = {}
            id = 1

        if isinstance(mat, str):
            mat = material.load(mat_name=mat)
            self._model_definitions["materials"][mat.name] = (mat, id)
        elif isinstance(mat, Material):
            self._model_definitions["materials"][mat.name] = (mat, id)
        else: return


    def add_coupling(self, nset : str):

        pass

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

    def get(self):
        pass

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


class Node:

    def __init__(self, model : Model):

        self.id : int
        self.dof : int
        self.mapdl = model.mapdl

    def _create_free_node(self):
        """finds max node id, creates a new one with the next number."""
    





    