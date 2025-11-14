
import ansys.mapdl.core as core
from ansys.mapdl.core import launch_mapdl
from dataclasses import dataclass
from pathlib import Path
from enum import Enum
from typing import Literal
import time
import functools

from .meshing import mesh
from .materials import material 
from .materials.material import Material
from .constraints import constraint
from .constraints.constraint import (MPCtypes, MPCmethods)
from .solver import mapdl


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

@dataclass
class Surface:
    pass


class Model:

    def __init__(self, solver='MAPDL', name = '', cpus = 4, interactive = True, work_folder = ''):
        
        self.solver = solver
        self.modelname = name
        self.cpus = cpus
        self.interactive = interactive
        self.work_folder = work_folder
        _creation_time = time.time()

        self._element_types = []
        self._model_definitions = {} # {'materials' : {'material_name' : [material_obj, mat_id],}}
                                     # {'components' : {'comp_name' : [comp_obj, comp_id],}}
                                    # ...
        # Mesh
        # Element Types !!!
        # Sets
        # Sections?
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
        #self.mapdl.input(mesh_file_path)
        self.mapdl.cdread(option="DB", fname=mesh_file_path)

        self.mapdl.db.
        

        # Create "Node/Element sets" after import
        for name, type in self.mapdl.components.items():
            if type == "ELEM":
                comp = Elset(name, type, self.mapdl.components[name].items)
                self.add_set(comp)
            elif type == "NODE":
                comp = Nset(name, type, self.mapdl.components[name].items)
                self.add_set(comp)

        # Create "Element Types" after import

        # Create "Real Constants" after import


    def add_set(self, comp : Component):
        
        if "sets" in self._model_definitions:
            id = self._find_max_id("sets") + 1
        else: 
            self._model_definitions["sets"] = {}
            id = 1

        self._model_definitions["sets"][comp.name] = (comp, id)

    def add_material(self, mat : str | Material):

        if "materials" in self._model_definitions:
            id = self._find_max_id("materials") + 1
        else: 
            self._model_definitions["materials"] = {}
            id = 1

        if isinstance(mat, str):
            mat = material.load(mat_name=mat)
            self._model_definitions["materials"][mat.name] = (mat, id)
        elif isinstance(mat, Material):
            self._model_definitions["materials"][mat.name] = (mat, id)
        else: return

    
    def merge_nodes(self, nset : str):
        """
        Merges coincident nodes for nset
        Merge those nodes together via NUMMRG.
        """
        return


    def add_dof_coupling(self, name : str, nset : str, dof : str):
        """
        Adds a new dof coupling to the model. 
        dof: "UX/UY/UZ/ROTX/ROTY/ROTZ" or
             "ALL" or
             ("UX", "UY", "ROTX", ..)
        """
        # Get nset object
        nset_obj = self._get_item_object(group="sets", item=nset)
        min_node_num = _get_min_node(nset_obj) # For non-interactive sessions

        # create a dof-coupling obj
        coupling_obj = constraint.CoupledDOF(nset=nset_obj, dof=dof)

        if "couplings" in self._model_definitions:
            id = self._find_max_id("couplings") + 1
        else: 
            self._model_definitions["couplings"] = {}
            id = 1

        coupling_obj.id = id
        self._model_definitions["couplings"][name] = (coupling_obj, id)

        return min_node_num
    
    def add_MPC(self, 
                name: str,
                type: Literal["RigidLink",
                              "RigidBeam",
                              "Slider",
                              "Revolute",
                              "Universal",
                              "Slot",
                              "Point",
                              "Translational",
                              "Cylindrical",
                              "Planar",
                              "Weld",
                              "Orient",
                              "Spherical",
                              "General",
                              "Screw"], 
                independents: Nset,
                dependent: Nset,
                method: None | Literal["DirectElemination",
                                       "LagrangeMultiplication",
                                       "SurfaceBased"] = None):                
        
        # Find MPC Model Item ID
        if "Multi-Point-Constraints" in self._model_definitions:
            id = self._find_max_id("Multi-Point-Constraints") + 1
        else: 
            self._model_definitions["Multi-Point-Constraints"] = {}
            id = 1

        # Find max type id
        max_type_no = max((_type[1] for _type in self._element_types), default=0)
        
        mpc = constraint.MPC(name=name,
                             id=id, # MPC id
                             etype_id=max_type_no+1, # MPC184 Element id
                             type=type, 
                             method=method, 
                             independent=independents, 
                             dependent=dependent)
        
        # Add MPC184 "Element Type" definition to the model list.
        self._element_types.append((mpc._etype, max_type_no+1))

        # Add MPC Object as Model Item
        self._model_definitions["Multi-Point-Constraints"][name] = (mpc, id)

        # Run with mapdl
        mapdl._MPC184(self.mapdl, mpc)
        

    def add_symmetry(self):
        # TODO: symmetric / periodic
        pass   

    def add_interface(self):
        # Creates interface for various use cases.
        pass 

    def add_contact(self, master='', slave=''):
        pass

    def add_load_step(self):
        pass

    def add_bc(self):
        pass

    def add_load(self):
        pass

    def remove(self, group : str, item : str) -> None | str:
        '''Removes items from the model'''
        # TODO: remove items from mapdl for interactive sessions.
        try:
            self._model_definitions[group].pop(item)
        except KeyError:
            raise KeyError(f"Item '{group}-{item}' does not exist")
        return None

    def plot(self):
        # plot mesh & BC
        self.mapdl.eplot()

    def _get_item_object(self, group : str, item : str):
        """Returns object for any model item"""
        try:
            obj = self._model_definitions[group][item][0]
        except KeyError:
            raise KeyError(f"Item '{group}-{item}' does not exist")
        return obj

    def info(self, group : str,  item : str) -> str:
        """Returns detailed information string for model items."""
        pass

    def _find_max_id(self, group_name : str) -> int:
        """Returns max id number in a group"""
        if len(self._model_definitions[group_name]):
            max_id = max(v[1] for v in self._model_definitions[group_name].values())
            return max_id
        else:
            return 0
            

    def __str__(self):
        # List model items as a tree structure
        lines = []
        lines.append(f"  {'Model':22}{'ID'}")
        lines.append("-" * 40)

        for group_name, group in self._model_definitions.items():
            lines.append(f"• {group_name}")
            for item_name, item in group.items():
                lines.append(f"    {item_name:20} {item[1]}")

        return "\n".join(lines)

    def _execute_commands(self):
        pass

    def solve(self):
        pass

    """def solver(self, func):
        @functools.wraps(func)
        def wrapper_solver(*args, **kwargs):  
            val = func(*args, **kwargs)
            args[0].finish()
            return val
        return wrapper_solver"""

    
def _create_sets_from_components(mapdl : core.Mapdl):
    pass

def _get_min_node(nset : Nset) -> int:
    """
    In order to find independent node for non-interactive sessions.
    For interactive session we can get the value for coupled dofs.
    """
    return min(nset.items)

def import_mesh_file(mapdl : core.Mapdl, path):

    mapdl.input(path)

def show_mesh(mapdl : core.Mapdl):

    mapdl.eplot()




    





    