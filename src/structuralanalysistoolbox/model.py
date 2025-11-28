
from __future__ import annotations
import ansys.mapdl.core as core
from ansys.mapdl.core import launch_mapdl
from dataclasses import dataclass
from typing import Literal
from pathlib import Path
import numpy as np
import time

from structuralanalysistoolbox.mapdl import files
from structuralanalysistoolbox.materials import material 
from structuralanalysistoolbox.materials.material import Material
from structuralanalysistoolbox.constraints import constraint
from structuralanalysistoolbox.mapdl import command, element
from structuralanalysistoolbox.loading import loadstep

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from structuralanalysistoolbox.mapdl.mesh import Mesh
    

@dataclass
class set:
    name : str
    items : np.array

@dataclass
class Nset(set):
    name: str
    items : np.array
    midx = 0

@dataclass
class Elset(set):
    name : str
    items : np.array
    midx = 0
    properties = None

    def __str__(self):
        return (
            f"Id: {self.midx}\n" 
            f"Elements: {self.items}")

@dataclass
class Surface(Nset):
    name : str
    items : np.array

@dataclass
class LocalCoordinateSystem:
    midx : int

class Model:

    def __init__(self, name = '', cpus = 4, interactive = True, work_folder = '', solver : Literal["MAPDL"] = "MAPDL"):
        
        self.solver_name = solver
        self.modelname = name
        self.cpus = cpus
        self.interactive = interactive
        self.work_folder = work_folder

        self._command_map = []

        self._model_definitions = {"Element Types" : {},
                                   "Materials" : {},
                                   "Sets" : {"Node Sets" : {}, 
                                             "Element Sets": {}},
                                   "Surfaces" : {},
                                   "Sections" : {},
                                   "Constraints" : {"Linear Couplings" : {},
                                                    "Constraint Equations" : {},
                                                    "MPC" : {},
                                                    "Joints" : {}},
                                   "Contacts" : {},
                                   "Load Steps" : {"Parameters" : {},
                                                   "Boundary Conditions" : {},
                                                   "Loads" : {},
                                                   "Outputs" : {},
                                                   "Restart" : {}}
                                   } 

        if work_folder == '':
            _current_path = Path.cwd()
            self.creation_time = time.time()
            self.working_directory = _current_path / f'{name}_{self.creation_time}'
            self.working_directory.mkdir()
        
        if interactive:
            self._start()

    def nset(self, name : str):
        return self._model_definitions["Sets"]["Node Sets"][name]

    def elset(self, name : str):
        return self._model_definitions["Sets"]["Element Sets"][name]

    def surface(self, name : str):
        return self._model_definitions["Sets"]["Surfaces"][name]

    def import_mesh(self, file : str):
        """
        # Reads NBlock, EBlock, and Element Attributes (ET) from ansys block structured file format.
        # If there is already an element attribute defined with the same id in the model, 
          it is overriden.
        # In case of interactive modeling, the imported mesh file is not executed by solver. It is
          converted to PyVista Unstructuredmesh format. So, at the pre-processing time, this mesh
          format re-converted to mapdl mesh format again. To do that, an mapdl command mapping was added
          to the end of this method.

        >> Reads:
            # EBLOCK
            # NBLOCK
            # ELEMENT ATTRIBUTES:
                ET, R, MP, SECTYPE, SECDATA 
            # COMPONENTS    

        >> Does not read:
            Loads, BCs, Load Steps, constraints (CP, CE, MPC) etc. 
        """
        self.mesh : Mesh = files._ansys_mesh(file)

        for name, nlist in self.mesh.nset_dict.items():
            set = Nset(name, nlist)
            self.add_set(set)

        for name, elist in self.mesh.elset_dict.items():
            set = Elset(name, elist)
            self.add_set(set)

        # TODO:
        et_list = self.mesh.element_types() # [(etype, etype_idx) ...]
        for etype in et_list:
            self.add_element_type(etype)

        # Create "Element Types" after import
        #   Define "KeyOpts"
        # Create "Real Constants" after import,

        # Map solver command:
        self._add_command("_execute_mesh_data", self.mesh)
     
    def add_element_type(self, etype):
        """Add element type object to the model"""

        if etype.midx != 0:
            self._model_definitions["Element Types"][etype.name] = etype
        else:
            etype.midx = self._find_max_id(self._model_definitions["Element Types"]) + 1
            self._model_definitions["Element Types"][etype.name] = etype
        
    def export_model(self, fileName : str):
        """Exports file in blocked format"""
        self.mapdl.allsel()
        self.mapdl.cdwrite(option="ALL", fileName=fileName, fmat="BLOCKED")

    def add_set(self, set : set):
        """Add Nset, Elset or Surface object to the model"""
        if isinstance(set, Nset):
            set.midx = self._find_max_id(self._model_definitions["Sets"]["Node Sets"]) + 1
            self._model_definitions["Sets"]["Node Sets"][set.name] = set
        elif isinstance(set, Elset):
            set.midx = self._find_max_id(self._model_definitions["Sets"]["Element Sets"]) + 1
            self._model_definitions["Sets"]["Element Sets"][set.name] = set
        elif isinstance(set, Surface):
            set.midx = self._find_max_id(self._model_definitions["Surfaces"]) + 1
            self._model_definitions["Surfaces"][set.name] = set
            
    def add_material(self, mat : str | Material):

        if isinstance(mat, str):
            mat = material.load(mat_name=mat)
            mat.midx = self._find_max_id(self._model_definitions["Materials"]) + 1
            self._model_definitions["Materials"][mat.name] = mat
        elif isinstance(mat, Material):
            mat.midx = self._find_max_id(self._model_definitions["Materials"]) + 1
            self._model_definitions["Materials"][mat.name] = mat
        else: return

        self._add_command("_material", mat)
  
    def merge_nodes(self, nset : str):
        """
        Merges coincident nodes for nset
        Merge those nodes together via NUMMRG.
        """
        return

    def add_coupling(self, name : str, nset : str, dof : str):
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
    
    def add_MPC_Rigid(self, dependent : str | Nset, independent : str | Nset,
                    mpc: constraint.MPCRigidLink | constraint.MPCRigidBeam = None):                
        
        if isinstance(dependent, str):
            dependent = self.nset(dependent)
        if isinstance(independent, str):
            independent = self.nset(independent)

        mpc_obj = None
        if mpc == None or mpc == constraint.MPCRigidLink:
            rigid_link = constraint.MPCRigidLink(dependent=dependent.items, independent=independent.items)
            rigid_link.midx = self._find_max_id(self._model_definitions["Constraints"]["MPC"]) + 1
            mpc_obj = self._model_definitions["Constraints"]["MPC"][f"MPC-{rigid_link.midx}"] = rigid_link
            
        elif mpc == constraint.MPCRigidBeam:
            mpc = constraint.MPCRigidBeam(dependent=dependent.items, independent=independent.items)
            mpc.midx = self._find_max_id(self._model_definitions["Constraints"]["MPC"]) + 1
            mpc_obj = self._model_definitions["Constraints"]["MPC"][f"MPC-{mpc.midx}"] = mpc
            
        else: return

        # Create a new MPC184 Element and Element Type ID
        # We need this to distinguish MPC object id and MPC Element Type id
        mpc_obj.etype = type_mpc = element.EType(name="MPC184")
        self.add_element_type(type_mpc)

        # Map solver command
        self._add_command("_MPC184", mpc_obj)

    def add_joint(self, joint : constraint.MPCJoint):

        # Map solver command
        self._add_command("_create_joint", joint)
    
    def add_symmetry(self):
        # TODO: symmetric / periodic
        pass   

    def add_interface(self):
        # Creates interface for various use cases.
        pass 

    def add_contact(self, master='', slave=''):
        pass

    def add_load_step(self):
        
        ls = loadstep()

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
        # plot mesh 
        if self.mesh:
            self.mesh.grid.plot(show_edges=True)

    def _get_item_object(self, path : tuple[str]):
        """Returns object for any model item"""
        
        node = self._model_definitions
        for key in path:
            if not isinstance(node, dict) or key not in node:
                return None
            node = node[key]
        return node

    def _find_max_id(self, group : dict) -> int:
        """
        Returns max id of an item in a sub-group
        gets a dict such as (group_name[sub-group-1][sub-group-2]..)
        """
        max_id = 0
        if len(group.items()) > 0:
            for item in group.values():
                if item.midx > max_id: max_id = item.midx
        return max_id

    def _find_item_object(self, tree, target_key):
        "depth-first search to find the key"
        for key, value in tree.items():
            if key == target_key:
                return value
            
            if isinstance(value, dict):
                result = self._find_item_object(value, target_key)
                if result is not None:
                    return result
        return None

    def info(self, item : None | str = None) -> str:
        """Returns detailed information string for model item object."""
        if item == None:
            # List model items as a tree structure
            self._print_tree(self._model_definitions, indent=" ")
        else:
            obj = self._find_item_object(self._model_definitions, item)
            print(obj.__str__())
        
    def _print_tree(self, tree, indent=""):
        if not isinstance(tree, dict):
            return

        last_key = list(tree.keys())[-1] if tree else None

        for key in tree:
            is_last = (key == last_key)
            prefix = "└── " if is_last else "├── "
            print(indent + prefix + key)

            next_indent = indent + ("    " if is_last else "│   ")
            if isinstance(tree[key], dict):
                self._print_tree(tree[key], next_indent)

    def _add_command(self, func_name : str, model_item_object):
        """
        adds command to << self._command_map >> list to be executed at the pre-processing time.
        func_name is corresponding mapdl commands function in the mapdl.py
        model_item_object is any kind of object that has information to convey to 
        mapdl commands function.
        """

        if self.solver_name == "MAPDL":
            solver_func = getattr(command, func_name) # Module must be pre-loaded!
        self._command_map.append((solver_func, model_item_object))

    def _execute_commands(self):

        if self.solver_name == "MAPDL":
            self.mapdl = self._start()

        for command in self._command_map:
            command[0](self.mapdl, command[1])

    def solve(self):
        pass
  
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




def _get_min_node(nset : Nset) -> int:
    """
    In order to find independent node for non-interactive sessions.
    For interactive session we can get the value for coupled dofs.
    """
    return min(nset.items)





    





    