
from __future__ import annotations
import ansys.mapdl.core as core
from ansys.mapdl.core import launch_mapdl
from dataclasses import dataclass
from typing import Literal
from pathlib import Path
import numpy as np
import time
from datetime import datetime
from typing import Dict, List, Any

from structuralanalysistoolbox.data import mesh
from structuralanalysistoolbox.materials import material 
from structuralanalysistoolbox.materials.material import Material
from structuralanalysistoolbox.constraints import constraint, contact
from structuralanalysistoolbox.mapdl import command, element
from structuralanalysistoolbox.loading.loadstep import LoadStep, Force, Moment, Pressure, Displacement
from structuralanalysistoolbox.plotting.plots import LoadStepPlot
from structuralanalysistoolbox.visualization.meshplotter import Plotter
from structuralanalysistoolbox.exceptions import ParameterError

@dataclass
class Nset:
    name : str
    items : np.ndarray | list
    midx = 0
    ignore_at_execution = False

@dataclass
class Elset:
    name : str
    items : np.ndarray | list
    midx = 0
    properties = {"Material" : None}
    ignore_at_execution = False

    def __str__(self):
        return (
            f"Id: {self.midx}\n" 
            f"Elements: {self.items}")

@dataclass
class Surface:
    """
    A surface consist of Nset and Elset
    1- Nset: defines the nodes on the surface
    2- Elset: defines the elements on the surface (For pressure -> element surface)
    3- EType: defines the element type of the surface elements
    """
    # This can be "Surface Elements" or "Element Faces"
    # There could be more than 1 (2) surface defined for the same nset. (Shell Elements)
    # !! But their etype ids must be different!!
    # Thus, distinguish them with etype id.
    nset : Nset | None = None
    #elset : Elset | None = None
    etype : element.EType | None = None 
    name : str = ""
    midx : int = 0
    application : Literal["Pressure", "Contact"] = "Pressure"
    local_csys_id : int | None = None
    ignore_at_execution = False
    
@dataclass
class LocalCoordinateSystem:
    midx : int = 10
    type = "Cartesian"
    origin : tuple = (0.0, 0.0, 0.0)
    rot_x : float = 0.0
    rot_y : float = 0.0
    rot_z : float = 0.0
    ignore_at_execution = False

@dataclass
class Vector:
    midx : int

class Analysis:
    
    def __init__(self):
        
        # What makes up a load history throughout load steps?
        # 1) If load object is applies to the same node set/surface/element set
        # 2) If load type is the same (Force/Moment/Pressure/Displacement)
        # 3) If direction is the same (X/Y/Z/RX/RY/RZ/ALL)
        # Then, we can track the load history for that load objects throughout load steps linking them together.
        self._load_history : List[LoadStep] = []
        self._load_history_table : Dict[LoadStep, List[int]] = {}   

    def _add_bc(self, bc : Force | Moment | Pressure | Displacement):
        """Adds a new bc row to the load history table with a new bc object.
        If step number > 1, it sets previous loadstep values zero for the new bc."""
        pass

    def _update_bc(self, bc, loadstep : LoadStep):
        """Updates an existing bc in the load history table for the given loadstep due to NEW/ADD/DELETE logic."""
        pass

    def _add_loadstep(self, loadstep : LoadStep):
        """Adds and empty 'load step column' to the load history table.
        Automatically inherits the same values that all existing bcs have in the previous steps."""

        self._load_history.append(loadstep)
 
        """if len(self._load_history_table) == 0:
            self._load_history_table[name] = []
        else:
            last_key = next(reversed(self._load_history_table))
            self._load_history_table[name] = self._load_history_table[last_key].copy()"""

    def _get_force_history(self, nodes : str | int, direction : Literal["X", "Y", "Z"]) -> list:
        
        history = [0] # for "STEPPED" loading, it must be start from a positive value. Later work!

        for ls in self._load_history:
            for force in ls.forces:
                # Filter "nodes name" and "direction":
                if force[0] == nodes and force[1] == direction:
                    # "Add-New-Delete" Logic:
                    if force[3] == "NEW":
                        history.append(force[2])
                    elif force[3] == "ADD":
                        history.append(history[-1] + force[2])
                    elif force[3] == "DELETE":
                        history.append(0)
                else: history.append(history[-1])
        return history
        
    def _get_moment_history(self, nodes : str | int, direction : Literal["RX", "RY", "RZ"]) -> list:
        
        history = [0] # for "STEPPED" loading, it must be start from a positive value. Later work!

        for ls in self._load_history:
            for moment in ls.moments:
                if moment[0] == nodes and moment[1] == direction:
                    # Add-New-Delete Logic
                    if moment[3] == "NEW":
                        history.append(moment[2])
                    elif moment[3] == "ADD":
                        history.append(history[-1] + moment[2])
                    elif moment[3] == "DELETE":
                        history.append(0)
                else: history.append(history[-1])
        return history

    def _get_pressure_history(self, surface : str | int) -> list:
        
        history = [0] # for "STEPPED" loading, it must be start from a positive value. Later work!

        for ls in self._load_history:
            for pressure in ls.pressures:
                if pressure[0] == surface:
                    # Add-New-Delete Logic
                    if pressure[3] == "NEW":
                        history.append(pressure[2])
                    elif pressure[3] == "ADD":
                        history.append(history[-1] + pressure[2])
                    elif pressure[3] == "DELETE":
                        history.append(0)
                else: history.append(history[-1])
        return history

    """def get_boundary_history(self, nodes : str | int, direction : Literal["ALL", "X", "Y", "Z", "RX", "RY", "RZ"]) -> list:
        
        history = []
        for ls in self._load_steps:
            for disp in ls.displacements:
                if disp[0] == nodes and disp[1] == direction:
                    history.append(disp)
        return history"""

class Model:

    def __init__(self, name = '',
                       cpus = 4, 
                       work_folder = '', 
                       solver : Literal["MAPDL"] = "MAPDL"):
        
        self.solver_name = solver
        self.modelname = name
        self.cpus = cpus
        self.work_folder = work_folder

        self._analysis = Analysis()

        self._solution = None

        self._command_pipeline = []

        self._model_definitions = {"Element Types" : {},
                                   "Materials" : {},
                                   "Sets" : {"Node Sets" : {}, 
                                             "Element Sets": {}
                                            },
                                   "Surfaces" : {},
                                   "Local Coordinate Systems" : {},
                                   "Sections" : {},
                                   "Constraints" : {"Linear Couplings" : {},
                                                    "Constraint Equations" : {},
                                                    "MPC" : {},
                                                    "Joints" : {}
                                                    },
                                   "Contacts" : {},
                                   "Load Steps" : {}
                                    } 

        if work_folder == '':
            _current_path = Path.cwd()
            timestamp = datetime.now().strftime("%d.%m.%Y-%H.%M.%S")
            self.creation_time = time.time()
            folder_name = f"{name}-{timestamp}"
            self.working_directory = _current_path / folder_name
            self.working_directory.mkdir()

    def get_nset(self, name : str) -> Nset | None:
        if name in self._model_definitions["Sets"]["Node Sets"]:
            return self._model_definitions["Sets"]["Node Sets"][name]
        else: return None

    def get_elset(self, name : str):
        if name in self._model_definitions["Sets"]["Element Sets"]:
            return self._model_definitions["Sets"]["Element Sets"][name]
        else: return None

    def get_surface(self, name : str):
        if name in self._model_definitions["Sets"]["Surfaces"]:
            return self._model_definitions["Sets"]["Surfaces"][name]
        else: return None

    def import_mesh(self, file : str):
        """
        >> Reads:
            # EBLOCK
            # NBLOCK
            # ELEMENT ATTRIBUTES:
                ET, R, MP, SECTYPE, SECDATA 
            # COMPONENTS    

        >> Does not read:
            Loads, BCs, Load Steps, constraints (CP, CE, MPC) etc. 
        """

        self.mesh = mesh._ansys_mesh(file)

        for name, nlist in self.mesh.nset_dict.items():
            self.add_node_set(name, nlist, ignore_at_execution=True)

        for name, elist in self.mesh.elset_dict.items():
            self.add_element_set(name, elist, ignore_at_execution=True)

        # Create element types
        # TODO:
        # Define "KeyOpts"
        # Create "Real Constants" after import,
        etypes = self.mesh.get_mapdl_types() # [(etype, etype_idx) ...]
        if etypes is not None:
            for item in etypes:
                # item[0] -> type id
                # item[1] -> type
                # Create element type object from structuralanalysistoolbox.mapdl.element
                etype = element.VALID_ETYPES[item[1]]()   
                etype.midx = item[0]
                self.add_element_type(etype, ignore_at_execution=True)

        # Map solver command:
        #self._add_command("_execute_mesh_data", self.mesh)

    def export_model(self, fileName : str):
        """Exports file in blocked format"""
        self.mapdl.allsel()
        self.mapdl.cdwrite(option="ALL", fname=fileName, fmat="BLOCKED")
     
    def add_element_type(self, etype, **kwargs):
        """Add element type object to the model"""

        if etype.midx != 0:
            self._model_definitions["Element Types"][etype.name] = etype
        else:
            etype.midx = self._find_max_id(self._model_definitions["Element Types"]) + 1
            self._model_definitions["Element Types"][etype.name] = etype

        if kwargs: # ignore_at_execution
            etype.ignore_at_execution = kwargs.get("ignore_at_execution", False)

    def add_node_set(self, name : str, nodes : np.ndarray | list, **kwargs) -> Nset:
        if name is None:
            raise ParameterError("Node set name must be defined!")
        
        nset = Nset(name, nodes)
        nset.midx = self._find_max_id(self._model_definitions["Sets"]["Node Sets"]) + 1
        self._model_definitions["Sets"]["Node Sets"][nset.name] = nset

        if kwargs: # ignore_at_execution
            nset.ignore_at_execution = kwargs.get("ignore_at_execution", False)

        return nset

    def add_element_set(self, name : str, elements : np.ndarray | list, **kwargs) -> Elset:
        if name is None:
            raise ParameterError("Element set name must be defined!")
        
        elset = Elset(name, elements)
        elset.midx = self._find_max_id(self._model_definitions["Sets"]["Element Sets"]) + 1
        self._model_definitions["Sets"]["Element Sets"][elset.name] = elset

        if kwargs: # ignore_at_execution
            elset.ignore_at_execution = kwargs.get("ignore_at_execution", False)

        return elset

    def add_surface(self, nset : str, name : str | None = None, 
                    element_type : Any | None = None,
                    application : Literal["Pressure", "Contact"] = "Pressure",
                    local_csys_id : int | None = None,
                    **kwargs) -> Surface:
        """Create a surface object through given Nset, return surface reference"""

        if name is None:
            raise ParameterError("Surface name must be defined!")

        surf = Surface(nset=self.get_nset(nset), name=name, application=application, local_csys_id=local_csys_id)
        surf.midx = self._find_max_id(self._model_definitions["Surfaces"]) + 1

        if element_type is not None:
            self.add_element_type(element_type)
            surf.etype = element_type

        self._model_definitions["Surfaces"][surf.name] = surf

        if kwargs: # ignore_at_execution
            surf.ignore_at_execution = kwargs.get("ignore_at_execution", False)

        return surf

    def add_local_csys(self, name : str = '',
                       origin : tuple = (0.0, 0.0, 0.0), 
                       rot_x : float = 0.0, rot_y : float = 0.0, rot_z : float = 0.0, **kwargs) -> LocalCoordinateSystem:
        """Adds a local coordinate system to the model."""
        csys = LocalCoordinateSystem(origin=origin, rot_x=rot_x, rot_y=rot_y, rot_z=rot_z)
        csys.midx += (self._find_max_id(self._model_definitions["Local Coordinate Systems"]) + 1)

        if name == '':
            name = f"CSYS-{csys.midx}"

        self._model_definitions["Local Coordinate Systems"][name] = csys

        if kwargs: # ignore_at_execution
            csys.ignore_at_execution = kwargs.get("ignore_at_execution", False)

        return csys

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

        return mat
  
    def assign_material(self, elset : str, material : Material):

        elset_obj = self.get_elset(elset)
        elset_obj.properties["Material"] = material
        self._add_command("_assign_material", elset_obj)

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
    
    def add_MPC_Rigid(self, dependent : str, 
                            independent : str,
                            type : Literal["Link", "Beam"] = "Beam"):                
        
        _dependent = self.get_nset(dependent)
        _independent = self.get_nset(independent)

        mpc_obj : constraint.MPC 
        if type == "Link":
            rigid_link = constraint.MPCRigidLink(dependent=_dependent, independent=_independent)
            rigid_link.midx = self._find_max_id(self._model_definitions["Constraints"]["MPC"]) + 1
            mpc_obj = self._model_definitions["Constraints"]["MPC"][f"MPC-{rigid_link.midx}"] = rigid_link 
        elif type == "Beam":
            mpc = constraint.MPCRigidBeam(dependent=_dependent, independent=_independent)
            mpc.midx = self._find_max_id(self._model_definitions["Constraints"]["MPC"]) + 1
            mpc_obj = self._model_definitions["Constraints"]["MPC"][f"MPC-{mpc.midx}"] = mpc
        else: return

        # Create a new MPC184 Element and Element Type ID
        # We need this to distinguish MPC object id and MPC Element Type id
        mpc_obj.etype = type_mpc = element.MPC184()
        self.add_element_type(type_mpc)

    def add_joint(self, joint : constraint.MPCJoint):

        # Map solver command
        self._add_command("_create_joint", joint)
    
    def add_symmetry(self):
        # TODO: symmetric / periodic
        pass   

    def add_interface(self):
        # Creates interface for various use cases.
        pass 

    def add_contact(self, master : Nset | str = '', slave : Nset | str = '') -> contact.Contact:
        pass

    def add_line_element(self):
        """Add beam, link, spring-damper"""
        pass

    def add_point_element(self):
        """Add mass, matrix"""
        pass

    def add_loadstep(self, 
                     name : str,
                     analysis : Literal["STATIC", "BUCKLE", "MODAL", "HARMIC" ,"TRANS", "SUBSTR"] = "STATIC",
                     status : Literal["NEW", "RESTART"] = "NEW",
                     step : int = 0, substep : int = 0) -> LoadStep:
        
        ls = LoadStep(name=name,
                      model=self,
                      end_time = 1,
                      analysis=analysis, status=status, step=step, substep=substep)
        
        ls.step_number = self._find_max_id(self._model_definitions["Load Steps"]) + 1
        self._model_definitions["Load Steps"][name] = ls
        #self._solution._add_loadstep(ls)        

        return ls

    def remove(self, group : str, item : str) -> None | str:
        '''Removes items from the model'''
        try:
            self._model_definitions[group].pop(item)
        except KeyError:
            raise KeyError(f"Item '{group}-{item}' does not exist")
        return None

    def plot_mesh(self):
        
        if self.mesh:
            plotter = Plotter()
            plotter.show_mesh(self.mesh.unstructuredgrid)

    def plot_load_history(self, nodes : str, load_type : Literal["Force", "Moment", "Displacement"],
                  direction : Literal["X", "Y", "Z", "RX", "RY", "RZ"]):
        ls_plot = LoadStepPlot(self._analysis, nodes=nodes, load_type=load_type, direction=direction)
        ls_plot.plot()

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

    def _add_command(self, func_name : str, *args):
        """
        adds command to << self._command_map >> list to be executed at the pre-processing time.
        func_name is corresponding mapdl commands function in the mapdl.py
        model_item_object is any kind of object that has information to convey to 
        mapdl commands function.
        """
        if self.solver_name == "MAPDL":
            solver_func = getattr(command, func_name) # Module must be pre-loaded!
        else: raise NameError("command not found in the module!")
        self._command_pipeline.append((solver_func, *args))

    def _execute_commands(self, *args):
        """if args[0] == False, it skips SOLVE commands and only writes other mapdl commands"""
        if self.solver_name == "MAPDL":
            self.mapdl = self._start()

        """for command in self._command_pipeline:
            print(command)"""

        for command in self._command_pipeline:

            if args[0] == False and command[0].__name__ == "_solve":
                #self._add_command("_write_comment", "------ SOLVE COMMAND SKIPPED!!! -----") 
                continue # Skip "Solve" Command

            if len(command) > 1:
                command[0](self.mapdl, command[1])
            else: 
                command[0](self.mapdl)

    def _create_command_pipeline(self):
        """
        Modeling Pipeline:
        --------------------
        ---- Start Solver
        2) Import Mesh Data
        3) Element Types
        4) Local Coordinate Systems
        5) Materials
        6) Sets (Nsets, Elsets)
        7) Surfaces
        8) Sections
        9) Constraints (CP, CE, MPC, Joints)
        10) Contacts    
        11) Load Steps
        ---- Solve
        """
        self._add_command("_execute_mesh_data", self.mesh)

        for etype in self._model_definitions["Element Types"].values():
            if not etype.ignore_at_execution: self._add_command("_create_element_type", etype)

        for cs in self._model_definitions["Local Coordinate Systems"].values():
            if not cs.ignore_at_execution: self._add_command("_local_csys", cs)
            
        for mat in self._model_definitions["Materials"].values():
            self._add_command("_material", mat)

        for nset in self._model_definitions["Sets"]["Node Sets"].values():
            if not nset.ignore_at_execution: self._add_command("_create_node_set", nset)

        for elset in self._model_definitions["Sets"]["Element Sets"].values():
            if not elset.ignore_at_execution: self._add_command("_create_element_set", elset)

        for surf in self._model_definitions["Surfaces"].values():
            if not surf.ignore_at_execution: self._add_command("_create_surface", surf)
            
        # Creating Sections Skipped

        # Creating CP Skipped
        # Creation CE Skipped
        # Creation Joints Skipped
        for mpc_rigid in self._model_definitions["Constraints"]["MPC"].values():
            self._add_command("_MPC184", mpc_rigid)

        # Creation Contacts Skipped

        for load_step in self._model_definitions["Load Steps"].values():
            self._add_command("_load_step", load_step)
            self._add_command("_solve")
            
    def solve(self):
        self._create_command_pipeline()
        self._add_command("_exit_mapdl")
        self._execute_commands(True)
        self._command_pipeline.clear()

    def write_input_file(self):
        self._create_command_pipeline()
        self._add_command("_write_input_file")
        self._add_command("_exit_mapdl")
        self._execute_commands(False) # Skip Solve Commands
        self._command_pipeline.clear()

    def _start(self):
        self.mapdl = launch_mapdl(run_location=self.working_directory,
                                  jobname=self.modelname,
                                  nproc=self.cpus,
                                  override=True)
        
        # First, reset the MAPDL database.
        self.mapdl.clear()
        return self.mapdl
    
    def _stop(self):
        # TODO : move it to the solver module
        self.mapdl.exit()

    def solution(self) -> Any:
        """Returns DPF solution object after solving the model."""
        _result_file_path = Path(self.working_directory / f"{self.modelname}.rst")

        solution = None
        if self.solver_name == "MAPDL":
            from structuralanalysistoolbox.post.mapdl import solutions
            solution = solutions.DpfSolution(_result_file_path)

        return solution.simulation








    





    