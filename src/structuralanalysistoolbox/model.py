
from __future__ import annotations
import ansys.mapdl.core as core
from ansys.mapdl.core import launch_mapdl
from dataclasses import dataclass, field
from typing import Literal, overload
from enum import Enum
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
from structuralanalysistoolbox.loading.loadstep import LoadStep
from structuralanalysistoolbox.plotting.plots import LoadStepPlot
from structuralanalysistoolbox.visualization.meshplotter import Plotter
from structuralanalysistoolbox.exceptions import ParameterError


@dataclass
class Nset:
    name : str
    items : np.ndarray | list | None
    midx : int = 0
    ignore_at_execution : bool = False

@dataclass
class Elset:
    name : str
    items : np.ndarray | list
    midx = 0
    properties : Dict[str, Material | None] = field(
        default_factory=lambda : {"Material" : None})
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
    etype : element.EType | None = None 
    name : str = ""
    midx : int = 0
    application : Literal["Pressure", "Contact"] = "Pressure"
    local_csys_id : int | None = None
    ignore_at_execution : bool = False
    
@dataclass
class LocalCoordinateSystem:
    midx : int = 10
    type = "Cartesian"
    origin : tuple = (0.0, 0.0, 0.0)
    rot_x : float = 0.0
    rot_y : float = 0.0
    rot_z : float = 0.0
    ignore_at_execution : bool = False

@dataclass
class PilotNode:
    """A node set without node ID. It has node set name and node coordinates
    in order to refer a node, which is going to be created during solver stage,
    with its name."""
    nset : Nset
    x : float
    y : float
    z : float
    local_cs : LocalCoordinateSystem | None = None

@dataclass
class Vector:
    midx : int

class Analysis:
    pass    

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

        # Attribute Lists
        self._material_list = []
        self._element_type_list = []
        self._real_constants_list = []
        self._local_cs_list = []
        self._section_list = []

        # Load Step List
        self._load_step_list : List[LoadStep] = []

        self._pilot_node_list = []

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
                                                    "Surface-Based" : {},
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

#####################
## IMPORT/EXPORT
#####################

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
                etype.ignore_at_execution=True
                self._add_attribute("Type", etype, etype.midx)
                self._model_definitions["Element Types"][etype.name] = etype

        # Map solver command:
        #self._add_command("_execute_mesh_data", self.mesh)

    def export_model(self, fileName : str):
        """Exports file in blocked format"""
        self.mapdl.allsel()
        self.mapdl.cdwrite(option="ALL", fname=fileName, fmat="BLOCKED")

#####################
## ATTRIBUTES
#####################    

    def add_material(self, mat : str | Material , elset : str) -> Material:
        """Adds a material to the model tree and assing to an elset."""
        
        elset_obj = self.get_elset(elset)
        if elset_obj is None:
            raise ValueError("Elset can not found!")

        if isinstance(mat, str):
            mat = material.load(mat_name=mat)
            # Add to both attribute list & model tree
            mat.midx = self._add_attribute(attribute_type="Mat", attribute=mat)
            self._model_definitions["Materials"][mat.name] = mat
        elif isinstance(mat, Material):
            # Add to both attribute list & model tree
            mat.midx = self._add_attribute(attribute_type="Mat", attribute=mat)
            self._model_definitions["Materials"][mat.name] = mat
        else: raise ValueError("mat value type is not str or Material")

        # Assign a material to an Elset
        if isinstance(elset_obj, Elset): 
            elset_obj.properties["Material"] = mat

        return mat

    def add_coordinate_system(self, name : str = '',
                       origin : tuple = (0.0, 0.0, 0.0), 
                       rot_x : float = 0.0, rot_y : float = 0.0, rot_z : float = 0.0, **kwargs) -> LocalCoordinateSystem:
        """Adds a local coordinate system to the model."""
        csys = LocalCoordinateSystem(origin=origin, rot_x=rot_x, rot_y=rot_y, rot_z=rot_z)

        if len(self._local_cs_list) == 0:
            csys.midx = 11
        else:
            csys.midx = self._add_attribute("Csys", csys)
            
        # Add a name for the model tree
        if name == '': name = f"CSYS-{csys.midx}"

        self._model_definitions["Local Coordinate Systems"][name] = csys

        if kwargs and kwargs["ignore_at_execution"]: # ignore_at_execution
            csys.ignore_at_execution = True

        return csys

    def add_unit_vector(self):
        pass

#####################
## SETS
#####################

    def add_node_set(self, name : str, nodes : np.ndarray | list | None, **kwargs) -> Nset:
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

        if kwargs and kwargs["ignore_at_execution"]: # ignore_at_execution
            surf.ignore_at_execution = True

        return surf

    def add_pilot_node(self, name : str, x : float, y : float, z : float, local_cs : LocalCoordinateSystem | None = None) -> PilotNode:

        nset = self.add_node_set(name=name, nodes=None, ignore_at_execution=True)
        pilot_node = PilotNode(nset=nset, x=x, y=y, z=z, local_cs=local_cs)
        self._pilot_node_list.append(pilot_node)

        return pilot_node

    def get_nset(self, name : str) -> Nset | None:
        if name in self._model_definitions["Sets"]["Node Sets"]:
            return self._model_definitions["Sets"]["Node Sets"][name]
        else: return None

    def get_elset(self, name : str) -> Elset | None:
        if name in self._model_definitions["Sets"]["Element Sets"]:
            return self._model_definitions["Sets"]["Element Sets"][name]
        else: return None

    def get_surface(self, name : str) -> Surface | None:
        if name in self._model_definitions["Sets"]["Surfaces"]:
            return self._model_definitions["Sets"]["Surfaces"][name]
        else: return None

#####################
## CONSTRAINTS
#####################

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
    
    def add_MPC_rigid_link_web(self):
        pass

    def add_MPC_rigid_beam_web(self):
        pass

    def add_force_dist_surf_constraint(self, pilot_node : str, contact_nodes : str,
                     dof : tuple = ("UX", "UY", "UZ", "ROTX", "ROTY", "ROTZ"),
                     weight : Literal["Contact Area", "RBE3"] | np.ndarray = "Contact Area",
                     contact_model : Literal["Surface-To-Surface", "Node-To-Surface"] = "Node-To-Surface",
                     contact_algorithm : Literal["MPC", "Lagrange"] = "MPC",
                     local_cs : LocalCoordinateSystem | None = None,
                     bonding_type : Literal["Bonded (always)", "Bonded (initial)"] = "Bonded (always)",
                     symmetry : Literal["Pilot XY Plane",
                                        "Pilot XZ Plane",
                                        "Pilot YZ Plane",
                                        "Pilot XY + XZ",
                                        "Pilot XY + YZ",
                                        "Pilot XZ + YZ"] | None = None) -> constraint.ForceDistributedConstraint:
        """MPC Deformable Surface
        * Weighting factor can be specified in the form of np.ndarray[node_id, weight] or 1, or None for internal
        * 
        
        """
        mpc = constraint.ForceDistributedConstraint(pilot_node=pilot_node, contact_nodes=contact_nodes,
                                             constrained_dof=dof,
                                             contact_model=contact_model,
                                             contact_algorithm=contact_algorithm,
                                             local_cs=local_cs,
                                             bonding_type=bonding_type,
                                             constrained_surface_symmetry=symmetry)
        
        # Add internally generated contact elements to the model
        mpc.contact_element_type.midx = self._add_attribute("Type", mpc.contact_element_type)
        mpc.target_element_type.midx = self._add_attribute("Type", mpc.target_element_type)

        # Create an empty material attribute
        mpc.material_midx = self._add_attribute(attribute_type="Mat", attribute=None)

        # Create and Associate the Real Constant Object
        mpc.real_constants = element.ContaRealConstants()
        mpc.real_constants.midx = self._add_attribute(attribute_type="Real" ,attribute=mpc.real_constants)

        # Define weighten logic
        if weight == "RBE3":
            mpc.weighting_factor = 1
            mpc.target_element_type.weighting_factor = "1"
        elif weight == "Contact Area":
            mpc.target_element_type.weighting_factor = "Internal"
        elif isinstance(weight, np.ndarray):
            mpc.weighting_factor = weight
            mpc.target_element_type.weighting_factor = "User Defined"
            # modify real constant fkn after Elements created
        
        mpc.midx = self._find_max_id(self._model_definitions["Constraints"]["Surface-Based"]) + 1
        self._model_definitions["Constraints"]["Surface-Based"][f"Force-Dist Surface-{mpc.midx} ({contact_algorithm})"] = mpc
        return mpc

    def add_rigid_surface_constraint(self, pilot_node : str, contact_nodes : str,
                     dof : tuple = ("UX", "UY", "UZ", "ROTX", "ROTY", "ROTZ"),
                     contact_model : Literal["Surface-To-Surface", "Node-To-Surface"] = "Node-To-Surface",
                     contact_algorithm : Literal["MPC", "Lagrange"] = "MPC",
                     local_cs : LocalCoordinateSystem | None = None,
                     bonding_type : Literal["Bonded (always)", "Bonded (initial)"] = "Bonded (always)") -> constraint.RigidSurfaceConstraint:

        mpc = constraint.RigidSurfaceConstraint(pilot_node=pilot_node, contact_nodes=contact_nodes,
                                             constrained_dof=dof,
                                             contact_model=contact_model,
                                             contact_algorithm=contact_algorithm,
                                             local_cs=local_cs,
                                             bonding_type=bonding_type)
        
        # Add internally generated contact elements to the model
        mpc.contact_element_type.midx = self._add_attribute("Type", mpc.contact_element_type)
        mpc.target_element_type.midx = self._add_attribute("Type", mpc.target_element_type)

        # Create an empty material attribute
        mpc.material_midx = self._add_attribute(attribute_type="Mat", attribute=None)

        # Create and Associate the Real Constant Object
        mpc.real_constants = element.ContaRealConstants()
        mpc.real_constants.midx = self._add_attribute(attribute_type="Real" ,attribute=mpc.real_constants)
        
        mpc.midx = self._find_max_id(self._model_definitions["Constraints"]["Surface-Based"]) + 1
        self._model_definitions["Constraints"]["Surface-Based"][f"Rigid Surface-{mpc.midx} ({contact_algorithm})"] = mpc
        return mpc

    def add_coupled_surface_constraint(self, pilot_node : str, contact_nodes : str,
                     dof : tuple = ("UX", "UY", "UZ", "ROTX", "ROTY", "ROTZ"),
                     contact_model : Literal["Surface-To-Surface", "Node-To-Surface"] = "Node-To-Surface",
                     local_cs : LocalCoordinateSystem | None = None,
                     bonding_type : Literal["Bonded (always)", "Bonded (initial)"] = "Bonded (always)") -> constraint.CoupledSurfaceConstraint:

        mpc = constraint.CoupledSurfaceConstraint(pilot_node=pilot_node, contact_nodes=contact_nodes,
                                             constrained_dof=dof,
                                             contact_model=contact_model,
                                             local_cs=local_cs,
                                             bonding_type=bonding_type)
        
        # Add internally generated contact elements to the model
        mpc.contact_element_type.midx = self._add_attribute("Type", mpc.contact_element_type)
        mpc.target_element_type.midx = self._add_attribute("Type", mpc.target_element_type)

        # Create an empty material attribute
        mpc.material_midx = self._add_attribute(attribute_type="Mat", attribute=None)

        # Create and Associate the Real Constant Object
        mpc.real_constants = element.ContaRealConstants()
        mpc.real_constants.midx = self._add_attribute(attribute_type="Real" ,attribute=mpc.real_constants)
        
        mpc.midx = self._find_max_id(self._model_definitions["Constraints"]["Surface-Based"]) + 1
        self._model_definitions["Constraints"]["Surface-Based"][f"Coupled Surface-{mpc.midx} (MPC)"] = mpc
        return mpc

    def add_joint(self, joint : constraint.MPCJoint):

        # Map solver command
        self._add_command("_create_joint", joint)
    
    def add_rigid_body(self):
        pass

    def add_symmetry(self):
        # symmetric / periodic
        pass   

#####################
## CONTACT
#####################

    def add_contact(self, master : Nset | str = '', slave : Nset | str = '') -> contact.Contact:
        pass

##########################
## CREATE SIMPLE ELEMENTS
##########################

    def add_line_element(self):
        """Add beam, link, spring-damper"""
        pass

    def add_point_element(self):
        """Add mass, matrix"""
        pass

#####################
## LOAD STEP
#####################

    def add_loadstep(self, 
                     name : str,
                     analysis : Literal["STATIC", "BUCKLE", "MODAL", "HARMIC" ,"TRANS", "SUBSTR"] = "STATIC",
                     status : Literal["NEW", "RESTART"] = "NEW",
                     step : int = 0, substep : int = 0) -> LoadStep:
        
        ls = LoadStep(name=name,
                      model=self,
                      end_time = 1,
                      analysis=analysis, status=status, step=step, substep=substep)
        
        self._load_step_list.append(ls)

        ls.step_number = self._find_max_id(self._model_definitions["Load Steps"]) + 1
        self._model_definitions["Load Steps"][name] = ls
        #self._solution._add_loadstep(ls)        

        return ls

#####################
## INFO
#####################

    def info(self, item : None | str = None) -> str:
        """Returns detailed information string for model item object."""
        if item == None:
            # List model items as a tree structure
            self._print_model_tree(self._model_definitions, indent=" ")
        else:
            obj = self._find_item_object(self._model_definitions, item)
            print(obj.__str__())

    def plot_mesh(self):
        
        if self.mesh:
            plotter = Plotter()
            plotter.show_mesh(self.mesh.unstructuredgrid)

    def plot_load_history(self, 
                          nset : str, 
                          load_type : Literal["Force", "Moment", "Pressure", "Displacement"],
                          direction : Literal["X", "Y", "Z", "RX", "RY", "RZ"]):
        
        history : list = []
        for ls in self._load_step_list:
            if load_type == "Force":
                for force in  ls.forces:
                    if force.set.name == nset and force.direction == direction:
                        history.append(force)
            elif load_type == "Moment":
                for moment in  ls.moments:
                    if moment.set.name == nset and moment.direction == direction:
                        history.append(moment)
            elif load_type == "Pressure":
                for pressure in  ls.pressures:
                    if pressure.surf.nset.name == nset and pressure.direction == direction:
                        history.append(pressure)
            elif load_type == "Displacement":
                for displacement in  ls.displacements:
                    if displacement.set.name == nset and displacement.direction == direction:
                        history.append(displacement)

        load_list = np.zeros(len(self._load_step_list), dtype=np.float64)
        load_list[0] = history[0].value
        for idx, load in enumerate(history[1:]):
            if load.operation == "NEW":
                load_list[idx+1] = load.value
            elif load.operation == "ADD":
                load_list[idx+1] = load_list[idx] + load.value
            elif load.operation == "DELETE":
                load_list[idx+1] = 0
            else:
                load_list[idx+1] = 0

        ls_plot = LoadStepPlot(load_list, nodes=nset, load_type=load_type, direction=direction)
        ls_plot.plot()

#####################
## MISC.
#####################

    def remove(self, group : str, item : str) -> None | str:
        '''Removes items from the model'''
        try:
            self._model_definitions[group].pop(item)
        except KeyError:
            raise KeyError(f"Item '{group}-{item}' does not exist")
        return None

#####################
## _private
#####################

    def _get_item_object(self, path : tuple[str]):
        """Returns object for any model item"""
        
        node = self._model_definitions
        for key in path:
            if not isinstance(node, dict) or key not in node:
                return None
            node = node[key]
        return node

    def _find_max_attribute_id(self, attribute_type : Literal["Mat", "Type", "Real", "Csys", "Secn"]) -> int:

        att_list : list
        if attribute_type == "Mat":
            att_list = self._material_list
        elif attribute_type == "Type":
            att_list = self._element_type_list
        elif attribute_type == "Real":
            att_list = self._real_constants_list
        elif attribute_type == "Csys":
            att_list = self._local_cs_list
        elif attribute_type == "Secn":
            att_list = self._section_list

        if len(att_list) == 0:
            return 0
        else:
            max_midx = 0
            for att in att_list:
                if att[0] > max_midx: max_midx = att[0]
            return max_midx

    def _add_attribute(self, attribute_type : Literal["Mat", "Type", "Real", "Csys", "Secn"], 
                       attribute : None | Any, 
                       midx : int | None = None) -> int:
        """An empty attribute can be add as a placeholder in the form of (id, None)
           Returns attribute id"""
        att_list : list
        if attribute_type == "Mat":
            att_list = self._material_list
        elif attribute_type == "Type":
            att_list = self._element_type_list
        elif attribute_type == "Real":
            att_list = self._real_constants_list
        elif attribute_type == "Csys":
            att_list = self._local_cs_list
        elif attribute_type == "Secn":
            att_list = self._section_list

        if midx:
            att_list.append((midx, attribute))
        else:
            midx = self._find_max_attribute_id(attribute_type) + 1
            att_list.append((midx, attribute))

        return midx

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

    def _print_model_tree(self, tree, indent=""):
        if not isinstance(tree, dict):
            return

        last_key = list(tree.keys())[-1] if tree else None

        for key in tree:
            is_last = (key == last_key)
            prefix = "└── " if is_last else "├── "
            print(indent + prefix + key)

            next_indent = indent + ("    " if is_last else "│   ")
            if isinstance(tree[key], dict):
                self._print_model_tree(tree[key], next_indent)

#####################
## EXECUTION
#####################

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
        1) Import Mesh Data
        2) Pilot Nodes
        3) Element Types
        4) Real Constants
        5) Local Coordinate Systems
        6) Materials
        7) Sets (Nsets, Elsets)
        8) Surfaces
        9) Sections
        10) Constraints (CP, CE, MPC, Surface-Based, Joints)
        11) Contacts    
        12) Load Steps
        ---- Solve
        """
        self._add_command("_execute_mesh_data", self.mesh)

        """for etype in self._model_definitions["Element Types"].values():
            if not etype.ignore_at_execution: self._add_command("_create_element_type", etype)"""
        for etype in self._element_type_list:
            if not etype[1].ignore_at_execution: self._add_command("_create_element_type", etype[1])

        for node in self._pilot_node_list:
            self._add_command("_create_pilot_node", node)

        for rc in self._real_constants_list:
            self._add_command("_set_real_constants", rc)

        for cs in self._local_cs_list:
            if not cs.ignore_at_execution: self._add_command("_local_csys", cs)
            
        for mat in self._model_definitions["Materials"].values():
            self._add_command("_material", mat)

        for nset in self._model_definitions["Sets"]["Node Sets"].values():
            if not nset.ignore_at_execution: self._add_command("_create_node_set", nset)

        for elset in self._model_definitions["Sets"]["Element Sets"].values():
            if not elset.ignore_at_execution: 
                self._add_command("_create_element_set", elset)
                self._add_command("_assign_material", elset)

        for surf in self._model_definitions["Surfaces"].values():
            if not surf.ignore_at_execution: self._add_command("_create_surface", surf)
            
        # Creating Sections Skipped

        # Creating CP Skipped
        # Creation CE Skipped
        # Creation Joints Skipped

        for mpc in self._model_definitions["Constraints"]["MPC"].values():
            pass            

        for surf_con in self._model_definitions["Constraints"]["Surface-Based"].values():
            if isinstance(surf_con, constraint.ForceDistributedConstraint):
                self._add_command("_create_force_distirbuted_constraint", surf_con)
            elif isinstance(surf_con, constraint.RigidSurfaceConstraint):
                self._add_command("_create_rigid_surface_constraint", surf_con)
            elif isinstance(surf_con, constraint.CoupledSurfaceConstraint):
                self._add_command("_create_coupled_surface_constraint", surf_con)

        # Creation Contacts Skipped

        for load_step in self._load_step_list:
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








    





    