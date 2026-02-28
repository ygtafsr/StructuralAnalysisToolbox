
from __future__ import annotations
import ansys.mapdl.core as core
from ansys.mapdl.core import launch_mapdl
from dataclasses import dataclass, field
from typing import Literal, overload, Dict, List, Any
from pathlib import Path
import numpy as np
import time
from datetime import datetime

import pandas as pd

from structuralanalysistoolbox.data import mesh
from structuralanalysistoolbox.materials import material 
from structuralanalysistoolbox.materials.material import Material
from structuralanalysistoolbox.constraints import constraint, contact
from structuralanalysistoolbox.mapdl import command, element, attributes
from structuralanalysistoolbox.loading.loadstep import LoadStep, Force, Moment, Pressure, Displacement
from structuralanalysistoolbox.plotting.plots import LoadStepPlot
from structuralanalysistoolbox.visualization.meshplotter import Plotter
from structuralanalysistoolbox.exceptions import ParameterError


"""
[1] An 'Add_' method adds ONLY corresponding model item to the tree.
    For example, add_surface only add surface not the internal element types etc.
    add_mpc_constraint do not adds surfaces and element types.

[2] Assigning Attributes::

    [User] :: Through element sets
    [Me] :: Set attirbutes inside any model item implementation

[3] Restart Analysis:
    Solve method looks up model._load_step_list starting from the last loadstep to
    find out the load step with a restart option and only run this loadstep.


"""


@dataclass
class Nset:
    name : str
    items : np.ndarray | list | None
    midx : int = 0
    ignore_at_execution : bool = False

@dataclass
class Elset:
    name : str | None = None
    items : np.ndarray | list | None = None
    element_type : attributes.Etype | None = None
    real_constants : attributes.Real | None = None
    material : attributes.Mat | None = None
    element_cs : int | None = None
    section : attributes.Section | None = None
    ignore_at_execution = False
    show_in_model_tree = True

@dataclass
class Surface:
    """
    Nset: defines the nodes on the surface
    EType: defines the element type of the surface elements
    """
    # There could be more than 1 (2) surface defined for the same nset. (Shell Elements)
    # !! But their etype ids must be different!!
    # Thus, distinguish them with etype id.
    nset : Nset | None = None
    name : str | None = None
    type : Literal["NODE SET", "ELEMENT FACES", "SURFACE ELEMENTS"] = "ELEMENT FACES"

    element_type : attributes.Etype | None = None
    real_constants : attributes.Real | None = None
    material : attributes.Mat | None = None
    csys : CoordinateSystem | None = None

    ignore_at_execution : bool = False
    
@dataclass
class CoordinateSystem:
    midx : int = 0
    type = "Cartesian"
    origin : tuple = (0.0, 0.0, 0.0)
    rot_x : float = 0.0
    rot_y : float = 0.0
    rot_z : float = 0.0
    ignore_at_execution : bool = False
    show_in_model_tree = True

@dataclass
class PilotNode:
    """A node set without node ID. It has node set name and node coordinates
    in order to refer a node, which is going to be created during solver stage,
    with its name."""
    nset : Nset
    x : float
    y : float
    z : float
    local_cs : CoordinateSystem | None = None

@dataclass
class Vector:
    midx : int

class Analysis:
    
    def __init__(self, model):

        self.model : Model = model
        self.ls_table : List[List[Force | Pressure | Moment | Displacement | None]] = []

        """
          LS-1    LS-2    LS-3
        [ None,   None,   None ...]
        [ None,   None,   None ...]
         .
         .
         .
        """

    def add_new_row(self, step_number : int):
        """Add new row full of 'None' up to new load step."""
        self.ls_table.append([None for i in range(0, step_number)])
        pass

    def add_new_column(self):
        for row in self.ls_table:
            row.append(None)
        pass

    def check_columns(self, bc : Force | Moment | Pressure | Displacement) -> tuple[int, int] | None:
        for row_idx, row in enumerate(self.ls_table):
            for col_idx, cell in enumerate(row):
                if cell == bc:
                    return (row_idx, col_idx)
        return None

    def add_bc(self, bc : Force | Moment | Pressure | Displacement, step_number : int):
        """ 
        Adds a bc to the last load step. If there is a same scope for bc which is
        already defined, then add the bc this row. Otherwise, adds a new scope row.
        """
        # Get bc scope and check its existence in ls_table
        # If it exist then add bc end of the ls list of this scope line.
        # Else, add_new_scope + add add bc end of the ls list of this scope line.
        # Add a new bc object in any case that the bc is a continuation of a previous one or not.
        # Continuity check is only important to decide using a new scope or a previous one.

        #step_number = len(self.model._load_step_list)

        if step_number == 1:
            self.add_new_row(step_number)
            self.ls_table[-1][-1] = bc
        else:
            if item := self.check_columns(bc):
                self.ls_table[item[0]][step_number-1] = bc 
            else:
                self.add_new_row(step_number)
                self.ls_table[-1][-1] = bc
        pass

    def get_bc_history(self) -> pd.DataFrame:

        self.ls_table.clear()

        for ls in self.model._load_step_list:
            self.add_new_column()
            for force in ls.forces:
                self.add_bc(force, ls.step_number)
            for press in ls.pressures:
                self.add_bc(press, ls.step_number)
            for disp in ls.displacements:
                self.add_bc(disp, ls.step_number)

        def _bc_to_str(bc_obj: Force | Moment | Pressure | Displacement | None) -> str | None:
            if bc_obj is None:
                return None
            if type(bc_obj) is Force:
                scope = bc_obj.set.name
                kind = "FORCE"
            elif type(bc_obj) is Moment:
                scope = bc_obj.set.name
                kind = "MOMENT"
            elif type(bc_obj) is Pressure:
                scope = bc_obj.set.nset.name
                kind = "PRESS"
            elif type(bc_obj) is Displacement:
                scope = bc_obj.set.name
                kind = "DISP"
            else:
                scope = ""
                kind = type(bc_obj).__name__.upper()
            return f"{bc_obj.direction} : {bc_obj.value}"

        if not self.ls_table:
            return pd.DataFrame()

        num_cols = max(len(row) for row in self.ls_table)
        columns = [f"LS-{ls_no+1}" for ls_no in range(num_cols)]

        data = []
        index = []
        for row_idx, row in enumerate(self.ls_table, start=1):
            padded = list(row) + [None] * (num_cols - len(row))
            data.append([_bc_to_str(cell) for cell in padded])

            first = next((cell for cell in padded if cell is not None), None)
            if first is None:
                index.append(f"SCOPE-{row_idx}")
            elif type(first) is Pressure:
                index.append(f"PRESS : {first.set.nset.name}")
            else:
                index.append(f"{type(first).__name__.upper()} : {first.set.name}")

        return pd.DataFrame(data=data, columns=columns, index=index)
       
class Model:

    def __init__(self, name = '',
                       cpus = 4, 
                       work_folder = '', 
                       solver : Literal["MAPDL"] = "MAPDL"):
        
        self.solver_name = solver
        self.modelname = name
        self.cpus = cpus
        self.work_folder = work_folder

        self._analysis = Analysis(self)

        self._solution = None

        self._command_pipeline = []

        # Attribute Lists
        self._element_type_list = []
        self._real_constants_list = []
        self._material_list = []
        self._cs_list = []
        self._section_list = []

        # Assign attributes within Elset
        # Some elsets do not show up in the model tree
        self._element_set_list = []
        self._surface_list = []

        # Load Step List
        self._load_step_list : List[LoadStep] = []

        self._pilot_node_list = []

        self._model_definitions = {"Element Types" : {},
                                   "Materials" : {},
                                   "Sets" : {"Node Sets" : {}, 
                                             "Element Sets": {}
                                           # "Surfaces" : {}
                                            },
                                   
                                   "Coordinate Systems" : {},
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
            self._current_path = Path.cwd()
            timestamp = datetime.now().strftime("%d.%m.%Y-%H.%M.%S")
            self.folder_name = f"{name}-{timestamp}"

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

        # Get element attributes and specify in the model
        # Element Types:
        etypes = self.mesh.get_mapdl_types() # [(etype, etype_idx) ...]
        if etypes is not None:
            for item in etypes: # item[0]::type id // item[1]::type
                etype = element.VALID_ETYPES[item[1]]()  
                self.add_element_type(etype, item[0], ignore_at_execution=True) 
            #self._add_attribute(attributes.Mat(ignore_at_execution=True), midx=1) # Mat(1)
            self._add_attribute(attributes.Real(ignore_at_execution=True), midx=1) # Real(1)
            self._add_attribute(attributes.Section(ignore_at_execution=True), midx=1) # Sec(1)
            self.add_coordinate_system(name="Cartesian", id=0, ignore_at_execution=True)
            self.add_coordinate_system(name="Cylindrical", id=1, ignore_at_execution=True)
            self.add_coordinate_system(name="Spherical", id=2, ignore_at_execution=True)

    def export_model(self, fileName : str):
        """Exports file in blocked format"""
        self.mapdl.allsel()
        self.mapdl.cdwrite(option="ALL", fname=fileName, fmat="BLOCKED")

#####################
## ATTRIBUTES
#####################    

    def add_element_type(self, etype : element.Etype, midx : int | None = None, **kwargs):
        """Adds 3D, Shell, 2D, 1D, 0D Elements to the model & Model Tree"""
        # 3D, Shell, 2D, 1D, 0D Elements are shown in the model tree
        # Others that are internally created are not shown in the model three
        # they are take place in model._element_type_list
        if midx: 
            self._add_attribute(etype, midx)
        else:
            etype.midx = self._add_attribute(etype)

        self._model_definitions["Element Types"][etype.name] = etype

        if kwargs and kwargs["ignore_at_execution"]: # ignore_at_execution
            etype.ignore_at_execution = kwargs["ignore_at_execution"]
        
    def add_material(self, mat : str | Material , elset : str) -> Material:
        """Adds a material to the model tree and assing to an elset."""
        
        elset_obj = self.get_elset(elset)
        if elset_obj is None:
            raise ValueError("Elset can not found!")

        if isinstance(mat, str):
            mat = material.load(mat_name=mat)
            # Add to both attribute list & model tree
            mat.midx = self._add_attribute(mat)
            self._model_definitions["Materials"][mat.name] = mat
        elif isinstance(mat, Material):
            # Add to both attribute list & model tree
            mat.midx = self._add_attribute(mat)
            self._model_definitions["Materials"][mat.name] = mat
        else: raise ValueError("mat value type is not str or Material")

        # Assign a material to an Elset
        if isinstance(elset_obj, Elset): 
            elset_obj.material = mat

        return mat

    def assign_attribute(self, elset : str | Elset, attribute : attributes.Etype | 
                                                                attributes.Real | 
                                                                attributes.Section |
                                                                CoordinateSystem):

        elset_ : Elset
        if isinstance(elset, Elset):
            elset_ = elset
        elif isinstance(elset, str) and self.get_elset(elset):
            elset_ = self.get_elset(elset)

        if isinstance(attribute, attributes.Etype):
            elset_.element_type = attribute
        elif isinstance(attribute, attributes.Real):
            elset_.real_constants = attribute
        elif isinstance(attribute, attributes.Section):
            elset_.section = attribute
        elif isinstance(attribute, CoordinateSystem):
            elset_.element_cs = attribute

    def add_coordinate_system(self, 
                              name : str = '',
                              origin : tuple = (0.0, 0.0, 0.0), 
                              rot_x : float = 0.0, rot_y : float = 0.0, rot_z : float = 0.0, **kwargs) -> CoordinateSystem:
        """Adds a coordinate system to the model."""
        csys = CoordinateSystem(origin=origin, rot_x=rot_x, rot_y=rot_y, rot_z=rot_z)

        if "id" in kwargs.keys(): # It's a global solver coordinate system
            csys.midx = kwargs["id"]
        else:
            if max(cs.midx for cs in self._cs_list) < 11: # Start from 11 for a local cs id
                csys.midx = 11
            else:
                csys.midx = self._add_attribute(csys)
        
        if "ignore_at_execution" in kwargs.keys(): # ignore_at_execution
            csys.ignore_at_execution = kwargs["ignore_at_execution"]

        # Add a name for the model tree
        if name == '': name = f"CSYS-{csys.midx}"

        self._cs_list.append(csys)

        if "show_in_model_tree" in kwargs.keys():
            if kwargs["show_in_model_tree"]:
                self._model_definitions["Coordinate Systems"][name] = csys
        else:
            self._model_definitions["Coordinate Systems"][name] = csys

        return csys
    
    def get_coordinate_system(self, name : str) -> CoordinateSystem | None:
        ''' Returns cs object from model tree'''
        for cs_name in self._model_definitions["Coordinate Systems"].keys():
            if cs_name == name:
                return self._model_definitions["Coordinate Systems"][cs_name]
            else: return None

    def add_unit_vector(self):
        pass

    def _add_attribute(self, attribute : attributes.Attribute | CoordinateSystem, midx : int | None = None) -> int:
        """
        Adds attribute to the corresponding attribute list
        Assigns its midx
        An empty attribute can be add as a placeholder in the form of (id, None)
        Returns attribute id
        """
 
        att_list = None
        if isinstance(attribute, attributes.Mat):
            att_list = self._material_list
        elif isinstance(attribute, attributes.Etype):
            att_list = self._element_type_list
        elif isinstance(attribute, attributes.Real):
            att_list = self._real_constants_list
        elif isinstance(attribute, CoordinateSystem):
            att_list = self._cs_list
        elif isinstance(attribute, attributes.Section):
            att_list = self._section_list

        def find_max_attribute_id() -> int:
            if not att_list:
                """if isinstance(attribute, CoordinateSystem):
                # Skip default attribute ids for mapdl
                    return 0 
                else: return 1"""
                return 0
            else: return max(att.midx for att in att_list)

        if midx is not None:
            attribute.midx = midx
        else:
            attribute.midx = find_max_attribute_id() + 1

        att_list.append(attribute)
        return attribute.midx

#####################
## SETS
#####################

    def add_node_set(self, name : str, nodes : np.ndarray | list | None, **kwargs) -> Nset:
        if name is None:
            raise ParameterError("Node set name must be defined!")
        
        nset = Nset(name, nodes)

        nset.midx = self._find_max_tree_item_id(self._model_definitions["Sets"]["Node Sets"]) + 1
        self._model_definitions["Sets"]["Node Sets"][nset.name] = nset

        if kwargs and kwargs["ignore_at_execution"]: # ignore_at_execution
            nset.ignore_at_execution = kwargs["ignore_at_execution"]

        return nset

    def add_element_set(self, name : str, elements : np.ndarray | list, **kwargs) -> Elset:
        """Adds Element Set to the Model & Model Tree"""
        if name == '':
            raise ParameterError("Element set name must be defined!")
        
        elset = Elset(name, elements)

        #self._element_set_list.append(elset)
        self._model_definitions["Sets"]["Element Sets"][elset.name] = elset

        if kwargs and kwargs["ignore_at_execution"]: # ignore_at_execution
            elset.ignore_at_execution = kwargs["ignore_at_execution"]

        return elset

    def add_pilot_node(self, name : str, x : float, y : float, z : float, local_cs : CoordinateSystem | None = None) -> PilotNode:

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

    def _add_surface(self, nset : str, 
                     name : str | None = None, 
                     surface_type : Literal["NODE SET", "ELEMENT FACES", "SURFACE ELEMENTS"] = "ELEMENT FACES",
                     etype : attributes.Etype | None = None,
                     real_constants : attributes.Real | None = None,
                     material : attributes.Mat | None = None,
                     coordinate_system : CoordinateSystem | None = None,
                     **kwargs) -> Surface:
        """Create a surface object through given Nset, return surface reference"""
        
        surf = Surface(self.get_nset(nset), name, surface_type)
        self._surface_list.append(surf)

        if surface_type == "SURFACE ELEMENTS":

            self._add_attribute(etype)
            self._add_attribute(real_constants)
            self._add_attribute(material)
            
            surf.element_type = etype
            surf.real_constants = real_constants
            surf.material = material
            surf.csys = coordinate_system

        if kwargs and kwargs["ignore_at_execution"]: # ignore_at_execution
            surf.ignore_at_execution = kwargs["ignore_at_execution"]

        return surf

    def _get_surface(self, nset : str, type : Literal["NODE SET", "ELEMENT FACES", "SURFACE ELEMENTS"]) -> Surface | None:

        _nset = self.get_nset(nset)
        if _nset:
            for surf in self._surface_list:
                if set(surf.nset.items) == set(_nset.items) and type == surf.type:
                    return surf
                else: return None
        else: ParameterError(f"Node set '{_nset}' does not exist in the model!")

        """if _nset in self._surface_list:
            return self._model_definitions["Sets"]["Surfaces"][_nset]
        else: return None"""

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

        nset_obj = self.get_nset(nset)

        # create a dof-coupling obj
        coupling_obj = constraint.CoupledDOF(nset=nset_obj, dof=dof)

        if "couplings" in self._model_definitions:
            id = self._find_max_tree_item_id("couplings") + 1
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

    def add_force_dist_surf_constraint(self, 
                                       pilot_node : str, contact_nodes : str,
                                       dof : tuple = ("UX", "UY", "UZ", "ROTX", "ROTY", "ROTZ"),
                                       weight : Literal["Contact Area", "RBE3"] | np.ndarray = "Contact Area",
                                       contact_model : Literal["Surface-To-Surface", "Node-To-Surface"] = "Node-To-Surface",
                                       contact_algorithm : Literal["MPC", "Lagrange"] = "MPC",
                                       local_cs : CoordinateSystem | None = None,
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
        cs : CoordinateSystem
        if not local_cs:
            cs = self.get_coordinate_system("Cartesian")
        else: cs = local_cs

        mpc = constraint.ForceDistributedConstraint(pilot_node=pilot_node, contact_nodes=contact_nodes,
                                                    constrained_dof=dof,
                                                    contact_model=contact_model,
                                                    contact_algorithm=contact_algorithm,
                                                    coordinate_system=cs,
                                                    bonding_type=bonding_type,
                                                    constrained_surface_symmetry=symmetry)

        # Add target and contact element types to the model
        self._add_attribute(mpc.contact_element_type)
        self._add_attribute(mpc.target_element_type)
        # Create an empty material attribute
        mpc.material = attributes.Mat(ignore_at_execution=True)
        self._add_attribute(mpc.material)
        # Create an empty Real Constant attribute
        mpc.real_constants = attributes.Real(ignore_at_execution=True)
        self._add_attribute(mpc.real_constants)

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
        
        mpc.midx = self._find_max_tree_item_id(self._model_definitions["Constraints"]["Surface-Based"]) + 1
        self._model_definitions["Constraints"]["Surface-Based"][f"Force-Dist Surface-{mpc.midx} ({contact_algorithm})"] = mpc
        return mpc

    def add_rigid_surface_constraint(self, pilot_node : str, contact_nodes : str,
                     dof : tuple = ("UX", "UY", "UZ", "ROTX", "ROTY", "ROTZ"),
                     contact_model : Literal["Surface-To-Surface", "Node-To-Surface"] = "Node-To-Surface",
                     contact_algorithm : Literal["MPC", "Lagrange"] = "MPC",
                     local_cs : CoordinateSystem | None = None,
                     bonding_type : Literal["Bonded (always)", "Bonded (initial)"] = "Bonded (always)") -> constraint.RigidSurfaceConstraint:

        mpc = constraint.RigidSurfaceConstraint(pilot_node=pilot_node, contact_nodes=contact_nodes,
                                             constrained_dof=dof,
                                             contact_model=contact_model,
                                             contact_algorithm=contact_algorithm,
                                             local_cs=local_cs,
                                             bonding_type=bonding_type)
        
        # Add internally generated contact elements to the model
        mpc.contact_element_type.midx = self._add_attribute(mpc.contact_element_type)
        mpc.target_element_type.midx = self._add_attribute(mpc.target_element_type)

        # Create an empty material attribute
        mpc.material_midx = self._add_attribute(attributes.Mat())

        # Create and Associate the Real Constant Object
        mpc.real_constants = element.ContaRealConstants()
        mpc.real_constants.midx = self._add_attribute(mpc.real_constants)
        
        mpc.midx = self._find_max_tree_item_id(self._model_definitions["Constraints"]["Surface-Based"]) + 1
        self._model_definitions["Constraints"]["Surface-Based"][f"Rigid Surface-{mpc.midx} ({contact_algorithm})"] = mpc
        return mpc

    def add_coupled_surface_constraint(self, pilot_node : str, contact_nodes : str,
                     dof : tuple = ("UX", "UY", "UZ", "ROTX", "ROTY", "ROTZ"),
                     contact_model : Literal["Surface-To-Surface", "Node-To-Surface"] = "Node-To-Surface",
                     local_cs : CoordinateSystem | None = None,
                     bonding_type : Literal["Bonded (always)", "Bonded (initial)"] = "Bonded (always)") -> constraint.CoupledSurfaceConstraint:

        mpc = constraint.CoupledSurfaceConstraint(pilot_node=pilot_node, contact_nodes=contact_nodes,
                                             constrained_dof=dof,
                                             contact_model=contact_model,
                                             local_cs=local_cs,
                                             bonding_type=bonding_type)
        
        # Add internally generated contact elements to the model
        mpc.contact_element_type.midx = self._add_attribute(mpc.contact_element_type)
        mpc.target_element_type.midx = self._add_attribute(mpc.target_element_type)

        # Create an empty material attribute
        mpc.material_midx = self._add_attribute(attributes.Mat())

        # Create and Associate the Real Constant Object
        mpc.real_constants = element.ContaRealConstants()
        mpc.real_constants.midx = self._add_attribute(mpc.real_constants)
        
        mpc.midx = self._find_max_tree_item_id(self._model_definitions["Constraints"]["Surface-Based"]) + 1
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
                     analysis : Literal["STATIC", "BUCKLE", "MODAL", "HARMIC" ,"TRANS", "SUBSTR"] = "STATIC") -> LoadStep:
        
        step_number = len(self._load_step_list) + 1
        ls = LoadStep(name= name,
                      model= self,
                      step_number= step_number,
                      end_time= 1,
                      analysis_type= analysis)
        
        self._load_step_list.append(ls)
        self._model_definitions["Load Steps"][name] = ls   
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
                    if pressure.set.nset.name == nset and pressure.direction == direction:
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

    def bc_history(self) -> pd.DataFrame:
        return self._analysis.get_bc_history()

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

    def _find_max_tree_item_id(self, group : dict) -> int:
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

        for command in self._command_pipeline:

            if args[0] == False and command[0].__name__ == "_solve":
                #self._add_command("_write_comment", "------ SOLVE COMMAND SKIPPED!!! -----") 
                continue # Skip "Solve" Command

            if len(command) > 1:
                command[0](self.mapdl, *command[1:])
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

        for node in self._pilot_node_list:
            self._add_command("_create_pilot_node", node)

        for etype in self._element_type_list:
            if not etype.ignore_at_execution: self._add_command("_create_element_type", etype)

        for rc in self._real_constants_list:
            if not rc.ignore_at_execution: self._add_command("_create_real_constants", rc)

        for mat in self._material_list:
            if not mat.ignore_at_execution: self._add_command("_material", mat)

        for cs in self._cs_list:
            if not cs.ignore_at_execution: self._add_command("_local_csys", cs)

             # Creating Sections Skipped

        for nset in self._model_definitions["Sets"]["Node Sets"].values():
            if not nset.ignore_at_execution: self._add_command("_create_node_set", nset)

        for elset in self._model_definitions["Sets"]["Element Sets"].values():
            if not elset.ignore_at_execution: 
                self._add_command("_create_element_set", elset)
                self._add_command("_assign_element_type", elset)
                self._add_command("_assign_real_constants", elset)
                self._add_command("_assign_material", elset)

        for surf in self._surface_list:
            if not surf.ignore_at_execution: self._add_command("_create_surface", surf)

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

        """self._add_command("_enter_solution") # Enter once for multistep analysis !!!
        for load_step in self._load_step_list:
            self._add_command("_load_step", load_step)
            self._add_command("_solve")"""
            
    @overload
    def solve(self, status : Literal["NEW", "RESTART"] = "NEW"): ...
    @overload
    def solve(self, status : Literal["NEW", "RESTART"] = "RESTART", 
              step : int = 1, substep : int = 1): ...
    def solve(self,
              status : Literal["NEW", "RESTART"] = "NEW",
              step : int = 0, substep : int = 0):
        
        

        if status == "RESTART":

            # !!! Be careful that for _load_step_list[step] indicate loadstep-2 for step parameter 1 !!!
            self._load_step_list[step].status = "RESTART"
            self._load_step_list[step].substep = substep          

            timestamp = datetime.now().strftime("%d.%m.%Y-%H.%M.%S")
            folder_name = f"Restart-LS{self._load_step_list[step].step_number}-{timestamp}"
            self._load_step_list[step].filename = f'Restart-{self._load_step_list[step].name}'
            self.target_dir = self.working_dir / folder_name
            self.target_dir.mkdir(exist_ok=True)

            self._add_command("_enter_solution") # Enter once for multistep analysis !!!
            self._add_command("_restart", step, substep)
            self._add_command("_solve")
            self._add_command("_change_job_name", self._load_step_list[step].filename)
            self._add_command("_enter_solution") # Enter once for multistep analysis !!!
            for load_step in self._load_step_list[step:]:
                self._add_command("_load_step", load_step)
                self._add_command("_solve")
            self._add_command("_exit_mapdl")

            # Move files
            if self.target_dir != self.working_dir: # This is not equal in case of a restart analysis
                self._move_solution_files(self._load_step_list[step].filename, self.target_dir)
        else:
            self.working_dir = self._current_path / self.folder_name
            self.working_dir.mkdir()
            self.target_dir = self.working_dir # For possible restart files movement

            self._create_command_pipeline()
            self._add_command("_enter_solution") # Enter once for multistep analysis !!!
            for load_step in self._load_step_list:
                self._add_command("_load_step", load_step)
                self._add_command("_solve")
            self._add_command("_exit_mapdl")
        
        self._execute_commands(True)
        self._command_pipeline.clear()


        

    def _move_solution_files(self, filename, target_dir):
        from pathlib import Path
        import shutil

        source_dir = Path(self.working_dir)
        target_dir.mkdir(exist_ok=True)

        for file in source_dir.glob(f"{filename}*"):
            shutil.move(str(file), target_dir / file.name)

    def write_input_file(self):
        self._create_command_pipeline()
        self._add_command("_write_input_file")
        self._add_command("_exit_mapdl")
        self._execute_commands(False) # Skip Solve Commands
        self._command_pipeline.clear()

    def _start(self):

        #try:
        self.mapdl = launch_mapdl(run_location=self.working_dir,
                                  jobname=self.modelname,
                                  nproc=self.cpus,
                                  override=True)
        #except: self.mapdl.exit()
        
        # First, reset the MAPDL database.
        self.mapdl.clear()
        return self.mapdl

    def solution(self) -> Any:
        """Returns DPF solution object after solving the model."""
        _result_file_path = Path(self.working_dir / f"{self.modelname}.rst")

        solution = None
        if self.solver_name == "MAPDL":
            from structuralanalysistoolbox.post.mapdl import solutions
            solution = solutions.DpfSolution(_result_file_path)

        return solution.simulation








    





    
