
from __future__ import annotations
import ansys.mapdl.core as core
import numpy as np
from enum import Enum
from dataclasses import asdict

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from structuralanalysistoolbox.materials import material
    from structuralanalysistoolbox import model
    from structuralanalysistoolbox.constraints import constraint

"""from materials import material"""


def _material(mapdl : core.Mapdl, mat : material.Material, mat_id : int) -> None:
    """Runs mapdl material definition commands"""

    mapdl.prep7()
    
    for model in mat.material_models.values():

        if isinstance(model, material.Physical):
            # Density
            if isinstance(model.density, (float, int)):
                mapdl.mp(lab="DENS", mat=f"{mat_id}", c0=f"{model.density}")
            else:  
                label = "density"
                mapdl.load_table(name=f"{label}_{mat_id}", array=np.array(model.density), var1="TEMP")
                mapdl.mp(lab="DENS", mat=f"{mat_id}", c0=f"%{label}_{mat_id}%")
        elif isinstance(model, material.IsotropicElasticity):
            # Elastic Modulus
            if isinstance(model.elastic_modulus, (float, int)):
                mapdl.mp(lab="EX", mat=f"{mat_id}", c0=f"{model.elastic_modulus}")
            else:
                label = "elastic_modulus"
                mapdl.load_table(name=f"{label}_{mat_id}", array=np.array(model.elastic_modulus), var1="TEMP")
                mapdl.mp(lab="EX", mat=f"{mat_id}", c0=f"%{label}_{mat_id}%")
            # Poisson's Ratio
            if isinstance(model.poisson_ratio, (float, int)):
                mapdl.mp(lab="NUXY", mat=f"{mat_id}", c0=f"{model.poisson_ratio}")
            else:
                label = "poissons_ratio"
                mapdl.load_table(name=f"{label}_{mat_id}", array=np.array(model.poisson_ratio), var1="TEMP")
                mapdl.mp(lab="NUXY", mat=f"{mat_id}", c0=f"%{label}_{mat_id}%")
        elif isinstance(model, material.MultilinearIsotropicHardening):
            # Multilinear Isotropic Hardening
            if model.temperature == None:
                data_count = len(model.data_table)
                mapdl.tb(lab="PLASTIC", mat=f"{mat_id}", ntemp="1", 
                                                            npts=f"{data_count}")
                mapdl.tbtemp(temp="")
                for data in model.data_table:
                    mapdl.tbpt(oper="DEFI", x1=f"{data[0]}", x2=f"{data[1]}")
            elif isinstance(model.temperature, (float, int)):
                data_count = len(model.data_table)
                mapdl.tb(lab="PLASTIC", mat=f"{mat_id}", ntemp="1", 
                                                            npts=f"{data_count}")
                mapdl.tbtemp(temp=f"{model.temperature}")
                for data in model.data_table:
                    mapdl.tbpt(oper="DEFI", x1=f"{data[0]}", x2=f"{data[1]}")

    mapdl.finish()

############################
## GET 
############################ 

# GET minimum node number for a given node set
# GET number of nodes currently selected
# GET next number after a given one in a set ??

class Get(Enum):

    Max_Node_Number = 1
    Min_Node_Number = 2


def _get(mapdl : core.Mapdl, set : model.Component, value : Get):
    # *GET, Par, Entity, Item1, IT1NUM, Item2, ITNUM2
    # mapdl.get(par="", entity="", item1="", it1num="", ...) -> (float, str)

    # Select items
    if type(set).__name__ == "Nset": # Type check without referencing actual module
        mapdl.cmsel(type_="S", name=set.name, entity="NODE")
        if value == Get.Max_Node_Number:
            num = mapdl.get(par="NMAX", entity="NODE", item1="NUM", it1num="MAX")
            mapdl.allsel()  # Reset selection before return
            return int(num)
        elif value == Get.Min_Node_Number:
            num = mapdl.get(par="NMAX", entity="NODE", item1="NUM", it1num="MIN")
            mapdl.allsel()  # Reset selection before return
            return int(num)

############################
## CONSTRAINTS
############################

def _coupling_dof(mapdl : core.Mapdl, coupling : constraint.CouplingDOF) -> int:
    """Create a dof coupling between nodes and returns primary node number 
    as an integer for coupled set."""
    mapdl.prep7()

    mapdl.cp(nset=f"{coupling.id}", lab=coupling.dof, node1=coupling.nset.name)
    
    mapdl.finish()
    primary_node_id = _get(mapdl, coupling.nset, Get.Min_Node_Number)
    return primary_node_id

def _coupling_interface(mapdl : core.Mapdl, id : int, nset : model.Nset, tolerance : float):
    pass

                            
def _create_MPC184_rigid(mapdl : core.Mapdl, dependent_nodes : str, independent_node : str) -> None:
    
    # Check current processor ???
    # Check whether components are exist ???
    
    mapdl.prep7()

    # get the components content
    _dependent_nodes = mapdl.components[dependent_nodes].items
    _independent_node = mapdl.components[independent_node].items[0]

    ## ????
    _current_max_type_number = int(mapdl.get('max_type_number','ETYP','','NUM','MAX'))
 
    # create the MPC184 elements
    mapdl.et(_current_max_type_number+1,'MPC184',1,1)
    mapdl.type(_current_max_type_number+1)

    for node in _dependent_nodes:
        mapdl.e(node, _independent_node)

    #mapdl.allsel()  

    mapdl.finish()                 