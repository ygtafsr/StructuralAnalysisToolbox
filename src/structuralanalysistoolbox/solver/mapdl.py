
from __future__ import annotations
import functools
import ansys.mapdl.core as core
import numpy as np
from enum import Enum
from dataclasses import asdict

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from structuralanalysistoolbox.materials import material
    from structuralanalysistoolbox import model
    from structuralanalysistoolbox.constraints import constraint

def prep7(func):
    @functools.wraps(func)
    def wrapper_prep7(*args, **kwargs):  
        args[0].prep7()
        val = func(*args, **kwargs)
        args[0].finish()
        return val
    return wrapper_prep7


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

@prep7
def _coupled_dof(mapdl : core.Mapdl, coupling : constraint.CoupledDOF) -> int:
    """Create a dof coupling between nodes and returns primary node number 
    as an integer for coupled set."""
    
    if coupling.dof == "ALL":
        # Using "ALL" directly results new ids for each dof. 
        # So, this method prefered.
        mapdl.cp(nset=f"{coupling.id}", lab="UX", node1=coupling.nset.name)
        mapdl.cp(nset=f"{coupling.id}", lab="UY", node1=coupling.nset.name)
        mapdl.cp(nset=f"{coupling.id}", lab="UZ", node1=coupling.nset.name)
        mapdl.cp(nset=f"{coupling.id}", lab="ROTX", node1=coupling.nset.name)
        mapdl.cp(nset=f"{coupling.id}", lab="ROTY", node1=coupling.nset.name)
        mapdl.cp(nset=f"{coupling.id}", lab="ROTZ", node1=coupling.nset.name)
    elif isinstance(coupling.dof, tuple):
        for dof in coupling.dof:
            mapdl.cp(nset=f"{coupling.id}", lab=dof, node1=coupling.nset.name)
    else:
        mapdl.cp(nset=f"{coupling.id}", lab=coupling.dof, node1=coupling.nset.name)
    #primary_node_id = _get(mapdl, coupling.nset, Get.Min_Node_Number)
    #return primary_node_id

@prep7
def _coupled_interface(mapdl : core.Mapdl, coupling : constraint.CoupledInterface):
    pass

@prep7                           
def _MPC184(mapdl : core.Mapdl, mpc: constraint.MPC) -> None:

    mapdl.et(itype=mpc._etype.id,
             ename=mpc._etype.type.value,
             kop1=mpc._etype.get_keyop(1),
             kop2=mpc._etype.get_keyop(2)
             )
    
    mapdl.type(mpc._etype.id)

    for node in mpc.dependent:
        mapdl.e(node, mpc.independent)





