
from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum

from structuralanalysistoolbox.meshing.mesh import ElementType, Etype

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from structuralanalysistoolbox import model
    

class Coupling(Enum):
    DOF = 1
    Interface = 2
    Periodic = 3

# COMMANDS:
# CP
# CPINTF vs CEINTF
# CPCYC
# NUMMRG, MERGE
# EINTF
# 

########################
## Couplings
########################  

@dataclass
class CoupledDOF:
    """
    nset: nodal component name
    The first node in the set is the primary node.
    id: A unique integer represents coupling
    dof (lab in MAPDL syntax):
    # UX, UY, UZ, ROTX, ROTY, ROTZ
    # HDSP (Hydrostatic pressure)
    # ALL
    When dof is "ALL", sets are generated for each
    active dof, and NSET is incremented auto to prevent 
    overwriting existing sets. 
    """
    nset : model.Nset
    dof : str  
    id : int = 0


@dataclass
class CoupledInterface:
    """
    Couple Coincident Nodes

    CPINTF: couples coincident nodes in a model by generating 
    one coupled set for each specified dof label at every pair of
    coincident nodes. This operation is useful for "buttoning"
    together several pair of nodes (such as at a seam).

    • If all degrees of freedom are to be coupled for coincident nodes, it is usually more efficient to
      simply merge those nodes together via NUMMRG. ++

    • You can connect coincident pairs of nodes by creating two-node elements between them via
      EINTF.

    • To tie together two regions having dissimilar mesh patterns, issue CEINTF. This operation
      generates constraint equations that connect the selected nodes of one region to the selected
      elements of the other region.
    """
    type : str  # Bonded(Coincident)/Sliding/Two-Node Element
    nset : model.Nset
    tolerance : float = 0.001
                
@dataclass
class CoupledPeriodic:
    """Defines a Periodic Boundary Condition"""
    type : str  # Bonded(Coincident)/Sliding/Two-Node Element
    nset : model.Nset
    pass

########################
## Constraint Equations
########################

@dataclass
class Tie:
    """
    CEINTF
    """
    pass

@dataclass
class Rigid:
    # Defines a rigid region.
    # CERIG
    dependent : model.Nset
    independent : model.Node
    dof : str

########################
## MPC184
########################
"""
To define various contact assemblies and kinematic constraints, use the internal multipoint constraint
(MPC) approach (KEYOPT(2) = 2) with certain bonded and no-separation contact definitions (KEYOPT(12)
= 4, 5, or 6).
Supported contact elements are CONTA172, CONTA174, CONTA175, and CONTA177.

MAIN COMMANDS:
CERIG
RBE3

"""


MPCtypes = { "RigidLink" : 0,
            "RigidBeam" : 1,
            "Slider" : 3,
            "Revolute" : 6,
            "Universal" : 7,
            "Slot" : 8,
            "Point" : 9,
            "Translational" : 10,
            "Cylindrical" : 11,
            "Planar" : 12,
            "Weld" : 13,
            "Orient" : 14,
            "Spherical" : 15,
            "General" : 16,
            "Screw" : 17}

MPCmethods = {"DirectElemination" : 0,
             "LagrangeMultiplication" : 1,
             "SurfaceBased" : 2}
               

@dataclass
class MPC:

    name : str
    id : int  # MPC object id
    type : str
    dependent : model.Nset
    independent : model.Nset
    etype_id : int  # MPC184 Element Type id
    method : None | str = None # Default method
    axis : None | int = None
    stress_stiffening: None | int = None
    _etype : ElementType = None
    
    def __post_init__(self):

        method = 0 if self.method == None else MPCmethods[self.method]

        self._etype = ElementType(id=self.etype_id, 
                                  type=Etype.MPC184, 
                                  keyopt=[
                                         (1,MPCtypes[self.type]),
                                         (2,method),
                                         (4,self.axis),
                                         (5,self.stress_stiffening)])

@dataclass
class MPCRigidSurface:
    """
    Rigid Surface Constraint
    (Similar to CERIG/RBE2)
    """
    dependent : model.Nset
    independent : model.Nset

@dataclass
class MPCDistributed:
    """
    Force-Distributed Constraint
    (Similar to RBE3)
    """   

@dataclass
class MPCCoupling:
    """
    Coupling Constraint
    (Similar to CP)
    """
    pass

@dataclass
class MPCRigidBody:
    pass


######################################
## Connections with Two Point Elements
######################################

# EINTF !!


@dataclass
class Spring:
    pass

@dataclass
class Beam:
    pass