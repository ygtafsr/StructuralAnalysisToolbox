
from __future__ import annotations
from dataclasses import dataclass
from enum import Enum

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from structuralanalysistoolbox import model

class Coupling(Enum):
    DOF = 1
    Interface = 2
    Periodic = 3

@dataclass
class CouplingDOF:
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
    id : int 
    nset : model.Nset
    dof : str  

class CouplingInterface:

    id : int
    type : str  # Bonded/Sliding
    nset : str
    tolerance : float
                
class CouplingPeriodic:
    """Defines a Periodic Boundary Condition"""
    pass

@dataclass
class Constraint:
    pass

@dataclass
class MultiPointConstraint:
    pass

@dataclass
class Spring:
    pass

@dataclass
class Beam:
    pass