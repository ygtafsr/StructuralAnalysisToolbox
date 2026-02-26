
"""
GENERAL CONSIDERATIONS ABOUT STBOX ELEMENT DEFINITIONS & IMPLEMENTATIONS

# Sometimes element types, which are defined internally as a result of another implementations
  such as MPC surface-based constraints, contact elements, are not put into model tree. Their
  properties (keyopts, real constants etc) can be editted within actual implementation object.

# Element types object attributes with a default value as 'None' skipped implementing keyopts and
  real constants.

"""




from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum
from typing import Literal, overload
from abc import abstractmethod
import numpy as np

#from structuralanalysistoolbox.mapdl.mesh import ElementType
from structuralanalysistoolbox.mapdl import element, attributes

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from structuralanalysistoolbox import model
    from structuralanalysistoolbox.materials.material import Material
    

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
    independent : model.Nset
    dof : str

###########################
## Multi Point Constraints
###########################

@dataclass
class MPC:

    etype : element.MPC184 | None = None
    midx : int = 0   # model index
    dependent : model.Nset | None = None
    independent : model.Nset | None = None  

@dataclass
class MPCRigidLink(MPC):
    
    behaviour = 0
    method : Literal["Direct Elimination", "Lagrange Multiplier"] = "Lagrange Multiplier"
    contribution : Literal["Include", "Exlude"] = "Include"
    
@dataclass
class MPCRigidBeam(MPC):
    
    behaviour = 1
    method : Literal["Direct Elimination", "Lagrange Multiplier"] = "Lagrange Multiplier"
    contribution : Literal["Include", "Exlude"] = "Include"
    
@dataclass
class MPCSlider(MPC):
    
    behaviour = 3
    method : Literal["Lagrange Multiplier", "Penalty-based"] = "Lagrange Multiplier"
    
##################
### MPC184 JOINTS
##################
"""
    MPC184:
    https://ansyshelp.ansys.com/public/account/secured?returnurl=//Views/Secured/corp/v242/en/ans_elem/Hlp_E_MPC184.html?q=MPC184

    SECTYPE:
    https://ansyshelp.ansys.com/public/account/secured?returnurl=///Views/Secured/corp/v242/en/ans_cmd/Hlp_C_SECTYPE.html

    SECDATA:
    https://ansyshelp.ansys.com/public/account/secured?returnurl=////Views/Secured/corp/v242/en/ans_cmd/Hlp_C_SECDATA.html

    SECJOINT DATA:
    https://ansyshelp.ansys.com/public/account/secured?returnurl=//Views/Secured/corp/v242/en/ans_cmd/Hlp_C_SECJOINT.html
"""
@dataclass
class MPCJoint(MPC):

    material : str = None                           # Join material option on the TB command
    penalty_factors : str = None                    # CHECK SECJOINT DATA
    bc : str = None                                 # DJ or FJ boundary conditions
    local_cs : model.CoordinateSystem = None   # SECJOINT DATA
    stop : list = None                              # SECSTOP (Stops or limits)
    lock : list = None                              # SECLOCK
    method : Literal["Lagrange Multiplier", "Penalty-Based"] = "Lagrange Multiplier"
    sectiontype = "JOINT"

@dataclass
class RevoluteXJoint(MPCJoint):

    angle : float =  None    # SECDATA
    friction : list  = None # SECJOINT DATA  :: friction: (Outer radius, Inner radius, Effective length)
    configuration = "x-axis"
    behaviour = 6
    section_sub_type = "REVO"

@dataclass
class RevoluteZJoint(MPCJoint):

    angle : float = None    # SECDATA
    friction : list = None  # SECJOINT DATA :: friction: (Outer radius, Inner radius, Effective length)
    configuration = "z-axis"
    behaviour = 6
    section_sub_type = "REVO"

@dataclass
class UniversalJoint(MPCJoint):

    angle_1 : float = None # SECDATA
    angle_2 : float = None # SECDATA
    behaviour = 7
    section_sub_type = "UNIV"

@dataclass
class SlotJoint(MPCJoint):

    length : float = None   # SECDATA
    friction : float = None # SECJOINT DATA :: friction: 
    behaviour = 8
    section_sub_type = "SLOT"

@dataclass
class PointJoint(MPCJoint):

    lenght_1 : float = None # SECDATA
    length_2 : float = None # SECDATA
    behaviour = 9
    section_sub_type = "PINP"

@dataclass
class TranslationalJoint(MPCJoint):

    length : float = None # SECDATA
    friction : list = None # SECJOINT DATA :: friction: (Effective length, Effective radius)
    behaviour = 10
    section_sub_type = "PRIS"

@dataclass
class CylindricalXJoint(MPCJoint):

    length : float = None # SECDATA
    angle : float = None # SECDATA
    configuration = "x-axis"
    behaviour = 11
    section_sub_type = "CYLI"

@dataclass
class CylindricalZJoint(MPCJoint):

    length_3 : float = None # SECDATA
    angle_3 : float = None # SECDATA
    configuration = "z-axis"
    behaviour = 11
    section_sub_type = "CYLI"

@dataclass
class PlanarXJoint(MPCJoint):

    lenght_2 : float = None # SECDATA
    length_3 : float = None # SECDATA
    angle_1 : float = None # SECDATA
    configuration = "x-axis"
    behaviour = 12
    section_sub_type = "PLAN"

@dataclass
class PlanarYJoint(MPCJoint):
    lenght_1 : float = None # SECDATA
    length_2 : float = None # SECDATA
    angle_3 : float = None # SECDATA
    configuration = "x-axis"
    behaviour = 12
    section_sub_type = "PLAN"

@dataclass
class WeldJoint(MPCJoint):
    behaviour = 13
    section_sub_type = "WELD"

@dataclass
class OrientJoint(MPCJoint):
    behaviour = 14
    section_sub_type = "ORIE"

@dataclass
class SphericalJoint(MPCJoint):
    angle_1 : float = None # SECDATA
    angle_2 : float = None # SECDATA
    angle_3 : float = None # SECDATA
    friction : int = None  # SECJOINT DATA :: friction: (Effective radius)
    behaviour = 15
    section_sub_type = "SPHE"

@dataclass
class GeneralJoint(MPCJoint):
    lenght_1 : float = None # SECDATA
    lenght_2 : float = None # SECDATA
    lenght_3 : float = None # SECDATA
    angle_1 : float = None # SECDATA
    angle_2 : float = None # SECDATA
    angle_3 : float = None # SECDATA
    configuration : Literal["displacement and rotation",
                            "only displacement"] = None
    relative_dof : Literal["UX", "UY", "UZ", "ROTX", "ROTY", "ROTZ"] = None # SECJOINT DATA
    behaviour = 16
    section_sub_type = "GENE"

@dataclass
class ScrewJoint(MPCJoint):
    lenght_3 : float = None # SECDATA
    angle_3 : float = None # SECDATA
    #The pitch is defined as the ratio of relative axial displacement (length units) to relative rotation (in radians).
    pitch : int = None # SECJOINT DATA
    behaviour = 17
    section_sub_type = "SCRE"

@dataclass
class SpotWeldJoint(MPCJoint):
    location : float = None # spring-damper location # SECJOINT DATA
    behaviour = 18
    section_sub_type = "SPWE"

@dataclass
class GenbJoint(MPCJoint):
    lenght_1 : float = None # SECDATA
    lenght_2 : float = None # SECDATA
    lenght_3 : float = None # SECDATA
    angle_1 : float = None # SECDATA
    angle_2 : float = None # SECDATA
    angle_3 : float = None # SECDATA
    configuration : Literal["displacement and rotation",
                            "only displacement"] = None
    relative_dof : Literal["UX", "UY", "UZ", "ROTX", "ROTY", "ROTZ"] = None # SECJOINT DATA
    behaviour = 19
    section_sub_type = "GENB"


##############################
## Surface Based Constraints
##############################

@dataclass
class RigidSurfaceConstraint:
    """
    Rigid Surface Constraint
    Constrained Surface Remains Rigid
    (Similar to CERIG)

    * Rigid surface constraint can be based on either MPC approach (keopt2 = 2)
      or Lagrange multiplier method (keyopt2 = 3)

    * The pilot node is an independent(retained) node. The contact nodes are the
      dependent (removed) nodes. Do not apply additional constraints to the dependent nodes.
    """  
    pilot_node : str

    contact_nodes : str

    midx : int = 0   # model index

    constrained_dof : tuple = ("UX", "UY", "UZ", "ROTX", "ROTY", "ROTZ")

    contact_algorithm : Literal["MPC", "Lagrange & Penalty"] = "MPC"

    contact_model : Literal["Surface-To-Surface",
                            "Node-To-Surface"] = "Node-To-Surface"

    local_cs : model.CoordinateSystem | None = None

    bonding_type : Literal["Bonded (always)", "Bonded (initial)"] = "Bonded (always)"

    material_midx : int = 0

    contact_element_type : element.Conta172 | element.Conta174 | element.Conta175 | element.Conta177 | None = None

    target_element_type : element.Targe169 | element.Targe170 | None = None

    real_constants : element.ContaRealConstants | None = None

    def __post_init__(self):
        """Create underlying contact and target element objects"""
        # Post-Correction :: 
        if self.contact_algorithm == "Lagrange": self.contact_algorithm = "Lagrange & Penalty"

        if self.contact_model == "Node-To-Surface":
            self.contact_element_type = element.Conta175(contact_algorithm=self.contact_algorithm,
                                                        surface_based_constraint = "Rigid Surface",
                                                        behaviour=self.bonding_type)
        
            self.target_element_type = element.Targe170(dof=self.constrained_dof, 
                                                        bc_for_rigid_target_nodes="User Definition")
        elif self.contact_model == "Surface-To-Surface":
            pass

@dataclass
class ForceDistributedConstraint:
    """
      Force-Distributed Constraint (Similar to RBE3)

    * Force-distributed surface constraint can be based on either MPC approach (keopt2 = 2)
      or Lagrange multiplier method (keyopt2 = 3) 

    * Pilot node is a dependent node (meaning the degrees of freedom for this node are removed)
      The contact nodes are independent nodes (the degrees of freedom are retained).
      If the pilot node has constraints applied to it, internally-generated MPC equations are
      rewritten so that the degrees of freedom of the pilot node are no longer dependant DOF.

    * By default, the program computes weighting factors automatically by summing the contact area
      of each contact node.
      !!! When you set keyopt7 = 1 on the target element (targe169 and targe170), the program uses
      a constant weighting factor of 1.0 for each contact node. (Similar to default behaviour of RBE3)
    """  
    pilot_node : str

    contact_nodes : str

    midx : int = 0   # model index

    constrained_dof : tuple = ("UX", "UY", "UZ", "ROTX", "ROTY", "ROTZ")

    contact_algorithm : Literal["MPC", "Lagrange & Penalty"] = "MPC"

    contact_model : Literal["Surface-To-Surface",
                            "Node-To-Surface"] = "Node-To-Surface"

    weighting_factor : np.ndarray | None = None 

    coordinate_system : model.CoordinateSystem | None = None

    bonding_type : Literal["Bonded (always)", "Bonded (initial)"] = "Bonded (always)"

    material : attributes.Mat | None = None

    constrained_surface_symmetry : Literal["Pilot XY Plane",
                                            "Pilot XZ Plane",
                                            "Pilot YZ Plane",
                                            "Pilot XY + XZ",
                                            "Pilot XY + YZ",
                                            "Pilot XZ + YZ"] | None = None 
    
    """stress_stiffening_effect : Literal["Include", "Not Include"] = "Include" # Target keyopt10 only if contact keyopt2 = mpc

    relaxion_method : Literal["Include", "Not Include"] = "Include" # Target keyopt11 only if contact keyopt2 = mpc"""

    contact_element_type : element.Conta172 | element.Conta174 | element.Conta175 | element.Conta177 | None = None

    target_element_type : element.Targe169 | element.Targe170 | None = None

    real_constants : attributes.Real | None = None

    """@property
    def stress_stiffening_effect(self):
        return self._stress_stiffening_effect
    
    @stress_stiffening_effect.setter
    def stress_stiffening_effect(self, value : Literal["Include", "Not Include"] = "Include"):"""


    def __post_init__(self):
        """Create underlying contact and target element objects"""
        # Post-Correction :: 
        if self.contact_algorithm == "Lagrange": self.contact_algorithm = "Lagrange & Penalty"

        if self.contact_model == "Node-To-Surface":
            self.contact_element_type = element.Conta175(contact_algorithm=self.contact_algorithm,
                                                         surface_based_constraint = "Force Distribution",
                                                         behaviour=self.bonding_type)
        
            self.target_element_type = element.Targe170(dof=self.constrained_dof, 
                                                        constrained_surface_symmetry=self.constrained_surface_symmetry,
                                                        bc_for_rigid_target_nodes="User Definition")
        elif self.contact_model == "Surface-To-Surface":
            pass

        

@dataclass
class CoupledSurfaceConstraint:
    """
    Coupling Constraint
    (Similar to CP)
    """
    pilot_node : str

    contact_nodes : str

    midx : int = 0   # model index

    constrained_dof : tuple = ("UX", "UY", "UZ", "ROTX", "ROTY", "ROTZ")

    contact_model : Literal["Surface-To-Surface",
                            "Node-To-Surface"] = "Node-To-Surface"

    local_cs : model.CoordinateSystem | None = None

    bonding_type : Literal["Bonded (always)", "Bonded (initial)"] = "Bonded (always)"

    material_midx : int = 0

    contact_element_type : element.Conta172 | element.Conta174 | element.Conta175 | element.Conta177 | None = None

    target_element_type : element.Targe169 | element.Targe170 | None = None

    real_constants : element.ContaRealConstants | None = None

    def __post_init__(self):
        """Create underlying contact and target element objects"""

        if self.contact_model == "Node-To-Surface":
            self.contact_element_type = element.Conta175(contact_algorithm= "MPC",
                                                        surface_based_constraint = "Coupling",
                                                        behaviour=self.bonding_type)
        
            self.target_element_type = element.Targe170(dof=self.constrained_dof, 
                                                        bc_for_rigid_target_nodes="User Definition")
        elif self.contact_model == "Surface-To-Surface":
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