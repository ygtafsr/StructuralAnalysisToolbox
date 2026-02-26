
from structuralanalysistoolbox.mapdl.attributes import Etype, Real, Section, Mat

from dataclasses import dataclass, fields
from typing import Literal


SOLID185_KEYOPT_MAPPING = {
    2 : ("technology", ("FULL INTEGRATION", "REDUCED INTEGRATION", "ENHANCED STRAIN", "SIMPLIFIED ENHANCED STRAIN")),
    3 : ("layer", ("LAYERED", "HOMOGENEOUS")),
    6 : ("formulation", ("PURE DISPLACEMENT", "MIXED UP")),
    15 : ("pml_absorbing", ("NOT INCLUDED", "INCLUDED")),
    16 : ("steady_state", ("DISABLED", "ENABLED")),
    17 : ("extra_surface_output", ("BASIC", "EXTRA")),
}

@dataclass
class Solid185(Etype):
    """8-Node Structural Solid"""

    name : str = "SOLID185"

    technology : Literal["FULL INTEGRATION",
                         "REDUCED INTEGRATION",
                         "ENHANCED STRAIN",
                         "SIMPLIFIED ENHANCED STRAIN"] = "FULL INTEGRATION"

    formulation : Literal["PURE DISPLACEMENT",
                          "MIXED UP"] = "PURE DISPLACEMENT"
    
    layering : Literal["LAYERED",
                       "HOMOGENEOUS"] = "HOMOGENEOUS"

    pml_absorbing : Literal["NOT INCLUDED",
                            "INCLUDED"] = "NOT INCLUDED"
    
    steady_state : Literal["DISABLED",
                           "ENABLED"] = "DISABLED"      
    
    extra_surface_output : Literal["BASIC",
                                   "EXTRA"] = "BASIC"

    def get_keyopts(self) -> tuple: # ((keyopt_number, keyopt_value), ...)
        """Returns a tuple of keyopt number and corresponding value based on the current attribute settings."""

        return ((2 ,SOLID185_KEYOPT_MAPPING[2][1].index(self.technology)),
                (3, SOLID185_KEYOPT_MAPPING[3][1].index(self.layering)),
                (6, SOLID185_KEYOPT_MAPPING[6][1].index(self.formulation)),
                (15, SOLID185_KEYOPT_MAPPING[15][1].index(self.pml_absorbing)),
                (16, SOLID185_KEYOPT_MAPPING[16][1].index(self.steady_state)),
                (17, SOLID185_KEYOPT_MAPPING[17][1].index(self.extra_surface_output)))

    def set_keyopt(self, keyopt_number: int, value: str):
        """Sets the attribute corresponding to the keyopt_number to the given value."""

        # raise error if keyopt_number is invalid
        if keyopt_number not in SOLID185_KEYOPT_MAPPING:
            raise ValueError(f"Invalid keyopt number: {keyopt_number}")
        # raise error if value is invalid
        if value not in SOLID185_KEYOPT_MAPPING[keyopt_number][1]:
            raise ValueError(f"Invalid value for keyopt {keyopt_number}: {value}")  
        # set the attribute
        setattr(self, SOLID185_KEYOPT_MAPPING[keyopt_number][0], 
                      SOLID185_KEYOPT_MAPPING[keyopt_number][1][value])

    def __repr__(self):
        return super().__repr__()

@dataclass
class Solid186(Etype):
    """20-Node Structural Solid"""

    name : str = "SOLID186"

    technology : Literal["FULL INTEGRATION",
                         "REDUCED INTEGRATION"] = "REDUCED INTEGRATION"
    
    formulation : Literal["PURE DISPLACEMENT",
                          "MIXED UP"] = "PURE DISPLACEMENT"
    
    layering : Literal["LAYERED",
                       "HOMOGENEOUS"] = "HOMOGENEOUS"
    
    pml_absorbing : Literal["INCLUDED",
                            "NOT INCLUDED"] = "NOT INCLUDED"
    
    steady_state : Literal["ENABLED",
                          "DISABLED"] = "DISABLED"
    
    extra_surface_output : Literal["BASIC",
                                   "EXTRA"] = "BASIC"
    

    def __repr__(self):
        return super().__repr__()
    
@dataclass
class Solid187(Etype):

    name : str = "SOLID187"

    formulation : Literal["PURE DISPLACEMENT",
                          "MIXED UP"] = "MIXED UP"
    
    pml_absorbing : Literal["INCLUDED",
                            "NOT INCLUDED"] = "NOT INCLUDED"
    
    steady_state : Literal["ENABLED",
                          "DISABLED"] = "DISABLED"
    
    extra_surface_output : Literal["BASIC",
                                   "EXTRA"] = "BASIC"

    def __repr__(self):
        return super().__repr__()

@dataclass
class Solid285(Etype):

    name : str = "SOLID285"

    formulation : Literal["PURE DISPLACEMENT",
                          "MIXED UP WITH CONSTANT PRESSURE",
                          "MIXED UP WITH LINEAR PRESSURE"] = "PURE DISPLACEMENT"
    
    pml_absorbing : Literal["INCLUDED",
                            "NOT INCLUDED"] = "NOT INCLUDED"
    
    steady_state : Literal["ENABLED",
                          "DISABLED"] = "DISABLED"
    
    extra_surface_output : Literal["BASIC",
                                   "EXTRA"] = "BASIC"
    def __repr__(self):
        return super().__repr__()

@dataclass
class Plane182(Etype):

    name : str = "PLANE182"

    technology : Literal["FULL INTEGRATION",
                         "REDUCED INTEGRATION",
                         "ENHANCED STRAIN",
                         "SIMPLIFIED ENHANCED STRAIN"] = "FULL INTEGRATION"
    
    shape : Literal["QUAD 4",
                    "TRIA 3"] = "QUAD 4"

    behaviour : Literal["PLANE STRESS",
                       "AXISYMMETRIC",
                        "PLANE STRAIN",
                        "PLANE STRESS WITH THICKNESS",
                        "GENERALIZED PLANE STRAIN",
                        "AXISYMMETRIC WITH TORSION"] = "PLANE STRESS"
    
    formulation : Literal["PURE DISPLACEMENT",
                          "MIXED UP"] = "PURE DISPLACEMENT"
    
    pml_absorbing : Literal["INCLUDED",
                            "NOT INCLUDED"] = "NOT INCLUDED"
    
    extra_surface_output : Literal["BASIC",
                                   "EXTRA"] = "BASIC"

    def __repr__(self):
        return super().__repr__()

@dataclass
class Plane183(Etype):

    name : str = "PLANE183"

    shape : Literal["QUAD 8",
                    "TRIA 6"] = "QUAD 8"
    
    behaviour : Literal["PLANE STRESS",
                       "AXISYMMETRIC",
                        "PLANE STRAIN",
                        "PLANE STRESS WITH THICKNESS",
                        "GENERALIZED PLANE STRAIN",
                        "AXISYMMETRIC WITH TORSION"] = "PLANE STRESS"
    
    formulation : Literal["PURE DISPLACEMENT",
                          "MIXED UP"] = "PURE DISPLACEMENT"
    
    pml_absorbing : Literal["INCLUDED",
                            "NOT INCLUDED"] = "NOT INCLUDED"
    
    extra_surface_output : Literal["BASIC",
                                   "EXTRA"] = "BASIC"

    def __repr__(self):
        return super().__repr__()
    
##################################
### Surface Effect Elements
##################################

@dataclass
class Surf153(Etype):
    """2-D Structural Surface Effect"""
    name : str = "SURF153"
    pass

#####################
## SURF154 & INPUTS
#####################

SURF154_KEYOPT_MAPPING = {
    "pressure_cs": (2, ("ELEMENT CS", "LOCAL CS")),
    "midside_nodes": (4, ("INCLUDE", "EXCLUDE")),
    "normal_pressure_direction": (6, ("CALCULATED", "POSITIVE ONLY", "NEGATIVE ONLY")),
    "large_deflection_area": (7, ("USE NEW AREA", "ORIGINAL AREA")),
    "ocean_loading": (8, ("IGNORE", "APPLY")),
    "press_vector_orient": (11, ("PRJ AREA + INCLUDE TANGENT",
                                 "PRJ AREA + EXCLUDE TANGENT",
                                 "FULL AREA + INCLUDE TANGENT")),
    "effect_of_direction": (12, ("PRESSURE APPLIED", "PRESSURE SUPPRESSED"))
}

@dataclass
class Surf154RealConstants(Real):

    efs : float | str = ''  # Foundation stiffness
    surt : float | str = '' # Surface Tension
    admsua : float | str = ''   # Added mass/unit area
    tki : float | str = ''  # Thickness at node i (defaults to tki)
    tkj : float | str = ''  # Thickness at node j (defaults to tki)
    tkk : float | str = ''  # Thickness at node k (defaults to tki)
    tkl : float | str = ''  # Thickness at node l (defaults to tki)

    def get_real_constants(self) -> tuple:
        """Create a 8x6 tuple of real constants for CONTA174 element."""
        return ((self.efs, self.surt, self.admsua, self.tki, self.tkj, self.tkk),
                (self.tkl, "", "", "", "", ""))

@dataclass
class Surf154Material(Mat):
    dens : float = 0.0
    visc : float = 0.0

@dataclass
class Surf154(Etype):
    """3-D Structural Surface Effect"""
    name : str = "SURF154"

    pressure_cs : Literal["ELEMENT CS", "LOCAL CS"] = "ELEMENT CS" # Element cs follows large deflection effects
    midside_nodes : Literal["INCLUDE", "EXCLUDE"] = "INCLUDE"
    normal_pressure_direction : Literal["CALCULATED", "POSITIVE ONLY", "NEGATIVE ONLY"] | None = None
    large_deflection_area : Literal["USE NEW AREA", "ORIGINAL AREA"] = "USE NEW AREA"
    ocean_loading : Literal["IGNORE", "APPLY"] = "IGNORE"
    press_vector_orient : Literal["PRJ AREA + INCLUDE TANGENT",
                                  "PRJ AREA + EXCLUDE TANGENT",
                                  "FULL AREA + INCLUDE TANGENT"] | None = None
    effect_of_direction : Literal["PRESSURE APPLIED", "PRESSURE SUPPRESSED"] | None = None
    
    def get_keyopts(self) -> list:
        keyopts = []
        for key, value in SURF154_KEYOPT_MAPPING.items():
            value = getattr(self, key)
            if value is not None: 
                keyopt_number, options = SURF154_KEYOPT_MAPPING[key]
                keyopt_value = options.index(value)
                keyopts.append((keyopt_number, keyopt_value))
        return keyopts
                    
    def __repr__(self):
        return super().__repr__()

#####################
### MPC184 & INPUTS
#####################

MPC184_KEYOPT_MAPPING = {
    1 : ("behaviour", ("RigidLink", "RigidBeam", "Slider", "Revolute", "Universal", "Slot", "Point", "Translational", "Cylindrical", "Planar", "Weld", "Orient", "Spherical", "General", "Screw")),
    2 : ("method", ("Direct Elimination", "Lagrange Multiplier", "Penalty Based Method")),
    3 : ("configuration", ("x-axis revolute", "y-axis revolute", "displacement and rotational", "only displacement")),
    5 : ("stress_stiffness", ("Include", "Exclude"))}

@dataclass
class MPC184(Etype):
    
    name : str = "MPC184"

    behaviour : int = 0
    method : Literal["Direct Elimination", "Lagrange Multiplier", "Penalty Based Method"] = "Lagrange Multiplier"
    contribution : Literal["Include", "Exclude"] = "Include"

    def get_keyopts(self) -> tuple:
        return ((1, self.behaviour),
                (2, ["Direct Elimination", "Lagrange Multiplier", "Penalty Based Method"].index(self.method)),
                (3, ["Include", "Exclude"].index(self.contribution)))

    def __repr__(self):
        return super().__repr__()

##################################
### CONTACT and TARGET Elements
##################################
@dataclass
class ContaRealConstants(Real):
    """REAL CONSTANTS for CONTA174 element."""

    r1 : float | str = ''  # Target radius for cylinder, cone, or sphere
    r2 : float | str = ''  # Target radius at second node of cone
    fkn : float | str = '' # Normal penalty stiffness factor / weightenin factor for MPC_RBE3
    ftoln : float | str = '' # Penetration tolerance factor
    icont : float | str = '' # Initial contact closure
    pinb : float | str = '' # Pinball region
    pzer : float | str = '' # Pressure at zero penetration
    czer : float | str = '' # Initial contact clearance
    taumax : float | str = '' # Maximum friction stress
    cnof : float | str = ''  # Contact surface offset
    fkop : float | str = ''  # Contact openning stiffness
    fkt : float | str = '' # Tangential penalty stiffness factor
    cohe : float | str = '' # Contact cohesion
    fact : float | str = '' # Static/Dynamic ratio
    dc : float | str = '' # Exponential decay coefficient
    slto : float | str = '' # Allowable elastic slip
    tnop : float | str = '' # Maximum allowable tensile contact pressure
    tols : float | str = '' # Target edge extension factor
    ppcn : float | str = '' # Pressure penetration criteria
    fpat : float | str = '' # Fluid penetration acting time
    cor : float | str = '' # Coefficient of restitution
    strm : float | str = '' # Load step number for ramping penetration OR 
                        # Starting time for contact stiffness ramping
    fdmn : float | str = '' # Normal stabilization damping factor
    fdmt : float | str = '' # Tangential stabilization damping factor
    fdmd : float | str = '' # Destabilization squeal damping factor
    fdms : float | str = '' # Stabilization squeal damping factor
    bsrl : int | str = '' # Original contact pair real constant ID (after contact splitting)
    ksym : int | str = '' # Real constant ID of the associated companion pair for symmetric contact or self contact (after contact splitting)
    tfor : float | str = '' # Pair-based force tolerance
    tend : float | str = '' # End time for ramping contact stiffness

    def get_real_constants(self) -> tuple:
        """Create a 8x6 tuple of real constants for CONTA174 element."""
        return ((self.r1, self.r2, self.fkn, self.ftoln, self.icont, self.pinb),
                (self.pzer, self.czer, self.taumax, self.cnof, self.fkop, self.fkt),
                (self.cohe, '', '', '', '', ''),
                ('', '', self.fact, self.dc, self.slto, self.tnop),
                (self.tols, '', self.ppcn, self.fpat, self.cor, self.strm),
                (self.fdmn, self.fdmt, self.fdmd, self.fdms, '', ''),
                ('', '', '', '', '', ''),
                ('', '', self.bsrl, self.ksym, self.tfor, self.tend))

@dataclass
class Conta172(Etype):
    "2-D 3-Node Surface-to-Surface Contact"
    name : str = 'CONTA172'

    pass

#####################
## CONTA174 & INPUTS
#####################

CONTA174_KEYOPT_MAPPING = {
    "dof" : (1, ("UX UY UZ", 0), ("UX UY UZ TEMP", 1), ("UX UY UZ PRES", 8), ("UX UY UZ PRES TEMP", 9)),
    "contact_algorithm" : (2, ("Augmented Lagrange", 0), ("Penalty", 1), ("MPC", 2), ("Lagrange & Penalty", 3), ("Lagrange", 4)),
    "units_of_normal_contact_stiffness" : (3, ("FORCE/LENGTH^3", 0), ("FORCE/LENGTH", 1)),
    "contact_detection" : (4, ("At Gauss Point", 0), ("Node: Normal From Contact/RBE3", 1), ("Node: Normal To Target/CERIGID", 2), ("Node: Surface Projection/CP", 3), ("Node: Dual Shape Projection", 4)),
    "surface_based_constraint" : (4, ("Force Distribution", 1), ("Rigid Surface", 2), ("Coupling", 3)),
    "result_section" : (4, ("Moves with Contact Elements", 1), ("Moves with Underlying Elements", 2)),
    "auto_adjustment" : (5, ("No Auto Adjustment", 0), ("Close_Gap", 1), ("Reduce_Penetration", 2), ("Close Gap/Reduce Penetration", 3), ("Default ICONT", 4)),
    "auto_contact_stiffness_change" : (6, ("Standard", 0), ("Aggressive", 1), ("Very Aggressive", 2), ("Exponential", 3)),
    "contact_time_load_prediction" : (7, ("No Predictions", 0), ("Automatic Bisection", 1), ("Reasonable T/L increment", 2), ("Minimum T/L increment", 3), ("Impact Constraint", 4)),
    "asymmetric_contact_active" : (8, ("Yes", 2)),
    "initial_penetration_gap" : (9, ("Include",0), ("Exclude", 1), ("Include - Ramp",2), ("Offset Only", 3), ("Offset - Ramp", 4), ("Always Offset Only", 5), ("Always Offset - Ramp", 6)),
    "contact_stiffness_update" : (10, ("Each Iteration", 0), ("Each Load Step", 1), ("Each Iteration & Slip Limit", 2)),
    "shell_thickness_effect" : (11, ("Exclude", 0), ("Include", 1)),
    "behaviour" : (12, ("Standard", 0), ("Rough", 1), ("No Seperation", 2), ("Bonded", 3), ("No Seperation (always)", 4), ("Bonded (always)", 5), ("Bonded (initial)", 6)),
    "fluid_penetration_load" : (14, ("Iteration Based", 0), ("Substep Based", 1), ("Iteration & From Starting Points", 2), ("Substep & From Strating Points", 3)),
    "effect_of_stabilization_damping" : (15, ("Active in 1st Load Step", 0), ("Deactive Auto Damping", 1), ("Active For All Load Steps", 2), ("Active All the Time", 3)),
    "input_of_squeal_dampings" : (16, ("Damping Factor",0), ("MU-Slip velocity", 1), ("Damping Coefficient", 2)),
    "sliding_behaviour" : (18, ("Finite Sliding", 0), ("Small Sliding", 1), ("Adaptive Small Sliding", 2))
}

@dataclass
class Conta174Material(Mat):
    mu_1 : float = 0.2
    mu_2 : float = 0.0

@dataclass
class Conta174(Etype):
    "3-D 8-Node Surface-to-Surface Contact"
    name : str = "CONTA174"

    dof : Literal["UX UY UZ",
                  "UX UY UZ TEMP",
                  "UX UY UZ PRES",
                  "UX UY UZ PRES TEMP"] = "UX UY UZ" # keyopt1
    
    contact_algorithm : Literal["Augmented Lagrange", 
                                "Penalty", 
                                "MPC",
                                "Lagrange & Penalty", 
                                "Lagrange"] = "Augmented Lagrange" # keyopt2
    
    units_of_normal_contact_stiffness : Literal["FORCE/LENGTH^3", "FORCE/LENGTH"] = "FORCE/LENGTH^3" # keyopt3

    contact_detection : Literal["At Gauss Point", 
                                "Node: Normal From Contact/RBE3",
                                "Node: Normal To Target/CERIGID",
                                "Node: Surface Projection/CP",
                                "Node: Dual Shape Projection"] = "At Gauss Point" # keyopt4
    
    # Alternate keyopt4 usage with Surface-Based Constraints
    # To define surface based constraints :: For keyopt2 = MPC(2) or Lagrange_Multiplier(3)
    surface_based_constraint : Literal["Force Distribution",
                                       "Rigid Surface",
                                       "Coupling"] | None = None # keyopt4

    # Alternate keyopt4 usage with result sections
    # result section :: How motion of the result section is accounted for in large deflection analysis
    result_section : Literal["Moves with Contact Elements",
                             "Moves with Underlying Elements"] | None = None # keyopt4

    auto_adjustment : Literal["No Auto Adjustment",
                              "Close_Gap",
                              "Reduce_Penetration",
                              "Close Gap/Reduce Penetration",
                              "Default ICONT"] | None = None # keyopt5
    
    auto_contact_stiffness_change : Literal["Standard", 
                                            "Aggressive",
                                            "Very Aggressive",
                                            "Exponential"] | None = None # keyopt6
    
    contact_time_load_prediction : Literal["No Predictions",
                                           "Automatic Bisection",
                                           "Reasonable T/L increment",
                                           "Minimum T/L increment",
                                           "Impact Constraint"] | None = None # keyopt7
    
    asymmetric_contact_active : Literal["Yes"] | None = None # keyopt8
    
    initial_penetration_gap : Literal["Include",
                                      "Exclude",
                                      "Include - Ramp",
                                      "Offset Only",
                                      "Offset - Ramp"
                                      "Always Offset Only",
                                      "Always Offset - Ramp"] | None = None # keyopt9
    
    contact_stiffness_update : Literal["Each Iteration",
                                       "Each Load Step",
                                       "Each Iteration & Slip Limit"] | None = None # keyopt10
    
    shell_thickness_effect : Literal["Exclude", "Include"] | None = None # keyopt11

    behaviour : Literal["Standard", 
                        "Rough", 
                        "No Seperation", 
                        "Bonded", 
                        "No Seperation (always)",
                        "Bonded (always)",
                        "Bonded (initial)"] | None = None # keyopt12
    
    # dof_control_for_thermal_shells : keyopt13
    
    fluid_penetration_load : Literal["Iteration Based",
                                     "Substep Based",
                                     "Iteration & From Starting Points",
                                     "Substep & From Strating Points"] | None = None # keyopt14

    effect_of_stabilization_damping : Literal["Active in 1st Load Step",
                                              "Deactive Auto Damping",
                                              "Active For All Load Steps",
                                              "Active All the Time"] | None = None # keyopt15
    
    input_of_squeal_dampings : Literal["Damping Factor",
                                       "MU-Slip velocity",
                                       "Damping Coefficient"] | None = None # keyopt16
    
    sliding_behaviour : Literal["Finite Sliding",
                                "Small Sliding",
                                "Adaptive Small Sliding"] = "Finite Sliding" # keyopt18
    
    def get_keyopts(self) -> list | None:
        keyopts = []
        for key in CONTA174_KEYOPT_MAPPING.keys():
            if getattr(self, key) == None:
                continue # do not add keyopt into keyopts list if its value is 'None'
            else:
                _keyopt_id = CONTA174_KEYOPT_MAPPING[key][0]
                _keyopt_val = next(val[1] for val in CONTA174_KEYOPT_MAPPING[key][1:] if key == val[0])
                keyopts.append((_keyopt_id, _keyopt_val))
        return keyopts or None

#####################
## CONTA175 & INPUTS
#####################

CONTA175_KEYOPT_MAPPING = {
    "dof" : (1, ("UX UY UZ", 0), ("UX UY UZ TEMP", 1), ("UX UY UZ PRES", 8), ("UX UY UZ PRES TEMP", 9)),
    "contact_algorithm" : (2, ("Augmented Lagrange", 0), ("Penalty", 1), ("MPC", 2), ("Lagrange & Penalty", 3), ("Lagrange", 4)),
    "contact_model" : (3, ("Contact Force Based", 0), ("Contact Traction Based", 1)),
    "contact_normal_direction" : (4, ("Normal To Target Surface", 0), ("Normal From Contact Nodes", 1), ("Normal From Contact Nodes (Shell/Beam)", 2), ("Normal To Target Surface (Shell/Beam)", 3)),
    "surface_based_constraint" : (4, ("Rigid Surface", 0), ("Force Distribution", 1), ("Coupling", 3)),
    "auto_adjustment" : (5, ("No Auto Adjustment", 0), ("Close_Gap", 1), ("Reduce_Penetration", 2), ("Close Gap/Reduce Penetration", 3), ("Default ICONT", 4)),
    "auto_contact_stiffness_change" : (6, ("Standard", 0), ("Aggressive", 1), ("Very Aggressive", 2), ("Exponential", 3)),
    "contact_time_load_prediction" : (7, ("No Predictions", 0), ("Automatic Bisection", 1), ("Reasonable T/L increment", 2), ("Minimum T/L increment", 3), ("Impact Constraint", 4)),
    "asymmetric_contact_active" : (8, ("Yes", 2)),
    "initial_penetration_gap" : (9, ("Include",0), ("Exclude", 1), ("Include - Ramp",2), ("Offset Only", 3), ("Offset - Ramp", 4), ("Always Offset Only", 5), ("Always Offset - Ramp", 6)),
    "contact_stiffness_update" : (10, ("Each Iteration", 0), ("Each Load Step", 1), ("Each Iteration & Slip Limit", 2)),
    "shell_thickness_effect" : (11, ("Exclude", 0), ("Include", 1)),
    "behaviour" : (12, ("Standard", 0), ("Rough", 1), ("No Seperation", 2), ("Bonded", 3), ("No Seperation (always)", 4), ("Bonded (always)", 5), ("Bonded (initial)", 6)),
    "fluid_penetration_load" : (14, ("Iteration Based", 0), ("Substep Based", 1), ("Iteration & From Starting Points", 2), ("Substep & From Strating Points", 3)),
    "effect_of_stabilization_damping" : (15, ("Active in 1st Load Step", 0), ("Deactive Auto Damping", 1), ("Active For All Load Steps", 2), ("Active All the Time", 3)),
    "input_of_squeal_dampings" : (16, ("Damping Factor",0), ("MU-Slip velocity", 1), ("Damping Coefficient", 2)),
    "sliding_behaviour" : (18, ("Finite Sliding", 0), ("Small Sliding", 1), ("Adaptive Small Sliding", 2))
}

@dataclass
class Conta175Material(Mat):
    mu_1 : float = 0.2
    mu_2 : float = 0.0

@dataclass
class Conta175(Etype):
    """
    2-D/3-D Node-to-Surface Contact
    EXAMPLE :: 
    keyo,cid,12,5              ! Bonded Contact (Bonded Always)
    keyo,cid,4,1               ! Deformable RBE3 style load (Force Distribution)
    keyo,cid,2,2               ! MPC style contact (MPC)
    """

    name : str = "CONTA175"

    dof : Literal["UX UY UZ", # Default
                  "UX UY UZ TEMP",
                  "UX UY UZ PRES",
                  "UX UY UZ PRES TEMP"] | None = None # keyopt1
    
    contact_algorithm : Literal["Augmented Lagrange", # Default
                                "Penalty", 
                                "MPC",
                                "Lagrange & Penalty", 
                                "Lagrange"] | None = None # keyopt2
    
    contact_model : Literal["Contact Force Based", # Default
                            "Contact Traction Based"] | None = None # keyopt3

    contact_normal_direction : Literal["Normal To Target Surface",  # Default
                                       "Normal From Contact Nodes",
                                       "Normal From Contact Nodes (Shell/Beam)",
                                       "Normal To Target Surface (Shell/Beam)"] | None = None # keyopt4
    
    # Alternate keyopt4 usage with Surface-Based Constraints
    surface_based_constraint : Literal["Rigid Surface",
                                       "Force Distribution",
                                       "Coupling"] | None = None # keyopt4

    auto_adjustment : Literal["No Auto Adjustment",
                              "Close_Gap",
                              "Reduce_Penetration",
                              "Close Gap/Reduce Penetration",
                              "Default ICONT"] | None = None # keyopt5
    
    auto_contact_stiffness_change : Literal["Standard", 
                                            "Aggressive",
                                            "Very Aggressive",
                                            "Exponential"] | None = None # keyopt6
    
    contact_time_load_prediction : Literal["No Predictions",
                                           "Automatic Bisection",
                                           "Reasonable T/L increment",
                                           "Minimum T/L increment",
                                           "Impact Constraint"] | None = None # keyopt7
    
    asymmetric_contact_active : Literal["Yes"] | None = None # keyopt8
    
    initial_penetration_gap : Literal["Include",
                                      "Exclude",
                                      "Include - Ramp",
                                      "Offset Only",
                                      "Offset - Ramp"
                                      "Always Offset Only",
                                      "Always Offset - Ramp"] | None = None # keyopt9
    
    contact_stiffness_update : Literal["Each Iteration",
                                       "Each Load Step",
                                       "Each Iteration & Slip Limit"] | None = None # keyopt10
    
    shell_thickness_effect : Literal["Exclude", "Include"] | None = None # keyopt11

    behaviour : Literal["Standard", 
                        "Rough", 
                        "No Seperation", 
                        "Bonded", 
                        "No Seperation (always)",
                        "Bonded (always)",
                        "Bonded (initial)"] | None = None # keyopt12
    
    # dof_control_for_thermal_shells : keyopt13
    
    fluid_penetration_load : Literal["Iteration Based",
                                     "Substep Based",
                                     "Iteration & From Starting Points",
                                     "Substep & From Strating Points"] | None = None # keyopt14

    effect_of_stabilization_damping : Literal["Active in 1st Load Step",
                                              "Deactive Auto Damping",
                                              "Active For All Load Steps",
                                              "Active All the Time"] | None = None # keyopt15
    
    input_of_squeal_dampings : Literal["Damping Factor",
                                       "MU-Slip velocity",
                                       "Damping Coefficient"] | None = None # keyopt16
    
    sliding_behaviour : Literal["Finite Sliding", # Default
                                "Small Sliding",
                                "Adaptive Small Sliding"] | None = None # keyopt18
        
    def get_keyopts(self) -> list | None:
        keyopts = []
        for key in CONTA175_KEYOPT_MAPPING.keys():
            _value = getattr(self, key)
            if _value == None:
                continue # do not add keyopt into keyopts list if its value is 'None'
            else:
                _keyopt_id = CONTA175_KEYOPT_MAPPING[key][0]
                _keyopt_val = next(val[1] for val in CONTA175_KEYOPT_MAPPING[key][1:] if _value == val[0])
                keyopts.append((_keyopt_id, _keyopt_val))
        return keyopts or None

#####################
## CONTA177 & INPUTS
#####################

@dataclass
class Conta177(Etype):
    "3-D Line-to-Surface Contact"
    name : str = 'CONTA177'


@dataclass
class Targe169(Etype):

    name : str = 'TARGE169'
    pass

TARGE170_KEYOPT_MAPPING = {
    "element_order" : (1, ("Low Order", 0), ("High Order", 1)),
    "bc_for_rigid_target_nodes" : (2, ("Auto Constraint", 0), ("User Definition", 1)),
    "dof" : None,
    "type_of_constraint" : (5, ("Auto", 0), ("Proj Constraint U Dof Only", 1), ("Proj Constraint U+R Dof", 2), ("Disp Constraint (Norm Dir)", 3), ("Disp Constraint (All Dir)", 4), ("Disp Constraint (Anywhere)", 5)),
    "constrained_surface_symmetry" : (6, ("Pilot XY Plane", "100"), ("Pilot XZ Plane", "10"), ("Pilot YZ Plane", "1"), ("Pilot XY + XZ", "110"), ("Pilot XY + YZ", "101"), ("Pilot XZ + YZ", "11")),
    "weighting_factor" : (7, ("Internal", 0), ("1", 1), ("User Defined", 2)),
    "beam_to_beam_contact_type" : (9, ("External", 0), ("Internal", 1)),
    "stress_stiffening_effect" : (10, ("Do Not Include", 0), ("Include", 1)),
    "relaxion_method" : (11, ("Do Not Use", 0), ("Use", 1)),
    "thermal_expansion_effect" : (12, ("Do Not Include", 0), ("Include", 1))
}

@dataclass
class Targe170(Etype):
    """
    EXAMPLE
    keyo,tid,2,1               ! Don't fix the pilot node (User Definition) ????
    keyo,tid,4,0               ! Activate all DOF's due to large deformation (0 actives all!!!)
    """

    name : str = "TARGE170"

    element_order : Literal["Low Order", "High Order"] = "Low Order" #K1

    bc_for_rigid_target_nodes : Literal["Auto Constraint", "User Definition"] = "Auto Constraint" #K2 

    # Thermal contact behaviour #K3

    dof :  tuple = ("UX", "UY", "UZ", "ROTX", "ROTY", "ROTZ") #K4 

    # This keyopt(5) is not used for surface-based constraints.
    type_of_constraint : Literal["Auto", 
                                 "Proj Constraint U Dof Only",
                                 "Proj Constraint U+R Dof",
                                 "Disp Constraint (Norm Dir)",
                                 "Disp Constraint (All Dir)",
                                 "Disp Constraint (Anywhere)"] = "Auto" #K5

    constrained_surface_symmetry : Literal["Pilot XY Plane",
                                           "Pilot XZ Plane",
                                           "Pilot YZ Plane",
                                           "Pilot XY + XZ",
                                           "Pilot XY + YZ",
                                           "Pilot XZ + YZ"] | None = None #K6

    weighting_factor : Literal["Internal", "1", "User Defined"] | None = None #K7

    beam_to_beam_contact_type : Literal["External", "Internal"] | None = None #K9

    stress_stiffening_effect : Literal["Do Not Include", "Include"] | None = None #K10

    relaxion_method : Literal["Do Not Use", "Use"] | None = None  #K11

    thermal_expansion_effect : Literal["Do Not Include", "Include"] | None = None #K12


    def __post_init__(self):
        # post dof setting 
        ux = uy = uz = rotx = roty = rotz = 0
        if "UX" in self.dof:
            ux = 1
        if "UY" in self.dof:
            uy = 1
        if "UZ" in self.dof:
            uz = 1
        if "ROTX" in self.dof:
            rotx = 1
        if "ROTY" in self.dof:
            roty = 1
        if "ROTZ" in self.dof:
            rotz = 1
        TARGE170_KEYOPT_MAPPING["dof"]= (4, (self.dof, f'{ux}{uy}{uz}{rotx}{roty}{rotz}'))
        #TARGE170_KEYOPT_MAPPING["dof"]= (4, (self.dof, f'0'))

    def get_keyopts(self) -> list | None:
        keyopts = []
        for key in TARGE170_KEYOPT_MAPPING.keys():
            _value = getattr(self, key)
            if _value == None:
                continue # do not add keyopt into keyopts list if its value is 'None'
            else:
                _keyopt_id = TARGE170_KEYOPT_MAPPING[key][0]
                _keyopt_val = next(val[1] for val in TARGE170_KEYOPT_MAPPING[key][1:] if _value == val[0])
                keyopts.append((_keyopt_id, _keyopt_val))
        return keyopts or None
        

########################################################
### SECTION Definitons for the elements with two points
########################################################

@dataclass
class SectionData:
    """
    SECDATA, VAL1, VAL2, VAL3, VAL4, VAL5, VAL6, VAL7, VAL8, VAL9, VAL10, VAL11, VAL12
    Describes the geometry of a section.
    https://ansyshelp.ansys.com/public/account/secured?returnurl=///Views/Secured/corp/v242/en/ans_cmd/Hlp_C_SECDATA.html
    """
    pass

"""####################################
### MPC184 ELEMENT CONFIGURATION
####################################

MPCtypes = {
            "RigidLink" : 0,
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
            "Screw" : 17
            }

behaviour = Literal["RigidLink",
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
                    "Screw"]
                    
                
method = Literal["Direct Elimination",
                "Lagrange Multiplier",
                "Penalty Based Method"]
                

config = Literal["x-axis revolute",
                 "y-axis revolute",
                 "displacement and rotational",
                 "only displacement"]

stress_stiffness = Literal["Include",
                           "Exclude"]"""

"""@dataclass
class MPC184:
    
    keyopt1 : behaviour # Element behavior
    keyopt2 : method # Element constraint imposition method
    keyopt4 : None | config # Element configuration
    keyopt5 : None | stress_stiffness  # Geometric stress-stiffness contributions for rigid beams"""


VALID_ETYPES = {
    154: Surf154,
    185: Solid185,
    186: Solid186,
    187: Solid187,
    285: Solid285,
    182: Plane182,
    183: Plane183
}

name : Literal["VERTEX", # PyVista Cell Types
                   "LINE",
                   "TRIA",
                   "QUAD",
                   "TETRA",
                   "HEX",
                   "WEDGE",
                   "PYRAMID",
                   "SURFACE"]
_etype_map = {
    "Surf154" : "SURFACE",
    "Solid185" : "HEX",
    "Solid186" : "HEX",
    "Solid187" : "HEX",
    "Solid285" : ""
}



"""
# element type to VTK conversion function call map
# 0: skip
# 1: Point
# 2: Line (linear or quadratic)
# 3: Shell
# 4: 3D Solid (Hexahedral, wedge, pyramid, tetrahedral)
# 5: Tetrahedral
# 6: Line (always linear)

_etype_map = [
    0,
    2,  # LINK1
    3,  # PLANE2
    2,  # BEAM3
    2,  # BEAM4
    4,  # SOLID5
    0,  # UNUSED6
    1,  # COMBIN7
    2,  # LINK8
    0,  # INFIN9
    2,  # LINK10
    2,  # LINK11
    2,  # CONTAC12
    3,  # PLANE13
    2,  # COMBIN14
    0,  # FLUID15
    2,  # PIPE16
    2,  # PIPE17
    2,  # PIPE18
    0,  # SURF19
    2,  # PIPE20
    1,  # MASS21
    3,  # SURF22
    2,  # BEAM23
    2,  # BEAM24
    3,  # PLANE25
    0,  # UNUSED26
    0,  # MATRIX27
    3,  # SHELL28
    3,  # FLUID29
    4,  # FLUID30
    2,  # LINK31
    2,  # LINK32
    2,  # LINK33
    2,  # LINK34
    3,  # PLANE35
    0,  # SOURC36
    2,  # COMBIN37
    2,  # FLUID38
    2,  # COMBIN39
    2,  # COMBIN40
    3,  # SHELL41
    3,  # PLANE42
    3,  # SHELL43
    2,  # BEAM44
    4,  # SOLID45
    4,  # SOLID46
    0,  # INFIN47
    0,  # UNUSED48
    0,  # UNUSED49
    0,  # MATRIX50
    3,  # SHELL51
    0,  # CONTAC52
    3,  # PLANE53
    3,  # BEAM54
    3,  # PLANE55
    0,  # UNUSED56
    3,  # SHELL57
    0,  # UNUSED58
    2,  # PIPE59
    2,  # PIPE60
    2,  # SHELL61
    4,  # SOLID62
    3,  # SHELL63
    4,  # SOLID64
    4,  # SOLID65
    2,  # FLUID66
    3,  # PLANE67
    2,  # LINK68
    4,  # SOLID69
    4,  # SOLID70
    1,  # MASS71
    0,  # UNUSED72
    0,  # UNUSED73
    0,  # UNUSED74
    3,  # PLANE75
    0,  # UNUSED76
    3,  # PLANE77
    3,  # PLANE78
    3,  # FLUID79
    4,  # FLUID80
    3,  # FLUID81
    3,  # PLANE82
    3,  # PLANE83
    0,  # UNUSED84
    0,  # UNUSED85
    0,  # UNUSED86
    5,  # SOLID87
    3,  # VISCO88
    4,  # VISCO89
    4,  # SOLID90
    3,  # SHELL91
    5,  # SOLID92
    3,  # SHELL93
    0,  # CIRCU94
    4,  # SOLID95
    4,  # SOLID96
    4,  # SOLID97
    5,  # SOLID98
    3,  # SHELL99
    4,  # USER100
    4,  # USER101
    4,  # USER102
    4,  # USER103
    4,  # USER104
    4,  # USER105
    3,  # VISCO106
    4,  # VISCO107
    4,  # VISCO108
    0,  # TRANS109
    0,  # INFIN110
    0,  # INFIN111
    0,  # ROT112
    0,  # UNUSED113
    0,  # UNUSED114
    3,  # INTER115
    2,  # FLUID116
    4,  # EDGE117
    3,  # HF118
    5,  # HF119
    4,  # HF120
    3,  # PLANE121
    4,  # SOLID122
    5,  # SOLID123
    0,  # CIRCU124
    0,  # CIRCU125
    2,  # TRANS126
    0,  # UNUSED127
    0,  # UNUSED128
    2,  # FLUID129
    3,  # FLUID130
    3,  # SHELL131
    3,  # SHELL132
    0,  # UNUSED133
    0,  # UNUSED134
    0,  # TRANS135
    3,  # FLUID136
    0,  # FLUID137
    0,  # FLUID138
    0,  # FLUID139
    5,  # ROM140
    0,  # UNUSED141
    0,  # UNUSED142
    3,  # SHELL143
    0,  # ROM144
    0,  # UNUSED145
    0,  # UNUSED146
    0,  # UNUSED147
    0,  # UNUSED148
    0,  # UNUSED149
    0,  # UNUSED150
    2,  # SURF151
    3,  # SURF152
    2,  # SURF153
    3,  # SURF154
    3,  # SURF155
    2,  # SURF156
    3,  # SHELL157
    0,  # UNUSED158
    0,  # SURF159
    0,  # UNUSED160
    2,  # BEAM161
    0,  # UNUSED162
    3,  # SHELL163
    4,  # SOLID164
    0,  # UNUSED165
    0,  # UNUSED166
    0,  # UNUSED167
    5,  # SOLID168
    0,  # TARGE169
    3,  # TARGE170
    2,  # CONTA171
    2,  # CONTA172
    3,  # CONTA173
    3,  # CONTA174
    1,  # CONTA175
    2,  # CONTA176
    2,  # CONTA177
    2,  # CONTA178
    0,  # PRETS179
    2,  # LINK180
    3,  # SHELL181
    3,  # PLANE182
    3,  # PLANE183
    0,  # RBAR184
    4,  # SOLID185
    4,  # SOLID186
    5,  # SOLID187
    6,  # BEAM188
    2,  # BEAM189
    4,  # SOLSH190
    0,  # UNUSED191
    4,  # INTER192
    0,  # INTER193
    0,  # INTER194
    0,  # INTER195
    0,  # LAYER196
    0,  # LAYER197
    0,  # UNUSED198
    0,  # UNUSED199
    0,  # UNUSED200
    0,  # FOLLW201
    0,  # INTER202
    0,  # INTER203
    0,  # INTER204
    0,  # INTER205
    0,  # INTER206
    0,  # UNUSED207
    2,  # SHELL208
    2,  # SHELL209
    0,  # UNUSED210
    0,  # UNUSED211
    3,  # CPT212
    3,  # CPT213
    6,  # COMBI214
    4,  # CPT215
    6,  # CPT216
    6,  # CPT217
    3,  # FLUID218
    3,  # FLUID219
    4,  # FLUID220
    5,  # FLUID221
    3,  # PLANE222
    3,  # PLANE223
    0,  # UNUSED224
    0,  # SOLID225
    4,  # SOLID226
    5,  # SOLID227
    0,  # UNUSED228
    0,  # UNUSED229
    3,  # PLANE230
    4,  # SOLID231
    5,  # SOLID232
    3,  # PLANE233
    0,  # UNUSED234
    0,  # UNUSED235
    4,  # SOLID236
    5,  # SOLID237
    3,  # PLANE238
    4,  # SOLID239
    5,  # SOLID240
    0,  # HSFLD241
    0,  # HSFLD242
    0,  # UNUSED243
    0,  # UNUSED244
    0,  # UNUSED245
    0,  # UNUSED246
    0,  # UNUSED247
    0,  # UNUSED248
    0,  # UNUSED249
    2,  # COMBI250
    2,  # SURF251
    3,  # SURF252
    0,  # UNUSED253
    0,  # UNUSED254
    0,  # UNUSED255
    0,  # UNUSED256
    0,  # INFIN257
    0,  # INFIN258
    0,  # INFIN259
    0,  # UNUSED260
    4,  # UNUSED261
    0,  # GLINK262
    0,  # REINF263
    0,  # REINF264
    0,  # REINF265
    0,  # UNUSED266
    0,  # UNUSED267
    0,  # UNUSED268
    0,  # UNUSED269
    0,  # UNUSED270
    0,  # UNUSED271
    0,  # SOLID272
    0,  # SOLID273
    0,  # UNUSED274
    0,  # UNUSED275
    0,  # UNUSED276
    0,  # UNUSED277
    4,  # SOLID278
    4,  # SOLID279
    2,  # CABLE280
    3,  # SHELL281
    3,  # SHELL282
    3,  # SHELL283
    0,  # UNUSED284
    5,  # SOLID285
    0,  # UNUSED286
    0,  # UNUSED287
    2,  # PIPE288
    2,  # PIPE289
    2,  # ELBOW290
    5,  # SOLID291
    3,  # PLANE292
    3,  # PLANE293
    0,  # UNUSED294
    0,  # UNUSED295
    0,  # UNUSED296
    0,  # UNUSED297
    0,  # UNUSED298
    0,  # UNUSED299
    0,  # USER300
]
"""