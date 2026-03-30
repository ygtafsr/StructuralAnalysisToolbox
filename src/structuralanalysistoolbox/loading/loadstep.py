from __future__ import annotations
from dataclasses import dataclass, field
from typing import Literal, overload
from pathlib import Path

from structuralanalysistoolbox.mapdl import element, attributes
from structuralanalysistoolbox.exceptions import ParameterError


from typing import TYPE_CHECKING, overload

if TYPE_CHECKING:
    from structuralanalysistoolbox.model import Nset, Elset, Surface, CoordinateSystem, Vector, Analysis, Model

"""
## NONLINEAR OPTIONS:
NEQIT : Specifies the max num. of eq. iter. per substep.
CNVTOL : Specifies convergence tolerances.
NCNV : Provides options for terminating analysis.
CUTCONTROL,Lab,VALUE,Option
Controls time-step cutback during a nonlinear solution.

## OUTPUT CONTROLS:
OUTRES : Controls what the program writes to the OUTRES database and results file
         and how often it is written.
OUTPR : Controls what is printed (written to the solution OUTPR output file, (Jobname.OUT) and how
        often it is written.
OUTGEOM : Controls the type of geometry data the program writes to the results file.
ERESX : enables you to review element integration point values in the postprocessor.
"""

@dataclass
class Force:
    set : Nset | Surface | None = None
    value : float | None = None
    direction : Literal["X", "Y", "Z"] | None = None
    csys : CoordinateSystem | None = None
    applied_by : Literal["NODE SET", "ELEMENT FACES", "SURFACE ELEMENTS"] = "NODE SET"
    operation : Literal["ADD", "NEW", "DELETE"] = "NEW"
    fixed : bool = False

@dataclass
class Moment:
    set : Surface 
    value : float | None = None
    direction : Literal["RX", "RY", "RZ"] | None = None
    csys : CoordinateSystem | None = None
    behavior : Literal["Rigid", "Deformable", "Coupled", "Beam"] = "Rigid"
    operation : Literal["ADD", "NEW", "DELETE"] = "NEW"
    fixed : bool = False

@dataclass
class Pressure:
    set : Surface | None = None
    value : float | None = None
    direction : Literal["X", "Y", "Z", "NORMAL TO"] = "NORMAL TO"
    csys : CoordinateSystem | None = None
    applied_by : Literal["NODE SET", "ELEMENT FACES", "SURFACE ELEMENTS"] = "SURFACE ELEMENTS"
    operation : Literal["ADD", "NEW", "DELETE"] = "NEW"
    loaded_area : Literal["INITIAL", "DEFORMED"] = "DEFORMED"
    
@dataclass
class Displacement:
    set : Nset 
    value : float | None = None
    direction : Literal["ALL", "X", "Y", "Z", "RX", "RY", "RZ"] | None = None
    csys : CoordinateSystem | None = None
    operation : Literal["ADD", "NEW", "DELETE"] = "NEW"
    ramping : bool = False
       

class LoadStep:
    
    def __init__(self,
                 name : str,
                 model : Model,
                 step_number : int = 0,
                 end_time : int = 0, 
                 analysis_type : Literal["STATIC", "BUCKLE", "MODAL", "HARMIC" ,"TRANS", "SUBSTR"] = "STATIC"):
        
        self.name = name
        self.midx = 0
        self.model = model
        self.step_number = step_number
        self.analysis_type = analysis_type

        self.target_dir : Path | None = None # For RESTART
        self.filename : str | None = None # For RESTART
        self.status : str | None = None # For RESTART
        self.step : int | None = None # For RESTART
        self.substep : int | None = None # For RESTART

        #self.perturbation : Literal["YES", "NO"] = "NO" 
        self.nonlinear_geo : Literal["ON", "OFF"] = "ON" # NLGEOM can not be changed after the first SOLVE.
        self.loading_style : Literal["RAMPED", "STEPPED"] = "RAMPED"
        self.reference_temp : float = 0

        self.step_controls = StepControls(step_end_time=end_time)
        self.solver_controls = SolverControls()
        self.nonlinear_controls = NonlinearControls()
        self.restart_controls = RestartControls() # Define restart points for this LS
  
        self.resiudal_monitoring = IterationMonitoring(function="NRRE", characteristics="ON", maxfile=4)
        self.violation_monitoring = IterationMonitoring(function="EFLG", characteristics="ON", maxfile=4)
        self.contact_monitoring = IterationMonitoring(function="CONT", characteristics="ITER", maxfile=4)

        self.forces : list[Force] = []
        self.moments : list[Moment] = []
        self.pressures : list[Pressure] = []
        self.displacements : list[Displacement] = []
        self.outputs : list = []


    @overload
    def force(self, nset : str, 
                    value : float | None = None,
                    direction : Literal["X", "Y", "Z"] | None = None,
                    coordinate_system : CoordinateSystem | None = None,
                    applied_by : Literal["NODE SET"] = "NODE SET",
                    operation : Literal["ADD", "NEW", "DELETE"] = "NEW", 
                    fixed : bool = False, **kwargs): ...
    @overload
    def force(self, nset : str, 
                    value : float | None = None,
                    direction : Literal["X", "Y", "Z"] | None = None,
                    coordinate_system : CoordinateSystem | None = None,
                    applied_by : Literal["ELEMENT FACES", "SURFACE ELEMENTS"] = "ELEMENT FACES",
                    operation : Literal["ADD", "NEW", "DELETE"] = "NEW", 
                    **kwargs): ...
    def force(self, nset : str, 
                    value : float | None = None,
                    direction : Literal["X", "Y", "Z"] | None = None,
                    coordinate_system : CoordinateSystem | None = None,
                    applied_by : Literal["NODE SET", "ELEMENT FACES", "SURFACE ELEMENTS"] = "NODE SET",
                    operation : Literal["ADD", "NEW", "DELETE"] = "NEW", 
                    fixed : bool = False, **kwargs):
        """
            Force applied on the constant direction and onto initial surface area for 'surface forces'.
            NEW :: Replace previous Force, apply the new one.
            ADD :: Apply the new load in addition to the previous one.
            DELETE :: Deletes loads on nset
            FIXED :: Displacement on the specified dof is locked. (Only valid for the forces is applied by "NODE SET")
        """

        # Check whether nset is exist.
        _nset = self.model.get_nset(nset)
        if not _nset:
            raise ParameterError(f"Nset {nset} not found!")
        
        if applied_by == "NODE SET":

            # If there is no csys given, use global cartesian by default. 
            if coordinate_system:
                csys = coordinate_system
            else:
                csys = self.model.get_coordinate_system("Cartesian")

            force = Force(_nset, value, direction, csys, applied_by, operation, fixed)
            self.forces.append(force)
            
        elif applied_by == "ELEMENT FACES":
            # SFE (through "Element Faces")

            if coordinate_system:
                csys = coordinate_system
            else:
                csys = self.model.add_coordinate_system(show_in_model_tree=False)

            # If surface is already exist, use is. If not, create new one.
            surf = self.model._get_surface(nset, applied_by)
            if not surf:
                surf = self.model._add_surface(nset=nset, 
                                               surface_type="ELEMENT FACES",
                                               ignore_at_execution=True)
            
            force = Force(_nset, value, direction, csys, applied_by, operation, fixed)
            self.forces.append(force)
            
        elif applied_by == "SURFACE ELEMENTS":
            # Apply through Pressure on Surface Elements
            
            # Coordinate System Implementation::
            # 1) Element CS (Default):
            #       Element CS Rotates with the element
            #       Pressure direction rotates with the structure.
            #       Large deformation -> Follower pressure
            # 2) User Defined CS:
            #       The load does not rotate with the element.
            #       Pressure direction is fixed.
            #       Define local CS + assign it to the element using ESYS

            # If surface is already exist, use is. If not, create new one.
            surf = self.model._get_surface(nset, "SURFACE ELEMENTS")
            if not surf:
                # Element Type Selection::
                if kwargs and kwargs["element_type"]:
                    # Possible Types:
                    if kwargs["element_type"] == "Surf154":
                        surface_element = element.Surf154()
                    elif kwargs["element_type"] == "Surf153":
                        surface_element = element.Surf153()
                    elif kwargs["element_type"] == "Surf156":
                        surface_element = element.Surf156()
                    elif kwargs["element_type"] == "Surf159":
                        surface_element = element.Surf159()
                    elif kwargs["element_type"] == "Follw201":
                        surface_element = element.Follw201()
                    surface_element = element.Surf154() # Default
                else:
                    surface_element = element.Surf154() # Default

                surface_element.large_deflection_area = "ORIGINAL AREA"
                surface_element.midside_nodes = "INCLUDE"

                surface_element.pressure_cs = "LOCAL CS"
                if coordinate_system: # User defined cs -> Fixed pressure
                    csys = coordinate_system
                else:
                    csys = self.model.add_coordinate_system(show_in_model_tree=False)
                
                real = attributes.Real(ignore_at_execution=True) # Create empty attribute
                mat = attributes.Mat(ignore_at_execution=True) # Create empty attribute
                surf = self.model._add_surface(nset=nset,  
                                               surface_type="SURFACE ELEMENTS",
                                               etype=surface_element,
                                               real_constants=real,
                                               material=mat,
                                               coordinate_system=csys)
            else: csys = surf.csys

            force = Force(_nset, value, direction, csys, applied_by, operation, fixed)
            self.forces.append(force)
            
    def moment(self, nset : str | int, 
                     value : float,
                     direction : Literal["RX", "RY", "RZ"], 
                     behavior : Literal["Rigid", "Deformable", "Coupled", "Beam"] = "Rigid",
                     operation : Literal["ADD", "NEW", "DELETE"] = "NEW", 
                     fixed : bool = False):
        
        # Create coupling elements if moment is applied to a node set
        # Use either Rigid, Deformable, Coupled, Beam choices
        if direction in ("RX", "RY", "RZ"):
            self.moments.append((nset, direction, value, operation, fixed))

    @overload
    def pressure(self, 
                 nset : str, 
                 *,
                 operation : Literal["DELETE"] = "DELETE"): ...
    @overload
    def pressure(self, 
                 nset : str,  
                 value : float,
                 direction : Literal["X", "Y", "Z", "NORMAL TO"] = "NORMAL TO",
                 coordinate_system : CoordinateSystem | None = None,
                 applied_by : Literal["ELEMENT FACES", "SURFACE ELEMENTS"] = "ELEMENT FACES",
                 loaded_area : Literal["INITIAL", "DEFORMED"] = "DEFORMED",
                 operation : Literal["ADD", "NEW"] = "NEW"): ...
    def pressure(self, 
                 nset : str,  
                 value : float | None = None,
                 direction : Literal["X", "Y", "Z", "NORMAL TO"] = "NORMAL TO",
                 coordinate_system : CoordinateSystem | None = None,
                 applied_by : Literal["ELEMENT FACES", "SURFACE ELEMENTS"] = "ELEMENT FACES",
                 loaded_area : Literal["INITIAL", "DEFORMED"] = "DEFORMED",
                 operation : Literal["ADD", "NEW", "DELETE"] = "NEW", **kwargs):    
        
        # Check whether nset is exist.
        _nset = self.model.get_nset(nset)
        if not _nset:
            raise ParameterError(f"Nset {nset} not found!")

        if applied_by == "ELEMENT FACES":
            # SFE (through "Element Faces")
            
            if direction in ("X", "Y", "Z"):
                if coordinate_system:
                    csys = coordinate_system
                else:
                    csys = self.model.add_coordinate_system(show_in_model_tree=False)
            else: csys = None

            # If surface is already exist, use is. If not, create new one.
            surf = self.model._get_surface(nset, applied_by)
            if not surf:
                surf = self.model._add_surface(nset=nset, 
                                               surface_type="ELEMENT FACES",
                                               ignore_at_execution=True)

            pressure = Pressure(surf, value, direction, csys, applied_by, operation, loaded_area)
            self.pressures.append(pressure)

        elif applied_by == "SURFACE ELEMENTS": 

            # If surface is already exist, use is. If not, create new one.
            surf = self.model._get_surface(nset, "SURFACE ELEMENTS")
            if not surf:
                # Element Type Selection::
                if kwargs and kwargs["element_type"]:
                    # Possible Types:
                    if kwargs["element_type"] == "Surf154":
                        surface_element = element.Surf154()
                    elif kwargs["element_type"] == "Surf153":
                        surface_element = element.Surf153()
                    elif kwargs["element_type"] == "Surf156":
                        surface_element = element.Surf156()
                    elif kwargs["element_type"] == "Surf159":
                        surface_element = element.Surf159()
                    elif kwargs["element_type"] == "Follw201":
                        surface_element = element.Follw201()
                    else: surface_element = element.Surf154() # Default
                else:
                    surface_element = element.Surf154() # Default

                if loaded_area == "INITIAL":
                    surface_element.large_deflection_area = "ORIGINAL AREA"
                else:
                    surface_element.large_deflection_area = "USE NEW AREA"

                surface_element.midside_nodes = "INCLUDE"

                if direction in ("X", "Y", "Z"):
                    surface_element.pressure_cs = "LOCAL CS"
                    if coordinate_system:
                        csys = coordinate_system
                    else:
                        csys = self.model.add_coordinate_system(show_in_model_tree=False)
                else:
                    surface_element.pressure_cs = "ELEMENT CS"
                    csys = None

                real = attributes.Real(ignore_at_execution=True) # Create empty attribute
                mat = attributes.Mat(ignore_at_execution=True) # Create empty attribute
                surf = self.model._add_surface(nset=nset,  
                                               surface_type="SURFACE ELEMENTS",
                                               etype=surface_element,
                                               real_constants=real,
                                               material=mat,
                                               coordinate_system=csys)    
            else: csys = surf.csys
       
            pressure = Pressure(surf, value, direction, csys, applied_by, operation, loaded_area)
            self.pressures.append(pressure)
        
    def displacement(self, nset : str, 
                     value : float | None = None,
                     direction : Literal["ALL", "X", "Y", "Z", "RX", "RY", "RZ"] = "ALL", 
                     local_cs : CoordinateSystem | None = None,
                     operation : Literal["ADD", "NEW", "DELETE"] = "NEW",
                     ramping : bool = False):
        
        # Validate the nset:
        if not (_nset := self.model.get_nset(nset)):
            raise ParameterError(f"Nset {nset} not found!")

        cs_id = 0
        if local_cs: cs_id = local_cs.midx
        
        displacement = Displacement(set=_nset, 
                                    value=value,
                                    direction=direction,
                                    csys = cs_id,
                                    operation=operation,
                                    ramping=ramping)
        
        self.displacements.append(displacement)       

    def pretension(self):
        pass

    # Model Changes
    def activate(self):
        pass    # EALIVE

    def deactivate(self):
        pass    # EKILL

    def stabilization(self):
        pass

    def contact_stabilization(self):
        pass

    def change_friction(self):
        pass

    def output(self, item : Literal["ALL", "BASIC", "ERASE", "NODAL DOF", "REACTION LOADS",
                                    "ALL EL-SOL", "EL-NODAL LOAD", "EL-NODAL STRESS", "EL-ELASTIC STRAINS", "EL-PLASTIC STRAINS",
                                    "EL-CREEP STRAINS", "EL-NODAL GRADIENTS",
                                    "EL-INTEGRATION POINT LOCATIONS", "EL-ENERGIES", "EL-MISCELLANEOUS", "EL-NONLINEAR DATA",
                                    "EL-SURFACE STRESSES", "EL-CONTACT DATA", "EL-BACKSTRESSES", "EL-EULER ANGLES",
                                    "ALL NODAL-AVG-SOL", "NODAL-AVG STRESSES", "NODAL-AVG ELASTIC STRAINS", 
                                    "NODAL-AVG PLASTIC STRAINS", "NODAL-AVG CREEP STRAINS"],
                                    freq : int | Literal["NONE", "ALL", "LAST"] = "ALL", 
                                    set : str = ""):
        
        self.outputs.append((item, freq, set))

    def restart(self, frequency : Literal["NONE", "LAST"] | int = "LAST"):
        self.restart_controls.frequency = frequency

@dataclass
class StepControls:
    
    auto_time_stepping : Literal["ON", "OFF"] = "ON" # Both time step prediction and time step bisection will be used:
                                                     # Activates Predictor-Corrector Option automatically.
                                                     # Sets Line Search option auto
                                                     # Cutback control active
    define_by : Literal["Time", "Substeps"] = "Time"
    carry_over_time_steps : Literal["ON", "OFF"] = "OFF"
    step_end_time : float = 1
    initial_time_step : float = 1E-2
    minimum_time_step : float = 1E-5
    maximum_time_step : float = 1
    time_step : float = 1E-2
    initial_substep : int = 10
    minimum_substeps : int = 10
    maximum_substeps : int = 1E5
    number_of_substeps : int = 100

@dataclass
class SolverControls:

    solver_type : Literal["DIRECT", "ITERATIVE"] = "DIRECT"
    #weak_springs : Literal["ON", "OFF", "Program Controlled"]
    pivot_check : Literal["AUTO", "ERROR", "WARN", "OFF"] = "WARN"  #PIVCHECK,KEY,PRNTCNTRL
    pivot_check_print : Literal["ONCE", "EVERY"] = "EVERY"
    #inertia_relief : Literal["ON", "OFF"] = "OFF"
    #quasi_static : Literal["ON", "OFF"] = "OFF"

@dataclass
class NewtonRaphson:
    option : Literal["AUTO", "FULL", "MODI", "INIT", "UNSYM"] = "AUTO"
    adaptive_descent : bool | None = None
    creep_limiting : bool = False
    limit : float | None = None

@dataclass
class Convergence:
    active : bool = True # When false, solver decides ?
    value : float = ""
    tolerance : float = ""    # Must be between 0 and 1
    norm : int = ""
    min_reference : float = ""

@dataclass
class Stabilization:
    """
    For the energy dissipation ratio, specify VALUE = 1.0e-4 if you have no prior experience with the current
    model; if convergence problems are still an issue, increase the value gradually. The damping factor is
    mesh-, material-, and time-step-dependent; an initial reference value from the previous run (such as a
    run with the energy-dissipation ratio as input) should suggest itself.
    """
    state : bool = False
    control : Literal["CONSTANT", "REDUCE"] = "REDUCE"
    method : Literal["ENERGY", "DAMPING"] = "ENERGY"
    energy_dissipation_ratio : float = 1E-4
    damping_factor : float = "" # default : 1E-2
    first_step_activation : Literal["NO", "MINTIME", "ANYTIME"] = "NO"
    force_limit : float = 0.2 # To omit assign 0
   
@dataclass
class NonlinearControls:
    
    newton_raphson = NewtonRaphson(option="AUTO")
    force_convergence = Convergence(active=True, value=0.005)
    moment_convergence = Convergence(active=False, value=0.005)
    displacement_convergence = Convergence(active=True, value=0.05)
    rotation_convergence = Convergence(active=False, value=0.05)
    line_search : bool | None = None # When None, it is in auto mode
                                     # When Autots is on, it decides if line search active or not.   
    stabilization = Stabilization()
    # cutback_criteria

@dataclass
class RestartControls:
    action : Literal["DEFINE", "NORESTART", "LINEAR", "DELETE"] = "DEFINE"
    load_step : Literal["ALL", "LAST", "NONE"] | int = 1 # loadstep frequency
    frequency : Literal["NONE", "LAST"] | int = "LAST"  # substep frequency 
                                                        # negative integer for N equally spaced restart points
                                                        # positive integer for every Nth substep of a loadstep

"""@dataclass
class OutputControls:
    
    store_results : Literal["ALL TIME POINTS", "LAST TIME POINT", "EQUALLY SPACED", "SPECIFIED"]
    stress : Literal["ON", "OFF"] = "ON"
    surface_stress : Literal["ON", "OFF"] = "OFF"
    back_stress : Literal["ON", "OFF"] = "OFF"
    strain : Literal["ON", "OFF"] = "ON"
    contact_data : Literal["ON", "OFF"] = "ON"
    nonlinear_data : Literal["ON", "OFF"] = "ON"
    nodal_forces : Literal["ON", "OFF"] = "OFF"
    volume_and_energy : Literal["ON", "OFF"] = "ON"
    euler_angles : Literal["ON", "OFF"] = "OFF"
    general_miscellaneous : Literal["ON", "OFF"] = "OFF"
    contact_miscellaneous : Literal["ON", "OFF"] = "OFF"
    value : int = 1"""

@dataclass
class IterationMonitoring:
    """
    NLDIAG,Label,Key,MAXFILE

    Sets nonlinear diagnostics functionality
    writes iteration outputs (residuals, contact etc) to the files.

    Labels:
    NRRE :: Store the Newton-Raphson residuals information.
    EFLG :: Identify or display elements or nodes that violate the criteria.
    CONT :: Write contact information to a single Jobname.cnd diagnostic text file during solution.
 
    Diagnostic function characteristics:
    OFF or 0
    Suppresses writing of diagnostic information (default).
    ON or 1
    Writes diagnostic information to the Jobname.ndxxx, Jobname.nrxxx, or Jobname.cnd
    file. (If Label = CONT, this option is the same as the SUBS option described below.)
    ITER
    Writes contact diagnostic information at each iteration. Valid only when Label = CONT.
    SUBS
    Writes contact diagnostic information at each substep. Valid only when Label = CONT.
    LSTP
    Writes contact diagnostic information at each load step. Valid only when Label = CONT.
    STAT
    Lists information about the diagnostic files in the current working directory.
    DEL
    Deletes all diagnostic files in the current working directory.
    """
    function : Literal["NRRE", "EFLG", "CONT"] = "NRRE"
    characteristics : Literal["OFF", "ON", "ITER", "SUBS", "LSTP", "STAT", "DEL"] = "ON"
    maxfile : int = 4

"""# Include NLHIST
The NLHIST command is a nonlinear diagnostics tool that enables you to monitor diagnostics results
of interest in real time during a solution."""




    







