from __future__ import annotations
from dataclasses import dataclass
from typing import Literal

from structuralanalysistoolbox.mapdl import element
from structuralanalysistoolbox.exceptions import ParameterError


from typing import TYPE_CHECKING, overload

if TYPE_CHECKING:
    from structuralanalysistoolbox.model import Nset, Elset, Surface, LocalCoordinateSystem, Vector, Analysis, Model



"""
NONLINEAR OPTIONS:
NEQIT : Specifies the max num. of eq. iter. per substep.
CNVTOL : Specifies convergence tolerances.
NCNV : Provides options for terminating analysis.
CUTCONTROL,Lab,VALUE,Option
Controls time-step cutback during a nonlinear solution.
"""

"""
OUTPUT CONTROLS:
OUTRES : Controls what the program writes to the OUTRES database and results file
         and how often it is written.
OUTPR : Controls what is printed (written to the solution OUTPR output file, (Jobname.OUT) and how
        often it is written.
OUTGEOM : Controls the type of geometry data the program writes to the results file.
ERESX : enables you to review element integration point values in the postprocessor.
"""

"""
CREATING MULTIPLE LOAD STEPs

"""

@dataclass
class Force:
    set : Nset | Surface
    value : float
    direction : Literal["X", "Y", "Z"]
    local_cs_id : int | None = None
    applied_by : Literal["NODE SET", "ELEMENT FACES", "SURFACE ELEMENTS"] = "NODE SET",
    operation : Literal["ADD", "NEW", "DELETE"] = "NEW"
    fixed : bool = False

@dataclass
class Moment:
    set : Nset | Surface
    direction : Literal["RX", "RY", "RZ"]
    value : float
    behavior : Literal["Rigid", "Deformable", "Coupled", "Beam"] = "Rigid"
    operation : Literal["ADD", "NEW", "DELETE"] = "NEW"
    fixed : bool = False

@dataclass
class Pressure:
    surf : Surface | None
    value : float
    local_cs_id : int | None = None
    direction : Literal["X", "Y", "Z", "NORMAL TO"] = "NORMAL TO"
    applied_by : Literal["NODE SET", "ELEMENT FACES", "SURFACE ELEMENTS"] = "SURFACE ELEMENTS"
    operation : Literal["ADD", "NEW", "DELETE"] = "NEW"
    loaded_area : Literal["INITIAL", "DEFORMED"] = "DEFORMED"

@dataclass
class Displacement:
    set : Nset 
    direction : Literal["ALL", "X", "Y", "Z", "RX", "RY", "RZ"]
    value : float
    operation : Literal["ADD", "NEW", "DELETE"] = "NEW"
    ramping : bool = False

class LoadStep:
    
    def __init__(self,
                 name : str,
                 model : Model,
                 step_number : int = 0,
                 end_time : int = 0, 
                 analysis : Literal["STATIC", "BUCKLE", "MODAL", "HARMIC" ,"TRANS", "SUBSTR"] = "STATIC",
                 status : Literal["NEW", "RESTART"] = "NEW",
                 step : int = 0, substep : int = 0):
        
        self.name = name
        self.midx = 0
        self.model = model
        self.step_number = step_number
        self.analysis = analysis
        self.status = status
        self.step = step # For RESTART
        self.substep = substep # For RESTART
        self.action : Literal["CONTINUE", "ENDSTEP", "RSTCREATE", "PERTURB"] = "CONTINUE" # For RESTART
        #self.perturbation : Literal["YES", "NO"] = "NO" 
        self.nonlinear_geo : Literal["ON", "OFF"] = "ON" # NLGEOM can not be changed after the first SOLVE.
        self.loading_style : Literal["RAMPED", "STEPPED"] = "RAMPED"
        self.reference_temp : float = 0

        self.step_controls = StepControls()
        self.step_controls.step_end_time = end_time
        
        self.solver_controls = SolverControls()
        self.restart_controls = RestartControls()
        self.nonlinear_controls = NonlinearControls()

        self.resiudal_monitoring = IterationMonitoring(function="NRRE", characteristics="ON", maxfile=4)
        self.violation_monitoring = IterationMonitoring(function="EFLG", characteristics="ON", maxfile=4)
        self.contact_monitoring = IterationMonitoring(function="CONT", characteristics="ITER", maxfile=4)

        self.forces : list[Force] = []
        self.moments : list[Moment] = []
        self.pressures : list[Pressure] = []
        self.displacements : list[Displacement] = []
        self.outputs : list = []

    def force(self, nset : str, 
                    value : float,
                    direction : Literal["X", "Y", "Z"],
                    local_cs : LocalCoordinateSystem | None = None,
                    applied_by : Literal["NODE SET", "ELEMENT FACES", "SURFACE ELEMENTS"] = "NODE SET",
                    operation : Literal["ADD", "NEW", "DELETE"] = "NEW", 
                    fixed : bool = False):
        """
           NEW :: Replace previous Force, apply the new one.
           ADD :: Apply the new load onto the previous one.
        """
        
        # Validate the nset:
        _nset : Nset
        if isinstance(nset, str) and self.model.get_nset(nset): 
            _nset = self.model.get_nset(nset) # pyright: ignore[reportAssignmentType]
        elif isinstance(nset, Nset):
            _nset = nset
        else: raise ParameterError(f"Node set '{nset}' does not exist in the model!")
        
        if applied_by == "NODE SET":
            self.forces.append(Force(set=_nset, 
                                     direction=direction, 
                                     applied_by=applied_by,
                                     value=value, 
                                     operation=operation, 
                                     fixed=fixed))
        elif applied_by == "ELEMENT FACES":
            # Apply through Pressure on Element Faces
            _surf = self.model.add_surface(nset=_nset.name, name=f"Surf_Force_{_nset.name}", ignore_at_execution=True)
            _cs_id = 0
            if direction in ("X", "Y", "Z"):
                _cs_id = self.model.add_coordinate_system(origin=(0,0,0), rot_x=0, rot_y=0, rot_z=0).midx # Default CS

            self.forces.append(Force(set=_surf,
                                     value=value, 
                                     direction=direction, 
                                     local_cs_id=_cs_id,
                                     applied_by=applied_by, 
                                     operation=operation, 
                                     fixed=fixed))
        elif applied_by == "SURFACE ELEMENTS":
            # Apply through Pressure on Surface Elements
            _etype_154 = element.Surf154() 
            _surf = self.model.add_surface(nset=_nset.name, name=f"Surf_Force_{_nset.name}", 
                                           element_type=_etype_154,
                                           application="Pressure")

            if direction in ("X", "Y", "Z"):
                _etype_154.pressure_cs = "LOCAL CS"
                if local_cs:
                    _surf.local_csys_id = local_cs.midx
                else:
                    _surf.local_csys_id = self.model.add_coordinate_system(origin=(0,0,0), rot_x=0, rot_y=0, rot_z=0).midx # Default CS

            _etype_154.large_deflection_area = "ORIGINAL AREA"
            _etype_154.midside_nodes = "INCLUDE"

            self.forces.append(Force(set=_surf,
                                     value=value, 
                                     direction=direction, 
                                     local_cs_id=_surf.local_csys_id, 
                                     applied_by=applied_by,
                                     operation=operation, 
                                     fixed=fixed))
            
    def remote_force(self, 
                     pilot_node : str,
                     contact_nodes : str,
                     value : float, 
                     direction : Literal["X", "Y", "Z"],
                     applied_by : Literal["Rigid", "Deformable", "Coupled"] = "Deformable",
                     operation : Literal["ADD", "NEW", "DELETE"] = "NEW", 
                     fixed : bool = False,
                     local_cs : LocalCoordinateSystem | None = None,):
        
        if applied_by == "Rigid":
            self.model.add_rigid_surface_constraint(pilot_node=pilot_node, contact_nodes=contact_nodes)
            self.force(pilot_node, value=value, direction=direction, 
                       local_cs=local_cs, applied_by="NODE SET", operation=operation, fixed=fixed)
        elif applied_by == "Deformable":
            self.model.add_force_dist_surf_constraint(pilot_node=pilot_node, contact_nodes=contact_nodes)
            self.force(pilot_node, value=value, direction=direction, 
                       local_cs=local_cs, applied_by="NODE SET", operation=operation, fixed=fixed)
        elif applied_by == "Coupled":
            self.model.add_coupled_surface_constraint(pilot_node=pilot_node, contact_nodes=contact_nodes)
            self.force(pilot_node, value=value, direction=direction, 
                       local_cs=local_cs, applied_by="NODE SET", operation=operation, fixed=fixed)

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

    def pressure(self, nset : str | Nset,  
                 value : float,
                 direction : Literal["X", "Y", "Z", "NORMAL TO"] = "NORMAL TO",
                 local_cs : LocalCoordinateSystem | None = None,
                 applied_by : Literal["NODE SET", "ELEMENT FACES", "SURFACE ELEMENTS"] = "SURFACE ELEMENTS",  # Apply through "surface elements"
                 loaded_area : Literal["INITIAL", "DEFORMED"] = "DEFORMED",                                                                                            # Direct Apply through "element surfaces"   
                 operation : Literal["ADD", "NEW", "DELETE"] = "NEW"):    # ADD: add to the existing pressure (SFCUM::ADD)
                                                                          # NEW: overwrite the existing pressure (SFCUM::REPL)
                                                                          # DELETE: remove the pressure (SFDELE)
        # Validate the nset:
        _nset : Nset
        if isinstance(nset, str) and self.model.get_nset(nset): 
            _nset = self.model.get_nset(nset) # pyright: ignore[reportAssignmentType]
        elif isinstance(nset, Nset):
            _nset = nset
        else: raise ParameterError(f"Node set '{nset}' does not exist in the model!")

        # Create and store the pressure object
        _pressure : Pressure
        if applied_by == "NODE SET":
            _surf = self.model.add_surface(nset=_nset.name, name=f"Surf_Pres_{_nset.name}", 
                                           ignore_at_execution=True)
            _pressure = Pressure(surf=_surf, 
                                value=value, 
                                direction="NORMAL TO",
                                applied_by=applied_by, 
                                operation=operation, 
                                loaded_area="INITIAL") # Check this
            
        elif applied_by == "ELEMENT FACES":
            # SFE (through "Element Faces")
            _surf = self.model.add_surface(nset=_nset.name, name=f"Surf_Pres_{_nset.name}", 
                                           ignore_at_execution=True)
            _cs_id = 0
            if direction in ("X", "Y", "Z"):
                _cs_id = self.model.add_coordinate_system(origin=(0,0,0), rot_x=0, rot_y=0, rot_z=0).midx # Default CS
            
            _pressure = Pressure(surf=_surf, 
                                value=value, 
                                local_cs_id=_cs_id,
                                direction=direction,
                                applied_by=applied_by, 
                                operation=operation, 
                                loaded_area=loaded_area,)
            
        elif applied_by == "SURFACE ELEMENTS": 
            # SURF154 + SFE (through "Surface Elements")
            # Create surface elements and specify keyopts
            _etype_154 = element.Surf154() 
            _surf = self.model.add_surface(nset=_nset.name, name=f"Surf_Pres_{_nset.name}", 
                                           element_type=_etype_154,
                                           application="Pressure")

            if direction == "NORMAL TO":    # keyopt2=0 (use default element CS)
                _etype_154.pressure_cs = "ELEMENT CS"
                _etype_154.normal_pressure_direction = "CALCULATED"
            elif direction in ("X", "Y", "Z"): # keyopt2=1 (use local CS as Element CS)
                _etype_154.pressure_cs = "LOCAL CS"
                if local_cs:
                    _surf.local_csys_id = local_cs.midx
                else:
                    _surf.local_csys_id = self.model.add_coordinate_system(origin=(0,0,0), rot_x=0, rot_y=0, rot_z=0).midx # Default CS

            if loaded_area == "INITIAL":
                _etype_154.large_deflection_area = "ORIGINAL AREA"
            else:
                _etype_154.large_deflection_area = "USE NEW AREA"
                
            _etype_154.midside_nodes = "INCLUDE"

            _pressure = Pressure(surf=_surf, 
                                value=value, 
                                local_cs_id=_surf.local_csys_id,
                                direction=direction,
                                applied_by=applied_by, 
                                operation=operation, 
                                loaded_area=loaded_area) # defined in element type

        self.pressures.append(_pressure)
        
    def dof(self, 
            nodes : str | int, 
            value : float,
            direction : Literal["ALL", "X", "Y", "Z", "RX", "RY", "RZ"], 
            operation : Literal["ADD", "NEW", "DELETE"] = "NEW", 
            ramping : bool = False):
        
        self.displacements.append((nodes, direction, value, operation, ramping))

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
    
    load_step : Literal["ALL", "LAST", "SPECIFY"] = "LAST"
    load_step_number : int = 0
    sub_step : Literal["ALL", "LAST", "RECURRENCE RATE", "EQUALLY SPACED"] = "LAST"
    generate_restart_points : Literal["MANUEL", "OFF"] = "OFF"
    reccurrence_rate : int = 1 
    equal_points_value : int = 1
    max_points : int = 1 # Max points to save per step
    retain_files : Literal["YES", "NO"] = "NO"
    combine_files : Literal["YES", "NO"] = "NO"

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
    """
    function : Literal["NRRE", "EFLG", "CONT"] = "NRRE"
    characteristics : Literal["OFF", "ON", "ITER", "SUBS", "LSTP", "STAT", "DEL"] = "ON"
    maxfile : int = 4


"""ls = LoadStep(step_number=1)
ls.output(Output(nodal_solution="NODAL DOF"),
          Output(nodal_solution="REACTION LOADS"),
          Output(element_solution="ELASTIC STRAINS"),
          Output(element_solution="PLASTIC STRAINS"),
          Output(nodal_averaged_solution="ALL"))"""



    







