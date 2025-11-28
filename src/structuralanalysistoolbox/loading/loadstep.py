from dataclasses import dataclass
from typing import Literal

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

class LoadStep:
    """
    Load Steps defined in the "solu" processor.
    KEY COMMANDS:
    TIME : Time at the end of a load step.
    KBC : Ramped/Stepped Loading within a load step.
    DELTIM : Time step size.
    NSUBST : Number of Substeps within a load step.
    AUTOTS : Automatic time stepping
    """
    def __init__(self, step_number : int, 
                 analysis : Literal["STATIC", "BUCKLE", "MODAL", "HARMIC" ,"TRANS", "SUBSTR"] = "STATIC",
                 status : Literal["NEW", "RESTART"] = "NEW",
                 step : int = 0, substep : int = 0):
        
        self.step_number = step_number
        self.analysis = analysis
        self.status = status
        self.step = step # For RESTART
        self.substep = substep # For RESTART
        self.action : Literal["CONTINUE", "ENDSTEP", "RSTCREATE", "PERTURB"] = "CONTINUE" # For RESTART
        #self.perturbation : Literal["YES", "NO"] = "NO" 
        self.nonlinear_geo : Literal["ON", "OFF"] = "ON" # NLGEOM can not be changed after the first SOLVE.
        self.loading_style : Literal["RAMPED", "STEPPED"]
        self.reference_temp : float = 0

        self.step_controls = StepControls()
        self.solver_controls = SolverControls()
        self.restart_controls = RestartControls()
        self.nonlinear_controls = NonlinearControls()
        self.output_controls = OutputControls()

        self.resiudal_monitoring = IterationMonitoring(function="NRRE", characteristics="ON", maxfile=4)
        self.violation_monitoring = IterationMonitoring(function="EFLG", characteristics="ON", maxfile=4)
        self.contact_monitoring = IterationMonitoring(function="CONT", characteristics="ITER", maxfile=4)

        self._output_list = []

    def add_force(self):
        pass

    def add_pressure(self):
        pass

    def add_bc(self):
        pass

    # Model Changes
    def activate(self):
        pass    # EALIVE

    def deactivate(self):
        pass    # EKILL

    def change_friction(self):
        pass

    def output(self, *args):
        for out in args:
            self._output_list.append(out)

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
    large_deflection : Literal["ON", "OFF"] = "ON"
    #inertia_relief : Literal["ON", "OFF"] = "OFF"
    #quasi_static : Literal["ON", "OFF"] = "OFF"

@dataclass
class Convergence:
    active : bool = True # When false, solver decides ?
    value : float = ""
    tolerance : float = ""    # Must be between 0 and 1
    norm : int = ""
    min_reference : float = ""

@dataclass
class NewtonRaphson:
    option : Literal["AUTO", "FULL", "MODI", "INIT", "UNSYM"] = "AUTO"
    adaptive_descent : bool | None = None
    creep_limiting : bool = False
    limit : float | None = None

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
    force_convergence = Convergence(active=True)
    moment_convergence = Convergence(active=False)
    displacement_convergence = Convergence(active=True)
    rotation_convergence = Convergence(active=False)
    line_search : bool | None = None # When None, it is in auto mode
                                     # When Autots is on, it decides if line search active or not.   
    stabilization = Stabilization()
    # cutback_criteria

@dataclass
class RestartControls:
    
    generate_restart_points : Literal["MANUEL", "OFF"] = "OFF"
    load_step : Literal["ALL", "LAST", "SPECIFY"]
    load_step_number : int
    sub_step : Literal["ALL", "LAST", "RECURRENCE RATE", "EQUALLY SPACED"]
    reccurrence_rate : int = 1 
    equal_points_value : int = 1
    max_points : int = 1 # Max points to save per step
    retain_files : Literal["YES", "NO"] = "NO"
    combine_files : Literal["YES", "NO"] = "NO"

@dataclass
class OutputControls:
    
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
    store_results : Literal["ALL TIME POINTS", "LAST TIME POINT", "EQUALLY SPACED", "SPECIFIED"]
    value : int = 1

@dataclass
class Output:

    general : Literal["ALL", "BASIC", "ERASE"]
    nodal_solution : Literal["NODAL DOF", "REACTION LOADS"]
    element_solution : Literal["ALL",
                               "NODAL LOAD",
                               "NODAL STRESS",
                               "ELASTIC STRAINS",
                               "PLASTIC STRAINS",
                               "CREEP STRAINS",
                               "DIFFUSION STRAINS",
                               "NODAL GRADIENTS",
                               "NODAL FLUXES",
                               "INTEGRATION POINT LOCATIONS",
                               "ENERGIES",
                               "MISCELLANEOUS",
                               "NONLINEAR DATA",
                               "SURFACE STRESSES",
                               "CONTACT DATA",
                               "BACKSTRESSES",
                               "EULER ANGLES"] = "ALL"
    nodal_averaged_solution : Literal["ALL",
                                      "STRESSES",
                                      "ELASTIC STRAINS",
                                      "PLASTIC STRAINS",
                                      "CREEP STRAINS"] = "ALL"
    frequency = 1
    component : str = "ALL"


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






