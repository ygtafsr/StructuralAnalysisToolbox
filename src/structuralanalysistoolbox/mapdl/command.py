
from __future__ import annotations
import functools
import ansys.mapdl.core as core
import numpy as np
from enum import Enum

from structuralanalysistoolbox.materials import material

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from structuralanalysistoolbox.mapdl.mesh import Mesh

    from structuralanalysistoolbox import model
    from structuralanalysistoolbox.constraints import constraint
    from structuralanalysistoolbox.loading import loadstep

def prep7(func):
    @functools.wraps(func)
    def wrapper_prep7(*args, **kwargs):  
        args[0].prep7()
        val = func(*args, **kwargs)
        args[0].finish()
        return val
    return wrapper_prep7

def solu(func):
    @functools.wraps(func)
    def wrapper_solu(*args, **kwargs):  
        args[0].slashsolu()()
        val = func(*args, **kwargs)
        args[0].finish()
        return val
    return wrapper_solu

@prep7
def _material(mapdl : core.Mapdl, mat : material.Material) -> None:
    """Runs mapdl material definition commands"""
    
    for model in mat.material_models.values():
        if isinstance(model, material.Physical):
            # Density
            if isinstance(model.density, (float, int)):
                mapdl.mp(lab="DENS", mat=f"{mat.midx}", c0=f"{model.density}")
            else:  
                label = "density"
                mapdl.load_table(name=f"{label}_{mat.midx}", array=np.array(model.density), var1="TEMP")
                mapdl.mp(lab="DENS", mat=f"{mat.midx}", c0=f"%{label}_{mat.midx}%")
        elif isinstance(model, material.IsotropicElasticity):
            # Elastic Modulus
            if isinstance(model.elastic_modulus, (float, int)):
                mapdl.mp(lab="EX", mat=f"{mat.midx}", c0=f"{model.elastic_modulus}")
            else:
                label = "elastic_modulus"
                mapdl.load_table(name=f"{label}_{mat.midx}", array=np.array(model.elastic_modulus), var1="TEMP")
                mapdl.mp(lab="EX", mat=f"{mat.midx}", c0=f"%{label}_{mat.midx}%")
            # Poisson's Ratio
            if isinstance(model.poisson_ratio, (float, int)):
                mapdl.mp(lab="NUXY", mat=f"{mat.midx}", c0=f"{model.poisson_ratio}")
            else:
                label = "poissons_ratio"
                mapdl.load_table(name=f"{label}_{mat.midx}", array=np.array(model.poisson_ratio), var1="TEMP")
                mapdl.mp(lab="NUXY", mat=f"{mat.midx}", c0=f"%{label}_{mat.midx}%")
        elif isinstance(model, material.MultilinearIsotropicHardening):
            # Multilinear Isotropic Hardening
            if model.temperature == None:
                data_count = len(model.data_table)
                mapdl.tb(lab="PLASTIC", mat=f"{mat.midx}", ntemp="1", 
                                                            npts=f"{data_count}")
                mapdl.tbtemp(temp="")
                for data in model.data_table:
                    mapdl.tbpt(oper="DEFI", x1=f"{data[0]}", x2=f"{data[1]}")
            elif isinstance(model.temperature, (float, int)):
                data_count = len(model.data_table)
                mapdl.tb(lab="PLASTIC", mat=f"{mat.midx}", ntemp="1", 
                                                            npts=f"{data_count}")
                mapdl.tbtemp(temp=f"{model.temperature}")
                for data in model.data_table:
                    mapdl.tbpt(oper="DEFI", x1=f"{data[0]}", x2=f"{data[1]}")


############################
## GET 
############################ 

# GET minimum node number for a given node set
# GET number of nodes currently selected
# GET next number after a given one in a set ??

class Get(Enum):

    Max_Node_Number = 1
    Min_Node_Number = 2


def _get(mapdl : core.Mapdl, set : model.set, value : Get):
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
def _MPC184(mapdl : core.Mapdl, mpc) -> None:

    if mpc.method == "Direct Elimination":
        method = 0
    elif mpc.method == "Lagrange Multiplier":
        method = 1
    elif mpc.method == "Penalty-based":
        method = 2

    mapdl.et(itype=mpc.etype.midx,
             ename="184",
             kop1=mpc.behaviour,
             kop2=method)
    
    mapdl.type(mpc.etype.midx)

    for node in mpc.dependent:
        mapdl.e(node, mpc.independent[0])

#@prep7
def _create_joint(mapdl : core.Mapdl, model_item_obj : constraint.MPCJoint):
    
    print(model_item_obj)


@prep7
def _execute_mesh_data(mapdl : core.Mapdl, model_item_obj : Mesh):
    """Reads the mesh model into mapdl session."""
    arv_file_path = model_item_obj.archieve.pathlib_filename
    if  arv_file_path.exists():
        mapdl.cdread(option="DB", fname=arv_file_path.absolute())


@solu
def _load_step(mapdl: core.Mapdl, load_step : loadstep.LoadStep):

    mapdl.allsel() # !!!!

    # Solver Controls 
    # Solver Type
    if load_step.solver_controls.solver_type == "DIRECT":
        mapdl.eqslv(lab="SPARSE")
    elif load_step.solver_controls.solver_type == "ITERATIVE":
        mapdl.eqslv(lab="PCG")
    else: return # Raise "InvalidParameter" exception

    # Pivot Check
    mapdl.pivcheck(key=load_step.solver_controls.pivot_check,
                   prntcntrl=load_step.solver_controls.pivot_check_print)
    
    # Large Deflection
    mapdl.nlgeom(key=load_step.solver_controls.large_deflection)

    # Weak Springs
    # Inertia Relif
    # Quasi-Static Analysis

    # Specify the analysis type and restart status
    if load_step.status == "RESTART":
        
        mapdl.antype(antype=load_step.analysis,
                    status=load_step.status,
                    ldstep=load_step.step,
                    substep=load_step.substep,
                    action=load_step.action)
    elif load_step.status == "NEW":
        mapdl.antype(antype=load_step.analysis)
    else: return # Raise an exception

    # Specify non-linear geometry option
    mapdl.nlgeom(load_step.nonlinear_geo)

    # Load Step Time
    mapdl.time(load_step.step_controls.step_end_time)

    # Time Stepping
    mapdl.autots(load_step.step_controls.auto_time_stepping)
    if load_step.step_controls.auto_time_stepping == "ON": # AUTO TIME STEPPING
        if load_step.step_controls.define_by == "Time":
            if load_step.step_controls.carry_over_time_steps == "ON":
                mapdl.deltim(dtmin=load_step.step_controls.minimum_time_step,
                            dtmax=load_step.step_controls.maximum_time_step,
                            carry=load_step.step_controls.carry_over_time_steps)
            elif load_step.step_controls.carry_over_time_steps == "OFF":
                mapdl.deltim(dtime=load_step.step_controls.initial_time_step,
                            dtmin=load_step.step_controls.minimum_time_step,
                            dtmax=load_step.step_controls.maximum_time_step,
                            carry=load_step.step_controls.carry_over_time_steps)
            else: return # Raise an exception
        elif load_step.step_controls.define_by == "Substeps":
            if load_step.step_controls.carry_over_time_steps == "ON":
                mapdl.nsubst(nsbmn=load_step.step_controls.minimum_substeps,
                            nsbmx=load_step.step_controls.maximum_substeps,
                            carry=load_step.step_controls.carry_over_time_steps)
            elif load_step.step_controls.carry_over_time_steps == "OFF":
                mapdl.nsubst(nsbstp=load_step.step_controls.initial_substep,
                            nsbmn=load_step.step_controls.minimum_substeps,
                            nsbmx=load_step.step_controls.maximum_substeps,
                            carry=load_step.step_controls.carry_over_time_steps)
            else: return # Raise an exception
        else: return # Raise an exception
    elif load_step.step_controls.auto_time_stepping == "OFF": # FIXED TIME STEPPING
            mapdl.deltim(dtime=load_step.step_controls.time_step)
            mapdl.nsubst(nsbstp=load_step.step_controls.number_of_substeps)
    else: return # raise an exception

    # Loading Style
    if load_step.loading_style == "RAMPED":
        mapdl.kbc(key=0)
    elif load_step.loading_style == "STEPPED":
        mapdl.kbc(key=1)

    # Reference Temperature
    mapdl.tref(tref=load_step.reference_temp)

    ########################
    ### Nonlinear Controls
    ########################

    # Newton-Raphson Option
    if load_step.nonlinear_controls.newton_raphson.creep_limiting:
        if load_step.nonlinear_controls.newton_raphson.limit:
            mapdl.nropt(option1=load_step.nonlinear_controls.newton_raphson.option,
                        option2="CRPL",
                        optval=load_step.nonlinear_controls.newton_raphson.limit)
        else:
            mapdl.nropt(option1=load_step.nonlinear_controls.newton_raphson.option,
                        option2="CRPL")
    else:
        if load_step.nonlinear_controls.newton_raphson.adaptive_descent == True: # Adaptive descent is active
            mapdl.nropt(option1=load_step.nonlinear_controls.newton_raphson.option,
                        option2="",
                        optval="ON")
        elif load_step.nonlinear_controls.newton_raphson.adaptive_descent == False:# Adaptive descent is deactive
            mapdl.nropt(option1=load_step.nonlinear_controls.newton_raphson.option,
                        option2="",
                        optval="OFF")
        elif load_step.nonlinear_controls.newton_raphson.adaptive_descent == None: # Default
            mapdl.nropt(option1=load_step.nonlinear_controls.newton_raphson.option)


    # Convergence Controls
    # Force Convergence
    if load_step.nonlinear_controls.force_convergence.active:
        mapdl.cnvtol(lab="F",
                     value=load_step.nonlinear_controls.force_convergence.value,
                     toler=load_step.nonlinear_controls.force_convergence.tolerance,
                     norm=load_step.nonlinear_controls.force_convergence.norm,
                     minref=load_step.nonlinear_controls.force_convergence.min_reference)
    # Moment Convergence
    if load_step.nonlinear_controls.moment_convergence.active:
        mapdl.cnvtol(lab="M",
                     value=load_step.nonlinear_controls.moment_convergence.value,
                     toler=load_step.nonlinear_controls.moment_convergence.tolerance,
                     norm=load_step.nonlinear_controls.moment_convergence.norm,
                     minref=load_step.nonlinear_controls.moment_convergence.min_reference)
    # Volume Convergence
    # mapdl.cnvtol(lab="DVOL")
    # Displacement Convergence
    if load_step.nonlinear_controls.displacement_convergence.active:
        mapdl.cnvtol(lab="U",
                     value=load_step.nonlinear_controls.displacement_convergence.value,
                     toler=load_step.nonlinear_controls.displacement_convergence.tolerance,
                     norm=load_step.nonlinear_controls.displacement_convergence.norm,
                     minref=load_step.nonlinear_controls.displacement_convergence.min_reference)
    # Rotation Convergence
    if load_step.nonlinear_controls.rotation_convergence.active:
        mapdl.cnvtol(lab="ROT",
                     value=load_step.nonlinear_controls.rotation_convergence.value,
                     toler=load_step.nonlinear_controls.rotation_convergence.tolerance,
                     norm=load_step.nonlinear_controls.rotation_convergence.norm,
                     minref=load_step.nonlinear_controls.rotation_convergence.min_reference)
    # Hydrostatic Pressure
    # mapdl.cnvtol(lab="HDSP")
    # Joint Element Constraint Check
    # mapdl.cnvtol(lab="JOINT")
    # Volumetric Compability Check
    # mapdl.cnvtol(lab="COMP")

    # Line Search Option
    if load_step.nonlinear_controls.line_search:
        mapdl.lnsrch(key="ON")
    elif not load_step.nonlinear_controls.line_search:
        mapdl.lnsrch(key="OFF")
    elif load_step.nonlinear_controls.line_search == None:
        mapdl.lnsrch(key="AUTO")

    # Cutback Controls


    # Iteration Monitoring
    mapdl.nldiag(label=load_step.resiudal_monitoring.function,
                 key=load_step.resiudal_monitoring.characteristics,
                 maxfile=load_step.resiudal_monitoring.maxfile)
    
    mapdl.nldiag(label=load_step.violation_monitoring.function,
                 key=load_step.violation_monitoring.characteristics,
                 maxfile=load_step.violation_monitoring.maxfile)
    
    mapdl.nldiag(label=load_step.contact_monitoring.function,
                 key=load_step.contact_monitoring.characteristics,
                 maxfile=load_step.contact_monitoring.maxfile)
    
    # Stabilization
    if load_step.nonlinear_controls.stabilization.state:
        if load_step.nonlinear_controls.stabilization.method == "ENERGY":
            mapdl.stabilize(key=load_step.nonlinear_controls.stabilization.control,
                            method=load_step.nonlinear_controls.stabilization.method,
                            value=load_step.nonlinear_controls.stabilization.energy_dissipation_ratio,
                            substpopt=load_step.nonlinear_controls.stabilization.first_step_activation,
                            forcelimit=load_step.nonlinear_controls.stabilization.force_limit)
        elif load_step.nonlinear_controls.stabilization.method == "DAMPING":
            mapdl.stabilize(key=load_step.nonlinear_controls.stabilization.control,
                            method=load_step.nonlinear_controls.stabilization.method,
                            value=load_step.nonlinear_controls.stabilization.damping_factor,
                            substpopt=load_step.nonlinear_controls.stabilization.first_step_activation,
                            forcelimit=load_step.nonlinear_controls.stabilization.force_limit)
    else:
        mapdl.stabilize(key="OFF")


    
    




            









