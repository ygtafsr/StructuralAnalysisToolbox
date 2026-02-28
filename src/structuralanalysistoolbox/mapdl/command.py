
from __future__ import annotations
import functools
import numpy as np
import ansys.mapdl.core as core
from structuralanalysistoolbox.materials import material

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    
    from structuralanalysistoolbox.data.mesh import Mesh
    from structuralanalysistoolbox.mapdl import element
    from structuralanalysistoolbox import model
    from structuralanalysistoolbox.model import Nset, Elset, Surface, PilotNode
    from structuralanalysistoolbox.constraints import constraint
    from structuralanalysistoolbox.loading.loadstep import LoadStep, Force, Moment, Pressure, Displacement

def prep7(func):
    @functools.wraps(func)
    def wrapper_prep7(*args, **kwargs):  
        args[0].prep7()
        val = func(*args, **kwargs)
        args[0].finish()
        return val
    return wrapper_prep7

##########################################
#### MISCALENOUS HELPER FUNCTIONS
##########################################

def is_surface_node(node_id : int) -> bool:
    "TODO: later implement for VTK"
    # Get a node
    # Find all connected elements connected to this node
    # 
    
    return True

def get_face_ids(nset : Nset) -> np.ndarray | None:
    """
    Get a node list. 
    Select nodes in the list
    Find conncected elements
    
    get elements face node list:
    MAPDL COMMANDs:

        *DIM, ENODES, ARRAY, 20 
        *VGET, ENODES(1), ELEM, 1800, NODE, 1,,,4

        For Solid 185:  I - J - K - L - M - N - O - P
            array idx:  0   1   2   3   4   5   6   7

        Face1 (J-I-L-K) (1,0,3,2)
        Face2 (I-J-N-M) (0,1,5,4)
        Face3 (J-K-O-N) (1,2,6,5)
        Face4 (K-L-P-O) (2,3,7,6)
        Face5 (L-I-M-P) (3,0,4,7)
        Face6 (M-N-O-P) (4,5,6,7)
    
    Return:: Face IDs List.
    """

    face_ids = []


    return np.array(face_ids)

def get_surface_area(mapdl : core.Mapdl, surf : Surface) -> float:
    """
    Calculates surface area of surf using "ARFACE(E)" get function.
    Function returns the area of the face of element E containing 
    selected nodes.
    """
    # Select all surface nodes
    nlist = mapdl.nsel(_type="S", vmin=surf.nset.items) # type: ignore
    elist = mapdl.esln(type_='S', ekey='0', nodetype='ALL') # type: ignore

    total_area = 0.0
    # iterate over elements, find the nodes for each and check whether face nodes in nlist
    for elem_id in elist:  
        mapdl.esel(type_='S', item='ELEM', vmin=str(elem_id)) # type: ignore
        elem_nodes = mapdl.nsle(type_='S', nodetype='ALL')
        # Find the common nodes between the element nodes and all selected list
        intersection = list(set(elem_nodes) & set(nlist))
        if len(intersection) >= 3: # Minimum 3 nodes to define a face
            with mapdl.non_interactive:
                mapdl.run(f"_param = ARFACE({elem_id})")
            area = float(mapdl.parameters["_param"])
            total_area += area

    mapdl.allsel() # Reset selection after force application
    return total_area

@prep7
def _reset_attributes(mapdl : core.Mapdl):

    # Reset Attributes to the defaults
    mapdl.csys(0)
    mapdl.real(1)
    mapdl.mat(1)
    mapdl.type(1)

#################################################################################

def _local_csys(mapdl : core.Mapdl, local_csys : model.CoordinateSystem) -> None:
    """Creates a local coordinate system in MAPDL."""
    mapdl.local(kcn=str(local_csys.midx), kcs=local_csys.type, 
                xc=local_csys.origin[0], yc=local_csys.origin[1], zc=local_csys.origin[2],
                thxy=str(local_csys.rot_x), thyz=str(local_csys.rot_y), thzx=str(local_csys.rot_z))
    
    # reset to default csys
    mapdl.csys(0)

@prep7
def _create_pilot_node(mapdl: core.Mapdl, *args):

    pilot_node : model.PilotNode = args[0]
    
    mapdl.csys(0)
    if pilot_node.local_cs:
        mapdl.csys(pilot_node.local_cs.midx)

    # Create a node
    pilot_node_number = mapdl.n(x=str(pilot_node.x), y=str(pilot_node.y), z=str(pilot_node.z))

    # Create a component with the name same as pilot node name
    mapdl.components[pilot_node.nset.name] = "NODE", [pilot_node_number] # type: ignore

def _create_node_set(mapdl : core.Mapdl, nset : Nset):
    mapdl.components[nset.name] = "NODE", nset.items # type: ignore

@prep7
def _create_element_set(mapdl : core.Mapdl, elset : Elset):

    if elset.name and elset.items:
        mapdl.components[elset.name] = "ELEM", elset.items # type: ignore

@prep7
def _create_surface(mapdl : core.Mapdl, *args):
        
    surf : model.Surface = args[0]

    mapdl.type(str(surf.element_type.midx))  # type: ignore
    mapdl.real(str(surf.real_constants.midx)) # type: ignore
    mapdl.mat(str(surf.material.midx)) # type: ignore
    mapdl.esys(str(surf.csys.midx)) # type : ignore

    mapdl.nsel(type_= 'S', vmin = surf.nset.items) # type: ignore
    mapdl.esurf()

    _reset_attributes(mapdl)
    mapdl.allsel()

##########################################
#### ELEMENT ATTRIBUTES
##########################################

@prep7
def _create_element_type(mapdl : core.Mapdl, *args):
    """Create the element type & Asssign Keyopts"""
    etype = args[0]
    mapdl.et(etype.midx, ename=etype.name)

    for keyopt in etype.get_keyopts(): # keyopt :: (keyopt_id, value)
        mapdl.keyopt(itype= str(etype.midx), knum=keyopt[0], value=keyopt[1])

@prep7
def _create_real_constants(mapdl : core.Mapdl, rc):
    """Gets a real constant object and creates real constants in MAPDL."""

    real_constants = rc.get_real_constants()

    for idx, rc_group in enumerate(real_constants):
        if idx == 0:
            mapdl.r(f"{rc.midx}", # real constant midx
                    rc_group[0], rc_group[1], rc_group[2], rc_group[3], rc_group[4], rc_group[5])
        else:
            mapdl.rmore(rc_group[0], rc_group[1], rc_group[2], rc_group[3], rc_group[4], rc_group[5])
    
def _update_real_constants(mapdl : core.Mapdl, *args):
    pass

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

@prep7
def _assign_material(mapdl : core.Mapdl, elset : Elset):

    if elset.material:
        mapdl.emodif(iel=elset.name, stloc="MAT", i1=elset.material.midx)
        mapdl.allsel()

@prep7
def _assign_element_type(mapdl : core.Mapdl, elset : Elset):

    if elset.element_type:
        mapdl.emodif(iel=elset.name, stloc="TYPE", i1=elset.element_type.midx)
        mapdl.allsel()

@prep7
def _assign_real_constants(mapdl : core.Mapdl, elset : Elset):

    if elset.real_constants: 
        mapdl.emodif(iel=elset.name, stloc="REAL", i1=elset.real_constants.midx)
        mapdl.allsel()


#############################################
#### CONSTRAINTS
#############################################

## SURFACE BASED CONSTRAINTS ##
@prep7
def _create_rigid_surface_constraint(mapdl : core.Mapdl, *args):
        
    rigid_surface_constraint : constraint.RigidSurfaceConstraint = args[0]
    contact = rigid_surface_constraint.contact_element_type
    target = rigid_surface_constraint.target_element_type
    pilot_node = mapdl.components[rigid_surface_constraint.pilot_node]
    contact_nodes = mapdl.components[rigid_surface_constraint.contact_nodes]

    # Activate associated real constant and material
    mapdl.real(str(rigid_surface_constraint.real_constants.midx))
    mapdl.mat(str(rigid_surface_constraint.material_midx))

    # Local cs
    if rigid_surface_constraint.local_cs:
        mapdl.csys(rigid_surface_constraint.local_cs.midx)

    # Create contact elements
    mapdl.type(str(contact.midx))
    for node in contact_nodes:
        mapdl.e(node)

    # Create pilot node element
    mapdl.type(str(target.midx))
    # Define target segment elements
    mapdl.tshap(shape="PILO") 
    mapdl.e(pilot_node[0])
    mapdl.tshap()

    _reset_attributes(mapdl)
    mapdl.allsel()

@prep7
def _create_force_distirbuted_constraint(mapdl : core.Mapdl, *args):
    
    force_distributed_constraint : constraint.ForceDistributedConstraint = args[0]
    contact = force_distributed_constraint.contact_element_type
    target = force_distributed_constraint.target_element_type
    pilot_node = mapdl.components[force_distributed_constraint.pilot_node]
    contact_nodes = mapdl.components[force_distributed_constraint.contact_nodes]

    # First modify real constant fkn, if there is a "user defined" weightening
    if isinstance(force_distributed_constraint.weighting_factor, np.ndarray):
        # For a naming convention get first letters of pilot and contact node set names
        name = f"{force_distributed_constraint.pilot_node[:4]}_{force_distributed_constraint.contact_nodes[:4]}"
        mapdl.load_table(name=name,
                         array=force_distributed_constraint.weighting_factor,
                         var1="NODE")
        mapdl.rmodif(nset=str(force_distributed_constraint.real_constants.midx), stloc=3, v1=f"%{name}%")

    # Activate associated attributes
    mapdl.real(str(force_distributed_constraint.real_constants.midx))
    mapdl.mat(str(force_distributed_constraint.material.midx))
    mapdl.csys(force_distributed_constraint.coordinate_system.midx)

    # Create contact elements
    mapdl.type(str(contact.midx))
    for node in contact_nodes:
        mapdl.e(node)

    # Create pilot node element
    mapdl.type(str(target.midx))
    # Define target segment elements
    mapdl.tshap(shape="PILO") 
    mapdl.e(pilot_node[0])
    mapdl.tshap()

    _reset_attributes(mapdl)
    mapdl.allsel()

@prep7
def _create_coupled_surface_constraint(mapdl : core.Mapdl, *args):

    coupled_surface_constraint : constraint.CoupledSurfaceConstraint = args[0]
    contact = coupled_surface_constraint.contact_element_type
    target = coupled_surface_constraint.target_element_type
    pilot_node = mapdl.components[coupled_surface_constraint.pilot_node]
    contact_nodes = mapdl.components[coupled_surface_constraint.contact_nodes]

    # Activate associated real constant and material
    mapdl.real(str(coupled_surface_constraint.real_constants.midx))
    mapdl.mat(str(coupled_surface_constraint.material_midx))

    # Local cs
    if coupled_surface_constraint.local_cs:
        mapdl.csys(coupled_surface_constraint.local_cs.midx)

    # Create contact elements
    mapdl.type(str(contact.midx))
    for node in contact_nodes:
        mapdl.e(node)

    # Create pilot node element
    mapdl.type(str(target.midx))
    # Define target segment elements
    mapdl.tshap(shape="PILO") 
    mapdl.e(pilot_node[0])
    mapdl.tshap()

    _reset_attributes(mapdl)
    mapdl.allsel()

## RIGID MPC CONSTRAINTS ##
@prep7                           
def _MPC184(mapdl : core.Mapdl, *args) -> None:

    mpc = args[0]
    method = 0
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

    for node in mpc.dependent.items:
        mapdl.e(node, mpc.independent.items[0])

## LINEAR COUPLINGS ##
@prep7
def _coupled_dof(mapdl : core.Mapdl, *args) -> int:
    """Create a dof coupling between nodes and returns primary node number 
    as an integer for coupled set."""
    coupling : constraint.CoupledDOF = args[0]
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
def _create_joint(mapdl : core.Mapdl, *args):
    model_item_obj : constraint.MPCJoint = args[0]
    print(model_item_obj)

#############################################
#### IMPORTING EXTERNAL MESH DATA
#############################################

@prep7
def _execute_mesh_data(mapdl : core.Mapdl, *args):
    """Reads the mesh model into mapdl session."""
    model_item_obj : Mesh = args[0]
    arv_file_path = model_item_obj.archieve.pathlib_filename
    if  arv_file_path.exists():
        mapdl.cdread(option="DB", fname=arv_file_path.absolute())

output_mapping = {
    "ALL" : "ALL",
    "BASIC" : "BASIC",
    "ERASE" : "ERASE",
    "NODAL DOF" : "NSOL",
    "REACTION LOADS" : "RSOL",
    "ALL EL-SOL" : "ESOL",
    "EL-NODAL LOAD" : "NLOAD",
    "EL-NODAL STRESS" : "STRS",
    "EL-ELASTIC STRAINS" : "EPEL",
    "EL-PLASTIC STRAINS" : "EPPL",
    "EL-CREEP STRAINS" : "EPCR",
    "EL-NODAL GRADIENTS" : "FGRAD",
    "EL-INTEGRATION POINT LOCATIONS" : "LOCI",
    "EL-ENERGIES" : "VENG",
    "EL-MISCELLANEOUS" : "MISC",
    "EL-NONLINEAR DATA" : "NLDAT",
    "EL-SURFACE STRESSES" : "SRFS",
    "EL-CONTACT DATA" : "CONT",
    "EL-BACKSTRESSES" : "BKSTR",
    "EL-EULER ANGLES" : "EANGL",
    "ALL NODAL-AVG-SOL" : "NAR",
    "NODAL-AVG STRESSES" : "NDST",
    "NODAL-AVG ELASTIC STRAINS" : "NDEL",
    "NODAL-AVG PLASTIC STRAINS" : "NDPL",
    "NODAL-AVG CREEP STRAINS" : "NDCR"
}

#############################################
#### LOAD STEP 
#############################################


def step_controls(mapdl: core.Mapdl, load_step: LoadStep):
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

    # Load Step Time
    mapdl.time(str(load_step.step_controls.step_end_time))

def solver_controls(mapdl: core.Mapdl, load_step: LoadStep):
    # Solver Type
    if load_step.solver_controls.solver_type == "DIRECT":
        mapdl.eqslv(lab="SPARSE")
    elif load_step.solver_controls.solver_type == "ITERATIVE":
        mapdl.eqslv(lab="PCG")
    else: return # Raise "InvalidParameter" exception

    # Weak Springs
    # Inertia Relif
    # Quasi-Static Analysis

    # Pivot Check
    mapdl.pivcheck(key=load_step.solver_controls.pivot_check,
                   prntcntrl=load_step.solver_controls.pivot_check_print)
    
    # Specify non-linear geometry option
    if load_step.step_number == 1:  # The NLGEOM command cannot be changed after the first load step.
        mapdl.nlgeom(load_step.nonlinear_geo)

def restart_controls(mapdl: core.Mapdl, load_step: LoadStep):
    
    mapdl.rescontrol(action=load_step.restart_controls.action,
                    ldstep=str(load_step.restart_controls.load_step),
                    frequency=str(load_step.restart_controls.frequency))

def nonlinear_controls(mapdl: core.Mapdl, load_step: LoadStep):
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

def output_controls(mapdl: core.Mapdl, load_step: LoadStep):
    # OUTRES for post-processing
    for out in load_step.outputs:
        mapdl.outres(item=output_mapping[out[0]], freq=out[1], cname=out[2])

def nonlinear_diagnostics(mapdl: core.Mapdl, load_step: LoadStep):
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
    

def _load_step(mapdl: core.Mapdl, *args):

    load_step : LoadStep = args[0]

    mapdl.allsel() # !!!!
    mapdl.antype(antype=load_step.analysis_type, status="NEW")
    step_controls(mapdl, load_step)
    solver_controls(mapdl, load_step)
    restart_controls(mapdl, load_step)
    nonlinear_controls(mapdl, load_step)
    output_controls(mapdl, load_step)
    nonlinear_diagnostics(mapdl, load_step)

    # Cutback Controls
    # mapdl.time # Time at the end of a load step

    ########################
    ### Load - BC
    ########################

    # Reference Temperature
    mapdl.tref(tref=str(load_step.reference_temp))

    # Loading Style
    if load_step.loading_style == "RAMPED":
        mapdl.kbc(key=str(0))
    elif load_step.loading_style == "STEPPED":
        mapdl.kbc(key=str(1))

    for force in load_step.forces:
        # Adjust direction syntax from stbox to MAPDL
        _direction = f"F{force.direction}"

        if force.operation == "NEW":
            mapdl.fcum(oper="REPL") # stbox = new (replace previous force value)
            set_force(mapdl, force)
            mapdl.fcum() # reset force accumulation
        elif force.operation == "ADD":
            mapdl.fcum(oper="ADD") # stbox = add (add to previous force value)
            set_force(mapdl, force)
            mapdl.fcum() # reset force accumulation
        elif force.operation == "DELETE":
            if force.fixed:
                mapdl.fdele(node=str(force.set), lab=_direction, lkey="FIXED")
            else:
                mapdl.fdele(node=str(force.set), lab=_direction)

    for moment in load_step.moments:

        if moment[3] == "NEW":
            mapdl.fcum(oper="REPL")
            mapdl.f(node=moment[0], lab=moment[1].replace("R", "M"), value=moment[2])
        elif moment[3] == "ADD":
            mapdl.fcum(oper="ADD")    
            mapdl.f(node=moment[0], lab=moment[1].replace("R", "M"), value=moment[2])
        elif moment[3] == "DELETE":
            if moment[4]:
                mapdl.fdele(node=force[0], lab=f"F{force[1]}", lkey="FIXED")
            else:
                mapdl.fdele(node=force[0], lab=f"F{force[1]}")

    for pressure in load_step.pressures: 
        if pressure.operation == "NEW":
            mapdl.sfcum(lab="PRES", oper="REPL")
            set_pressure(mapdl, pressure)
            #mapdl.sfcum() # Reset
        elif pressure.operation == "ADD":
            mapdl.sfcum(lab="PRES", oper="ADD")
            set_pressure(mapdl, pressure)
            #mapdl.sfcum() # Reset
        elif pressure.operation == "DELETE":
            mapdl.sfdele(nlist=str(pressure.set.nset), lab="PRES")
            
    for displacement in load_step.displacements:
        # Adjust direction syntax from stbox to MAPDL
        _direction : str
        if displacement.direction == "ALL": _direction = displacement.direction
        elif displacement.direction in ("X", "Y", "Z"): _direction = f"U{displacement.direction}"
        elif displacement.direction in ("RX", "RY", "RZ"): _direction = f"ROT{displacement.direction[1]}"

        if displacement.operation == "NEW":
            mapdl.dcum(oper="REPL")
            set_displacement(mapdl, displacement)
        elif displacement.operation == "ADD":
            mapdl.dcum(oper="ADD")
            set_displacement(mapdl, displacement)
        elif displacement.operation == "DELETE":
            if displacement.ramping:
                mapdl.ddele(node=str(displacement.set.items), lab=_direction, rkey="ON")
            else:
                mapdl.ddele(node=str(displacement.set.items), lab=_direction, rkey="OFF")

def set_force(mapdl: core.Mapdl, force : Force):
    # Adjust direction syntax from stbox to MAPDL
    _direction = f"F{force.direction}"
    if force.applied_by == "NODE SET":
        mapdl.csys(force.csys.midx)
        mapdl.f(node=force.set.name, lab=_direction, value=str(force.value))
    elif force.applied_by == "ELEMENT FACES":
        area = get_surface_area(mapdl, force.set)
        pressure = force.value / area # Convert Force to Surface Traction (Pressure)
        mapdl.nsel(type_='S', vmin=force.set.nset.items) # type: ignore
        mapdl.esln(type_='S', ekey="0", nodetype="ALL") 
        _sfcontrol(mapdl=mapdl, direction=_direction, coordinate_system=force.csys, loaded_area="INITIAL")
        mapdl.sfe(elem="ALL", lab="PRES", val1=str(pressure))
        _sfcontrol(mapdl=mapdl, reset=True) # Reset SFCONTROL after SFE command
    elif force.applied_by == "SURFACE ELEMENTS":
        area = get_surface_area(mapdl, force.set)
        pressure = force.value / area # Convert Force to Surface Traction (Pressure)
        mapdl.nsel(type_='S', vmin=force.set.nset.items) # type: ignore
        mapdl.esel(type_='S', item="type", vmin=str(force.set.etype.midx)) # type: ignore
        mapdl.sfe(elem="ALL", lab="PRES", val1=str(force.value))

    #_reset_attributes(mapdl)
    mapdl.allsel() # Reset selection after force application

def _sfcontrol(mapdl: core.Mapdl, direction=None, coordinate_system=None, loaded_area=None, reset=False):
    """
    Sets the SFCONTROL parameters for directional pressure/traction application.
    """
    kcsys = lcomp = val1 = val2 = val3 = ktaper = kuse = karea = kproj = kfollow = ''

    if reset:
        with mapdl.non_interactive:
            mapdl.run(f"SFCONTROL, none") # Reset SFCONTROL
        return

    if direction == "NORMAL TO":
        kcsys = "0"
        lcomp = "0"
        val1 = "0" # use global CSYS
    else: # Apply Pressure with component directions
        kcsys = "1"
        if direction == "X":
            lcomp = "0" # X
        elif direction == "Y":
            lcomp = "1" # Y     
        elif direction == "Z":
            lcomp = "2" # Z
        #mapdl.csys(coordinate_system_id)
        val1 = str(coordinate_system.midx) # type: ignore

    if loaded_area == "DEFORMED":
        karea = "0"
    else:
        karea = "1" 

    with mapdl.non_interactive:
        mapdl.run(f"SFCONTROL,{kcsys},{lcomp},{val1},{val2},{val3},{ktaper},{kuse},{karea},{kproj},{kfollow}")

def set_pressure(mapdl: core.Mapdl, pressure : Pressure):

    # Pressure Application Methods ::
    # 1) Pressure with Node Selection :: SF Command
    # 2) Pressure with Surface Element Selection :: SFE Command after ESEL
    # 3) Pressure with Element Faces :: SFE Command
    if pressure.applied_by == "NODE SET":
        # Select surface nodes by nset
        mapdl.nsel(type_='S', vmin=pressure.set.nset.items) # type: ignore
        mapdl.sf(nlist="ALL", lab="PRES", value=str(pressure.value))
    elif pressure.applied_by == "ELEMENT FACES":
        # Select surface nodes by nset
        mapdl.nsel(type_='S', vmin=pressure.set.nset.items) # type: ignore
        # Select element only if all of its nodes are in the selected nodal set
        mapdl.esln(type_='S', ekey="0", nodetype="ALL") 
        _sfcontrol(mapdl=mapdl, direction=pressure.direction, coordinate_system=pressure.csys, loaded_area=pressure.loaded_area)
        mapdl.sfe(elem="ALL", lab="PRES", val1=str(pressure.value))
        _sfcontrol(mapdl=mapdl, reset=True) # Reset SFCONTROL after SFE command
    elif pressure.applied_by == "SURFACE ELEMENTS":
        # Selecet elements by type id
        mapdl.nsel(type_='S', vmin=pressure.set.nset.items) # type: ignore
        mapdl.esel(type_='S', item="type", vmin=str(pressure.set.element_type.midx)) # type: ignore
        #mapdl.csys(pressure.coordinate_system_id)
        mapdl.sfe(elem="ALL", lab="PRES", val1=str(pressure.value))

    mapdl.csys(0) # Reset to global CSYS
    mapdl.allsel() # Reset selection after surface load application

def set_displacement(mapdl: core.Mapdl, displacement : Displacement):
    # Adjust direction syntax from stbox to MAPDL
    _direction : str
    if displacement.direction == "ALL": _direction = displacement.direction
    elif displacement.direction in ("X", "Y", "Z"): _direction = f"U{displacement.direction}"
    elif displacement.direction in ("RX", "RY", "RZ"): _direction = f"ROT{displacement.direction[1]}"
    mapdl.d(node=displacement.set.name, lab=_direction, value=str(displacement.value))
    mapdl.allsel()

#############################################
#### SOLUTION 
#############################################

def _restart(mapdl: core.Mapdl, *args):
    mapdl.antype(status="RESTART", ldstep=args[0], substep=args[1], action="RSTCREATE")
    pass

def _change_job_name(mapdl: core.Mapdl, job_name : str):
    mapdl.finish()
    mapdl.filname(job_name)
    
def _change_directory(mapdl: core.Mapdl, path):
    mapdl.cwd(path.as_posix())

def _enter_solution(mapdl: core.Mapdl, *args):
    mapdl.slashsolu()

def _solve(mapdl: core.Mapdl, *args):  
    mapdl.allsel()
    mapdl.solve()
    mapdl.save()

def _write_input_file(mapdl: core.Mapdl, *args):
    mapdl.cdwrite(option="DB", fname="model.cdb")

def _exit_mapdl(mapdl: core.Mapdl, *args):
    #mapdl.cdwrite("DB", "COMMANDS.cdb")
    #mapdl.lgwrite("COMMANDS.log")
    mapdl.finish()
    mapdl.exit()

def _write_comment(mapdl: core.Mapdl, *args):
    mapdl.com(comment=args[0]) # type: ignore


            









