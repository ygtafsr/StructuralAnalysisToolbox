
from dataclasses import dataclass
from enum import Enum

##########################
## ELEMENT TECHNOLOGIES  
##########################

class Tech185(Enum): 
    full_integration = 0
    reduced_integration = 1
    enhanced_strain = 2
    simplified_enhanced_strain = 3

class Tech186(Enum): 
    reduced_integration = 0
    full_integration = 1

class Tech182(Enum): 
    full_integration = 0
    reduced_integration = 1
    enhanced_strain = 2
    simplified_enhanced_strain = 3

##########################
## ELEMENT FORMULATIONS
##########################

class Form185(Enum): 
    pure_displacement = 0
    mixed_up = 1

class Form186(Enum): 
    pure_displacement = 0
    mixed_up = 1

class Form187(Enum): 
    pure_displacement = 0
    mixed_up_cons_press = 1
    mixed_up_lin_press = 2

class Form182(Enum): 
    pure_displacement = 0
    mixed_up = 1

class Form183(Enum): 
    pure_displacement = 0
    mixed_up = 1

class Form285(Enum):
    mixed_up = 0
    pure_displacement = 1

###############################

class Behv(Enum):
    plane_stress = 0
    axisymmetric = 1
    plane_strain = 2
    plane_stress_with_thickness = 3
    generalized_plane_strain = 5
    axisymmetric_with_torsion = 6
    
class Layer(Enum): 
    non_layered = 0
    layered = 1

class Shape(Enum):
    quad_8 = 0
    tria_6 = 1

class PML(Enum) : 
    not_include = 0
    include = 1

class Steady(Enum):
    disabled = 0
    enabled = 1

class SurfOut(Enum):
    basic = 0
    extra = 4

@dataclass
class Solid185:
    """8-Node Structural Solid"""
    technology : Tech185 = Tech185.full_integration
    formulation : Form185 = Form185.pure_displacement
    layering : Layer = Layer.non_layered
    pml_absorbing : PML =  PML.not_include
    steady_state : Steady = Steady.disabled
    extra_surface_output : SurfOut = SurfOut.basic

@dataclass
class Solid186:
    """20-Node Structural Solid"""
    technology : Tech186 = Tech186.reduced_integration
    formulation : Form186 = Form186.pure_displacement
    layering : Layer = Layer.non_layered
    pml_absorbing : PML = PML.not_include
    steady_state : Steady = Steady.disabled
    extra_surface_output : SurfOut = SurfOut.basic

@dataclass
class Solid187:

    formulation : Form187 = Form187.pure_displacement
    pml_absorbing : PML = PML.not_include
    steady_state : Steady = Steady.disabled
    extra_surface_output : SurfOut = SurfOut.basic

@dataclass
class Solid285:

    formulation : Form285 = Form285.mixed_up
    pml_absorbing : PML = PML.not_include
    steady_state : Steady = Steady.disabled
    extra_surface_output : SurfOut = SurfOut.basic

@dataclass
class Plane182:

    technology : Tech182
    behaviour : Behv
    formulation : Form182 = Form182.pure_displacement
    pml_absorbing : PML = PML.not_include
    extra_surface_output : SurfOut = SurfOut.basic

@dataclass
class Plane183:

    shape : Shape
    behaviour : Behv
    formulation : Form183 = Form183.pure_displacement
    pml_absorbing : PML = PML.not_include
    extra_surface_output : SurfOut = SurfOut.basic