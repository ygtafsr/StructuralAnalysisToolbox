
# STBOX ELEMENT TYPES

from dataclasses import dataclass
from typing import Literal  


@dataclass
class HEX8:

    #Constant Attributes
    name = "HEX8"
    pyvista_type = "HEXAHEDRON"
    dimension = 3
    points = 8
    order = "Linear"
   
@dataclass
class WEDGE6:

    #Constant Attributes
    name = "WEDGE6"
    pyvista_type = "WEDGE"
    dimension = 3
    points = 6
    order = "Linear"

@dataclass
class PYRAMID5:

    #Constant Attributes
    name = "PYRAMID5"
    pyvista_type = "PYRAMID"
    dimension = 3
    points = 5
    order = "Linear"

@dataclass
class TETRA4:

    #Constant Attributes
    name = "TETRA4"
    pyvista_type = "TETRA"
    dimension = 3
    points = 4
    order = "Linear"
    
@dataclass
class TRIA3:

    #Constant Attributes
    name = "TRIA3"
    pyvista_type = "TRIANGLE"
    dimension = 2
    points = 3
    order = "Linear"
    
@dataclass
class QUAD4:

    #Constant Attributes
    name = "QUAD4"
    pyvista_type = "QUAD"
    dimension = 2
    points = 4
    order = "Linear"
    
ELEMENT_TYPE_MAPPING = {
    "HEXAHEDRON" : HEX8,
    "WEDGE"     : WEDGE6,
    "TETRA"     : TETRA4,
    "PYRAMID"   : PYRAMID5,
    "TRIANGLE"  : TRIA3,
    "QUAD"      : QUAD4
}