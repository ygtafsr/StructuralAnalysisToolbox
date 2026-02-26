from dataclasses import dataclass
from typing import Literal


'''
For attributes Real, Mat and Section, an empty attribute 
can be created for an attribute id placeholder. (Ansys allow this)
Therefore, during mapdl command implementation, these ids can be used
to element creation for dedicated attributes that are in their default 
(unassigned) forms.
'''

@dataclass
class Attribute:
    midx : int
    type : Literal["Mat", "Type", "Real", "Csys", "Secn"]
    
@dataclass
class Etype(Attribute):
    midx : int = 0  # Element type index (Same as in solver)
    type : Literal["Mat", "Type", "Real", "Csys", "Secn"] = "Type"
    name : str = ''
    ignore_at_execution : bool = False
    show_in_model_tree = True
    
@dataclass
class Real(Attribute):
    midx : int = 0
    type : Literal["Mat", "Type", "Real", "Csys", "Secn"] = "Real"
    ignore_at_execution : bool = False

@dataclass
class Mat(Attribute):
    midx : int = 0
    type : Literal["Mat", "Type", "Real", "Csys", "Secn"] = "Mat"
    ignore_at_execution : bool = False

@dataclass
class Section(Attribute):
    """
    SECTYPE, SECID, Type, Subtype, Name, REFINEKEY
    Associates section type information with a section ID number.
    https://ansyshelp.ansys.com/public/account/secured?returnurl=///Views/Secured/corp/v242/en/ans_cmd/Hlp_C_SECTYPE.html
    """
    midx : int = 0
    type : Literal["Mat", "Type", "Real", "Csys", "Secn"] = "Secn"
    ignore_at_execution : bool = False


