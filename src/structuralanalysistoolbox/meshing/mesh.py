
from __future__ import annotations
from dataclasses import dataclass
from enum import Enum

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from structuralanalysistoolbox import model

class Node:
    def __init__(self, model : model.Model):

        self.id : int
        self.dof : int
        self.mapdl = model.mapdl

    def _create_free_node(self):
        """finds max node id, creates a new one with the next number."""

class Etype(Enum):
    
    # SURFACE ELEMENTS
    CONTA172 = "CONTA172"
    CONTA174 = "CONTA174"
    CONTA175 = "CONTA175"
    CONTA177 = "CONTA177"
    TARGE169 = "TARGE169"
    TARGE170 = "TARGE170"
    # CONSTRAINTS
    MPC184 = "MPC184"
  
@dataclass
class ElementType():
    """
    Defines a local element type from the element library.
    """
    id: int
    type: Etype
    keyopt: list[tuple[int, int]]
    inopr: int = 0

    def get_keyop(self, key:int):
        keyop_val = next((t[1] for t in self.keyopt if t[0] == key), None)
        return keyop_val



class Mesh:
    
    x:str