
from enum import Enum
from dataclasses import dataclass, asdict, field
import numpy as np
from . import matlib

class MaterialModelCategory(Enum):
    PHYSICAL = 1
    ELASTICITY = 2
    PLASTICITY = 3
    HYPER_ELASTIC = 4
    STRENGTH = 5
    FATIGUE = 6
    
class DataAlreadyDefined(Exception):

    def __init_(self):
        super().__init__('Data has already defined!')

@dataclass
class Physical:     # --> material model

    model_category = MaterialModelCategory.PHYSICAL
    density : float | list = None # -> material label

    """def __post_init__(self):
        if isinstance(self.density, int):
            self.density = float(self.density)"""

@dataclass
class IsotropicElasticity:  

    model_category = MaterialModelCategory.ELASTICITY
    elastic_modulus : float | list = None  
    poisson_ratio : float | list = None  

@dataclass
class MultilinearIsotropicHardening:

    model_category = MaterialModelCategory.PLASTICITY
    temperature : float | list = 25
    data_table : list = None 

@dataclass
class BilinearIsotropicHardening:

    model_category = MaterialModelCategory.PLASTICITY
    temperature : float = 25
    yield_strength : float = None
    tangent_modulus : float = None

@dataclass
class Strength:

    model_category = MaterialModelCategory.STRENGTH
    tensile_yield_strength : float = None
    tensile_ultimate_strength : float = None

@dataclass
class StrainLife:

    model_category : MaterialModelCategory.FATIGUE
    strength_coefficient : float = None
    strength_exponent : float = None
    ductility_coefficient : float = None
    ductility_exponent : float = None

    cyclic_strength_coefficient : float = None
    cyclic_strength_hardening_exponent : float = None

class Material:

    def __init__(self, name, category = '', remarks = '', source = '', last_update = ''):
        
        # material labels
        self.midx : int = 0
        self.name : str  = name   # MATERIAL NAME MUST BE UNIQUE!
        self.category : str  = category
        self.remarks : str  = remarks
        self.source : str  = source
        self.last_update : str  = last_update

        self.material_models = {} # {"model category" : model object,}

    def add_model(self, data):

        if not data.model_category in [value.model_category for value in self.material_models.values()]:
            self.material_models[data.model_category.name] = data
        else: raise DataAlreadyDefined()

    def save(self):
        if matlib.exist(self.name):
            matlib._update(self)
        else:
            matlib._create(self) 

    def list(self):
        pass

    def plot(self):
        '''Plots Material Data'''
        pass

    def _get(self, mat_label : str) -> list:
        """Return material data for labels (density, elastic_modulus etc.)"""
        for category in self.material_models.values():
            for label, data in category:
                if mat_label == label:
                    return data
        return None

def load(mat_name : str) -> Material:

    if matlib.exist(mat_name):
        new_mat = Material(name=mat_name)
        database_mat = matlib._read(mat_name)
        """Converts database object to a material object"""
        for key, value in database_mat.items():
            if key == "category": new_mat.category = value
            elif key == "remarks": new_mat.remarks = value
            elif key == "source": new_mat.source = value
            elif key == "last update": new_mat.last_update = value
            elif key == "material models":
                for model_category in value.values(): # iterate over database model categories
                    for data_name, data in model_category.items():
                        if data_name == "model name":
                            mat_model_class = globals()[data]
                            mat_model_obj = mat_model_class()
                            new_mat.add_model(mat_model_obj)
                        else:
                            vars(mat_model_obj)[data_name] = data
        return new_mat

def import_csv_data(path : str):
    data_table = np.loadtxt(fname=path, delimiter=',')
    return data_table.tolist()


#######################################
## Calculations for Material Properties
#######################################
def convert_true_to_eng():
    pass

#######################################
## Material Data Plotting
#######################################
def plot(data: str = "True Stress_True Strain"):
    pass

