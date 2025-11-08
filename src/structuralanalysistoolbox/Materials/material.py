
from enum import Enum
from dataclasses import dataclass, asdict
import numpy as np

from .matlib import MatDatabase
library = MatDatabase()

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
    density : list[list] = None # material_data (-> material label)

@dataclass
class IsotropicElasticity:  

    model_category = MaterialModelCategory.ELASTICITY
    elastic_modulus : list[list] = None # [[temp, E],]  
    poisson_ratio : list[list]   = None # [[temp, nu],] 

@dataclass
class MultilinearIsotropicHardening:

    model_category = MaterialModelCategory.PLASTICITY
    data_table : list[list]  = None # [temp, [strain , stress]]

@dataclass
class BilinearIsotropicHardening:

    model_category = MaterialModelCategory.PLASTICITY
    yield_strength : list[list] = None
    tangent_modulus : list[list] = None

@dataclass
class Strength:

    model_category = MaterialModelCategory.STRENGTH
    tensile_yield_strength : float = 0.0
    tensile_ultimate_strength : float = 0.0

@dataclass
class StrainLife:

    model_category : MaterialModelCategory.FATIGUE
    strength_coefficient : float = 0.0
    strength_exponent : float = 0.0
    ductility_coefficient : float = 0.0
    ductility_exponent : float = 0.0

    cyclic_strength_coefficient : float = 0.0
    cyclic_strength_hardening_exponent : float = 0.0

class Material:

    def __init__(self, name, category = '', remarks = '', source = '', last_update = ''):
        
        # material labels
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

    def _generate_commands(self, solver='MAPDL'):
        '''Generate MAPDL material model commands'''
        pass

def load(mat_name : str) -> Material:

    if library._exist(mat_name):
        new_mat = Material(name=mat_name)
        database_mat = library._read(mat_name)
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

