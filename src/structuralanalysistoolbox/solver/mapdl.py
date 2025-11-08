# Material Model Base File
## MP: Defines a linear material property as a constant or a function of temperature.
# Material Property Label (lab)
## MPDATA: Defines property data associated with the temperature table
## MPTEMP: Defines a temperature table for material properties.
## MPLIST: Lists linear material properties.
# *****************************************************************************************
## TB (Create Material Data Table): Activates a data table for material properties or special element input.
## TBDATA: Defines data for the material data table.
## TBFIELD
## TBTEMP
# *****************************************************************************************
## MPWRITE: Writes linear material properties in the database to a file (if the LIB option is not specified)
#  or writes both linear and nonlinear material properties (if LIB is specified) from the database to a file.
## MPREAD: Reads a file containing material properties.
## MPLIB: Sets the default material library read and write paths.

import ansys.mapdl.core as core
import numpy as np
from dataclasses import asdict
#from .. import material
#from .. import model
from materials import material
#import model

# If a material displays nonlinear stress-strain behavior, use the TB family of commands to define the
# nonlinear material property relationships in terms of a data table: TB, TBTEMP, TBDATA, TBPT, TBCOPY,
# TBLIST, TBPLOT, TBDELE


def _material(mapdl : core.Mapdl, mat : material.Material, mat_id : int) -> None:
    """Runs mapdl material definition commands"""

    mapdl.prep7()
    
    for model_category, model_object in mat.material_models.items():
        for label, data in asdict(model_object).items():

            if label == "density":
                mapdl.load_table(name=f"{label}_{mat_id}", array=np.array(data), var1="TEMP")
                mapdl.mp(lab="DENS", mat=f"{mat_id}", c0=f"%{label}_{mat_id}%")            
            elif label == "elastic_modulus":
                mapdl.load_table(name=f"{label}_{mat_id}", array=np.array(data), var1="TEMP")
                mapdl.mp(lab="EX", mat=f"{mat_id}", c0=f"%{label}_{mat_id}%")
            elif label == "poisson_ratio":
                    mapdl.load_table(name=f"{label}_{mat_id}", array=np.array(data), var1="TEMP")
                    mapdl.mp(lab="NUXY", mat=f"{mat_id}", c0=f"%{label}_{mat_id}%")
            elif model_category == "PLASTICITY" and type(model_object).__name__ == "MultilinearIsotropicHardening":
                    data_table = mat.material_models[model_category].data_table # [25, [[a,b], [a,b], [a,b],...]]
                    temp_data_number = len(data_table)
                    plastic_data_number = len(data_table[1])
                    mapdl.tb(lab="PLASTIC", mat=f"{mat_id}", ntemp=temp_data_number, 
                                                             npts=plastic_data_number)
                    for temp_idx in data_table:
                         mapdl.tbtemp(temp=f"{data_table[temp_idx]}")
                         for data in data_table[temp_idx]:
                              mapdl.tbpt(oper="DEFI", x1=f"{data[0]}", x2=f"{data[1]}")

    mapdl.finish()
                            
                    