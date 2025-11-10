
import ansys.mapdl.core as core
import numpy as np
from dataclasses import asdict
from ..materials import material
"""from materials import material"""


def _material(mapdl : core.Mapdl, mat : material.Material, mat_id : int) -> None:
    """Runs mapdl material definition commands"""

    mapdl.prep7()
    
    for model in mat.material_models.values():

        if isinstance(model, material.Physical):
            # Density
            if isinstance(model.density, (float, int)):
                mapdl.mp(lab="DENS", mat=f"{mat_id}", c0=f"{model.density}")
            else:  
                label = "density"
                mapdl.load_table(name=f"{label}_{mat_id}", array=np.array(model.density), var1="TEMP")
                mapdl.mp(lab="DENS", mat=f"{mat_id}", c0=f"%{label}_{mat_id}%")
        elif isinstance(model, material.IsotropicElasticity):
            # Elastic Modulus
            if isinstance(model.elastic_modulus, (float, int)):
                mapdl.mp(lab="EX", mat=f"{mat_id}", c0=f"{model.elastic_modulus}")
            else:
                label = "elastic_modulus"
                mapdl.load_table(name=f"{label}_{mat_id}", array=np.array(model.elastic_modulus), var1="TEMP")
                mapdl.mp(lab="EX", mat=f"{mat_id}", c0=f"%{label}_{mat_id}%")
            # Poisson's Ratio
            if isinstance(model.poisson_ratio, (float, int)):
                mapdl.mp(lab="NUXY", mat=f"{mat_id}", c0=f"{model.poisson_ratio}")
            else:
                label = "poissons_ratio"
                mapdl.load_table(name=f"{label}_{mat_id}", array=np.array(model.poisson_ratio), var1="TEMP")
                mapdl.mp(lab="NUXY", mat=f"{mat_id}", c0=f"%{label}_{mat_id}%")
        elif isinstance(model, material.MultilinearIsotropicHardening):
            # Multilinear Isotropic Hardening
            if model.temperature == None:
                data_count = len(model.data_table)
                mapdl.tb(lab="PLASTIC", mat=f"{mat_id}", ntemp="1", 
                                                            npts=f"{data_count}")
                mapdl.tbtemp(temp="")
                for data in model.data_table:
                    mapdl.tbpt(oper="DEFI", x1=f"{data[0]}", x2=f"{data[1]}")
            elif isinstance(model.temperature, (float, int)):
                data_count = len(model.data_table)
                mapdl.tb(lab="PLASTIC", mat=f"{mat_id}", ntemp="1", 
                                                            npts=f"{data_count}")
                mapdl.tbtemp(temp=f"{model.temperature}")
                for data in model.data_table:
                    mapdl.tbpt(oper="DEFI", x1=f"{data[0]}", x2=f"{data[1]}")

    mapdl.finish()
                            
                    