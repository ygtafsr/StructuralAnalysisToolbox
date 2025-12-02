
import structuralanalysistoolbox as stbox
from pathlib import Path
from ansys.mapdl import reader as pymapdl_reader
from ansys.mapdl.reader import examples


def _add_material():

    physical = stbox.material.Physical(density=[[25, 8000], [50, 9000]])
    elastic = stbox.material.IsotropicElasticity(elastic_modulus=210000, poisson_ratio=0.3)
    multilinear_csv = stbox.material.import_csv_data(path=r'src\structuralanalysistoolbox\materials\plastic_data.csv')
    plastic = stbox.material.MultilinearIsotropicHardening(data_table=multilinear_csv)

    mat_1 = stbox.material.Material(name='My Steel-5',
                                    category='Steel', 
                                    remarks='It is a good steel',
                                    source='My Mind',
                                    last_update='Today')

    mat_1.add_model(physical)
    mat_1.add_model(elastic)
    mat_1.add_model(plastic)

    stbox.material.library.save(mat_1)

    
def main():

    model = stbox.Model(name="my-model", interactive=False)
    model.import_mesh(r'src\rod.cdb')

    mat = model.add_material(mat="My Steel")
    model.assign_material("MESH", mat)

    model.add_MPC_Rigid(dependent="NS_LOAD_DEPEN", independent="NS_LOAD")
    model.info()

    load_step_1 = model.add_load_step()

    load_step_1.force("NS_LOAD", "X", 8000)
    load_step_1.dof("NS_FIX", "ALL", 0)

    load_step_1.output("NODAL DOF", "ALL")
    load_step_1.output("REACTION LOADS", "ALL")
    load_step_1.output("NODAL-AVG PLASTIC STRAINS", "ALL")
    load_step_1.output("NODAL-AVG STRESSES", "ALL")

    load_step_2 = model.add_load_step()
    load_step_2.force("NS_LOAD", "X", -2000, operation="ADD")

    load_step_3 = model.add_load_step()
    load_step_3.force("NS_LOAD", "Y", 5000)

    load_step_4 = model.add_load_step()
    load_step_4.force("NS_LOAD", "X", 5000)
  
    

    model.plot_load("NS_LOAD", "Force", "X")
    #model.solve()
    pass

if __name__ == "__main__":
    main()