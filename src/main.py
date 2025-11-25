
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
        
    model.add_MPC_Rigid(dependent="NS_LOAD_DEPEN", independent="NS_LOAD")
    model.info()
    
    model._execute_commands()

if __name__ == "__main__":
    main()