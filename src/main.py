
import structuralanalysistoolbox as stbox

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

    model = stbox.Model(name="my-model")
    model.add_mesh(mesh_file_path=r'C:\Users\yigit\StructuralAnalysisToolbox\src\mesh_rod.cdb')

    # NSET for constraint node
    depen = model._model_definitions["sets"]["NS_LOAD_DEPEN"][0]
    indepen = model._model_definitions["sets"]["NS_LOAD"][0]

    model.add_MPC(name="MPC-1", type="RigidLink", dependent=depen, independents=indepen)
    print(model)    

if __name__ == "__main__":
    main()