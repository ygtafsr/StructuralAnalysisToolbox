
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
    model.import_mesh(r'C:\Users\yigit\StructuralAnalysisToolbox\examples\meshes\cantilevier_beam.cdb')
    model.add_material("Linear Steel", "MESH")

    model.add_pilot_node("Force_Node", -15, 50 ,-300)

    load_step_1 = model.add_loadstep(name="LS-1")
    load_step_1.remote_force("Force_Node", "LOADING_NODES_2", -5000, "Y")
    load_step_1.dof("NS_FIX", 0, "ALL")
    load_step_1.output("ALL")

    load_step_2 = model.add_loadstep(name="LS-2")
    load_step_2.remote_force("Force_Node", "LOADING_NODES_2", 3000, "Y", operation="ADD")
    load_step_1.dof("NS_FIX", 0, "ALL")
    load_step_2.output("ALL")

    load_step_3 = model.add_loadstep(name="LS-3")
    load_step_2.remote_force("Force_Node", "LOADING_NODES_2", 3000, "Y", operation="NEW")
    load_step_1.dof("NS_FIX", 0, "ALL")
    load_step_3.output("ALL")

    model.plot_load_history("Force_Node", "Force", "Y")

    model.info()

    model.write_input_file()
    #model.solve()


if __name__ == "__main__":
    main()