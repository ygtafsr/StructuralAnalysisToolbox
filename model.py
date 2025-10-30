
from ansys.mapdl.core import launch_mapdl

class start:

    def __init__(self, solver='MAPDL', work_folder = ''):
        
        self.solver = solver
        self.work_folder = work_folder
    

    def start(self, cpus = 4):

        mapdl = launch_mapdl(run_location=self.work_folder,
                             nproc=cpus,
                             override=True)
        # First, reset the MAPDL database.
        mapdl.clear()

        return mapdl
    