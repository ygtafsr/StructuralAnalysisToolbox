import ansys.mapdl.core as pymapdl
from ansys.mapdl.core import launch_mapdl
from pathlib import Path
import os


class model:

    def __init__(self, solver='MAPDL', name = '', cpus = 4, interactive = True, work_folder = ''):
        
        self.solver = solver
        self.modelname = name
        self.cpus = cpus
        self.interactive = interactive
        self.work_folder = work_folder

        if work_folder == '':
            _current_path = Path.cwd()
            self.working_directory = _current_path / name
            self.working_directory.mkdir()
        
        if interactive:
            self.start()

    def start(self):

        self.mapdl = launch_mapdl(run_location=self.working_directory,
                                  jobname=self.modelname,
                                  nproc=self.cpus,
                                  override=True)
        
        # First, reset the MAPDL database.
        self.mapdl.clear()

        return self.mapdl
    
    def stop(self):
        pymapdl.close_all_local_instances()
    

    


    