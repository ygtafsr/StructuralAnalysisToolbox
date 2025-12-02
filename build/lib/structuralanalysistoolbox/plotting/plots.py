from __future__ import annotations
from typing import Literal
import matplotlib.pyplot as plt
import numpy as np

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from structuralanalysistoolbox.model import Solution

class LoadStepPlot:

    def __init__(self, solution : Solution, nodes: str, 
                 load_type : Literal["Force", "Moment", "Displacement"],
                 direction: str):
        
        self.fig, self.ax = plt.subplots()
        
        self.ax.set_title(f"{nodes} load history - {direction}")

        self.ax.set_xlabel("Load Step")
        

        if load_type == "Force":
            self.ax.set_ylabel(f"{load_type} (N)")
            history = solution.get_force_history(nodes=nodes, direction=direction)
            xticks = np.arange(0, len(history), 1)
            self.ax.set_xticks(ticks=xticks)
            self.ax.plot(xticks, history)
            self.ax.grid(True)
            
        elif load_type == "Moment":
            self.ax.set_ylabel(f"{load_type} (Nmm)")
            self.ax.plot(solution.get_moment_history(nodes=nodes, direction=direction))
        elif load_type == "Displacement":
            self.ax.set_ylabel(f"{load_type} (mm)")

    def show(self):
        #self.fig.tight_layout()
        #self.fig.show()
        return self.ax