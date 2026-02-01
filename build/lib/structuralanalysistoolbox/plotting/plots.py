from __future__ import annotations
from typing import Literal
import matplotlib
#matplotlib.use(backend="inline")
import matplotlib.pyplot as plt
import numpy as np

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from structuralanalysistoolbox.model import Analysis

class LoadStepPlot:

    def __init__(self, solution : Analysis, nodes: str, 
                 load_type : Literal["Force", "Moment", "Displacement"],
                 direction: str):
        
        self.fig, self.ax = plt.subplots()
        
        # Axis
        self.ax.set_title(f"{nodes} load history - {direction}")
        self.ax.set_xlabel("Load Step")

        # Data
        if load_type == "Force":
            self.ax.set_ylabel(f"{load_type} (N)")
            history = solution._get_force_history(nodes=nodes, direction=direction)
        elif load_type == "Moment":
            self.ax.set_ylabel(f"{load_type} (Nmm)")
            history = solution._get_moment_history(nodes=nodes, direction=direction)
        elif load_type == "Displacement":
            self.ax.set_ylabel(f"{load_type} (mm)")

        xticks = np.arange(0, len(history), 1)
        self.ax.plot(xticks, history)
        self.ax.set_xticks(ticks=xticks)
        self.ax.grid(True)

    def plot(self):
        self.fig.tight_layout()
        self.ax.plot()