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

    def __init__(self, 
                 load_list : np.ndarray, 
                 nodes: str, 
                 load_type : Literal["Force", "Moment", "Pressure", "Displacement"],
                 direction: str):
        
        self.fig, self.ax = plt.subplots()
        
        # Axis
        self.ax.set_title(f"{nodes} load history - {direction}")
        self.ax.set_xlabel("Load Step")
        self.ax.set_ylabel(f"{load_type} (mm)")

        xticks = np.arange(0, len(load_list), 1)
        self.ax.plot(xticks, load_list)
        self.ax.set_xticks(ticks=xticks)
        self.ax.set_yticks(ticks=load_list)
        self.ax.grid(True)

    def plot(self):
        self.fig.tight_layout()
        self.ax.plot()