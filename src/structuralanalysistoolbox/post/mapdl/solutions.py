
from ansys.dpf import post
from pathlib import Path

class DpfSolution:
    def __init__(self, file_path: Path):
        """Initialize DpfSolution object."""

        self._simulation = post.load_simulation(file_path)     

    @property
    def simulation(self):
        """Get the DPF simulation instance."""
        return self._simulation

