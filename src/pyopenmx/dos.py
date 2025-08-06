# Script to plot DOS and PDOS from collinear calculations
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib import style

style_file = Path(__file__).parent / "defectpl.mplstyle"
style.use(style_file)

# Make all the fonts arial
plt.rcParams["font.family"] = "Arial"
plt.rcParams["text.usetex"] = False
colors = ["#2ca02c", "#1f77b4", "#ff7f0e"]


class PlotDOS:
    def __init__(self, dos_file, root_path="./", spinpolarization="on"):
        self.root_path = Path(root_path)
        self.dos_file = self.root_path / dos_file
        self.dos_data = self.load_dos(spinpolarization=spinpolarization)

    def load_dos(self, spinpolarization="on"):
        """Load DOS data from the specified file."""
        data = np.loadtxt(self.dos_file)
        energies = data[:, 0]  # First column is energy
        if spinpolarization == "on":
            dos_up = data[:, 1]
            dos_down = data[:, 2]
        else:
            dos_up = data[:, 1]
            dos_down = np.zeros_like(dos_up)
        data = {}
        data["energies"] = energies
        data["dos_up"] = dos_up
        data["dos_down"] = dos_down
        data["spinpolarization"] = spinpolarization
        return data

    def plot_dos(
        self, figsize=(6, 4), save_file="dos_plot.svg", energy_range=None, title=None
    ):
        """Plot the DOS and PDOS."""
        plt.figure(figsize=figsize)
        energies = self.dos_data["energies"]
        dos_up = self.dos_data["dos_up"]
        plt.plot(energies, dos_up, label="Up", color=colors[0])
        if self.dos_data["spinpolarization"] == "on":
            dos_down = self.dos_data["dos_down"]
            plt.plot(energies, dos_down, label="Down", color=colors[1])

        plt.xlabel("Energy (eV)")
        plt.ylabel("DOS (states/eV)")
        plt.title(title if title else "Density of States")
        if energy_range:
            plt.xlim(energy_range)
        else:
            plt.xlim(energies.min(), energies.max())
        plt.legend()
        plt.grid(alpha=0.5, linestyle="--")
        plt.tight_layout()
        if save_file:
            plt.savefig(save_file)
        else:
            plt.show()


if __name__ == "__main__":
    # Example usage
    dos_file = "f3gt.DOS.Tetrahedron"  # Replace with your DOS file

    plotter = PlotDOS(dos_file, spinpolarization="on")
    plotter.plot_dos(
        save_file="f3gt_dos_tetrahedron.svg",
        title="Fe$_3$GeTe$_2$",
        energy_range=(-5, 5),
        figsize=(4, 4),
    )  # Save the plot as SVG
