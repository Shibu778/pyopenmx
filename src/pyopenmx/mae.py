# Script to plot Magnetocrystalline Anisotropy Energy (MAE) from SVM GGA PBE calculations
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

Ha2eV = 27.2114  # Conversion factor from Hartree to eV
Ha2meV = Ha2eV * 1000  # Conversion factor from Hartree to meV


class PlotMAE:
    def __init__(self, filelist, labellist, root_path="./", number_of_atoms=3):
        self.root_path = Path(root_path)
        self.filelist = filelist
        self.labellist = labellist
        self.data = self.extract_all_data()
        self.minima = self.minimum_energy()
        self.number_of_atoms = number_of_atoms

    def extract_data(self, filename):
        """Extracts MAE data from a given file."""
        data = np.loadtxt(filename, skiprows=1)
        angles = data[:, 0]  # First column is angles in degrees
        values = data[:, 1]  # Second column is Uele in Hartree
        return data

    def extract_all_data(self):
        """Extracts data from all files and stores it in a dictionary."""
        data = {}
        for filename, label in zip(self.filelist, self.labellist):
            filepath = self.root_path / filename
            if filepath.exists():
                data[label] = self.extract_data(filepath)
            else:
                print(f"File {filepath} does not exist.")
        return data

    def plot_data(
        self,
        unit_conversion=Ha2meV,
        figsize=(4, 4),
        save_file="mae_plot_gga_pbe.svg",
        colors=colors,
        title=f"Fe$_3$GeTe$_2$",
        grid=False,
        legend_loc="lower center",
    ):
        """Plots the extracted MAE data."""
        plt.figure(figsize=figsize)
        for label, values in self.data.items():
            angles = values[:, 0]
            uele = values[:, 1] * unit_conversion
            min_uele = min(self.minima) * unit_conversion
            uele -= min_uele
            uele /= self.number_of_atoms  # Normalize by number of atoms
            plt.plot(angles, uele, marker="o", label=label, color=colors.pop(0))

        if unit_conversion == Ha2meV:
            plt.ylabel("MAE (meV/Fe atom)")
        else:
            plt.ylabel("MAE (Hartree/Fe atom)")
        plt.xlabel("$\\theta$ (deg)")
        plt.title(title)
        plt.legend(loc=legend_loc)
        if grid:
            plt.grid(alpha=0.5, linestyle="--")
        plt.xticks(np.arange(0, 181, 30))
        plt.tight_layout()
        if save_file:
            plt.savefig(save_file)
        else:
            plt.show()

    def minimum_energy(self):
        """Finds the minimum MAE value across all datasets."""
        min_values = []
        for label, values in self.data.items():
            uele = values[:, 1]
            min_values.append(np.min(uele))
        return min_values


if __name__ == "__main__":
    # Define the file paths and labels
    filelist = [
        "xy_plane_angle_uele.info",
        "xz_plane_angle_uele.info",
        "yz_plane_angle_uele.info",
    ]
    labellist = ["xy", "xz", "yz"]

    # Create an instance of PlotMAE
    plotter = PlotMAE(filelist, labellist)
    plotter.plot_data(title=f"Fe$_3$GeTe$_2$", save_file="f3gt_mae_plot_gga_pbe.svg")
