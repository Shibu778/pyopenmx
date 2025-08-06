# Script to plot partial DOS from collinear calculations
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib import style
import sys

style_file = Path(__file__).parent / "defectpl.mplstyle"
style.use(style_file)

# Make all the fonts arial
plt.rcParams["font.family"] = "Arial"
plt.rcParams["text.usetex"] = False
colors = [
    "#2ca02c",
    "#1f77b4",
    "#ff7f0e",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
    "#ffbb78",
    "#98df8a",
]


class PlotPDOS:
    def __init__(
        self,
        pdos_file,
        root_path="./",
        pdos_folder="pdos_tetrahedron",
        spinpolarization="on",
    ):
        self.root_path = Path(root_path)
        self.pdos_folder = self.root_path / pdos_folder
        self.pdos_filename, self.pdos_data = self.get_pdos_files()
        self.spinpolarization = spinpolarization
        self.pdos_data = self.extract_all_pdos()
        self.grouped_orbitals = self.group_orbitals()
        self.grouped_species_orbitals = self.group_species_and_orbitals()

    def get_pdos_files(self):
        """Get all PDOS files in the specified folder."""
        pdos_files = list(self.pdos_folder.glob("*.PDOS.Tetrahedron.*"))
        # Get the filenames without the path
        pdos_files = [file.name for file in pdos_files]
        pdos_data = {}
        for f in pdos_files:
            f = f.split(".")
            if len(f) == 4:
                pdos_data[f[3]] = {}
                pdos_data[f[3]]["atom"] = {}
            elif len(f) == 5:
                pdos_data[f[3]][f[4]] = {}
            else:
                raise ValueError(f"Unexpected file format: {f}")

        return pdos_files, pdos_data

    def load_pdos(self, pdos_file):
        """Load PDOS data from the specified file."""
        data = np.loadtxt(self.pdos_folder / pdos_file)
        energies = data[:, 0]  # First column is energy
        data_dict = {}
        if self.spinpolarization == "on":
            pdos_up = data[:, 1]
            pdos_down = data[:, 2]
        else:
            pdos_up = data[:, 1]
            pdos_down = np.zeros_like(pdos_up)

        data_dict["energies"] = energies
        data_dict["pdos_up"] = pdos_up
        data_dict["pdos_down"] = pdos_down
        data_dict["spinpolarization"] = self.spinpolarization
        return data_dict

    def extract_all_pdos(self):
        """Extract PDOS data for all files."""
        for pdos_file in self.pdos_filename:
            f = pdos_file.split(".")
            if len(f) == 4:
                atom = f[3]
                self.pdos_data[atom]["atom"] = self.load_pdos(pdos_file)
            elif len(f) == 5:
                atom = f[3]
                orbital = f[4]
                self.pdos_data[atom][orbital] = self.load_pdos(pdos_file)
            else:
                raise ValueError(f"Unexpected file format: {f}")

        return self.pdos_data

    def group_orbitals(self):
        """Function to group orbital wise data. For example, d1, d2, d3, d4, d5 will be summed and stored as d."""
        grouped_data = {}
        for atom, data in self.pdos_data.items():
            grouped_data[atom] = {"atom": data["atom"]}
            for orbital, orbital_data in data.items():
                if orbital != "atom":
                    group_key = orbital[0]
                    if group_key not in grouped_data[atom]:
                        grouped_data[atom][group_key] = {
                            "energies": orbital_data["energies"],
                            "pdos_up": 0,
                            "pdos_down": 0,
                        }
                    grouped_data[atom][group_key]["pdos_up"] += orbital_data["pdos_up"]
                    grouped_data[atom][group_key]["pdos_down"] += orbital_data[
                        "pdos_down"
                    ]
        return grouped_data

    def group_species_and_orbitals(
        self, list_of_atoms={"Fe": [1, 2, 3], "Ge": [4], "Te": [5, 6]}
    ):
        """Group PDOS data for a list of atoms."""
        grouped_data = {}
        for atom, indices in list_of_atoms.items():
            grouped_data[atom] = {
                "atom": {"energies": [], "pdos_up": [], "pdos_down": []},
                "s": {"energies": [], "pdos_up": [], "pdos_down": []},
                "p": {"energies": [], "pdos_up": [], "pdos_down": []},
                "d": {"energies": [], "pdos_up": [], "pdos_down": []},
                "f": {"energies": [], "pdos_up": [], "pdos_down": []},
            }
            for index in indices:
                atom_key = f"atom{index}"
                if atom_key in self.grouped_orbitals:
                    atom_data = self.grouped_orbitals[atom_key]
                    if not np.array(grouped_data[atom]["atom"]["energies"]).any():
                        grouped_data[atom]["atom"]["energies"] = atom_data["atom"][
                            "energies"
                        ]
                    grouped_data[atom]["atom"]["pdos_up"].append(
                        atom_data["atom"]["pdos_up"]
                    )
                    grouped_data[atom]["atom"]["pdos_down"].append(
                        atom_data["atom"]["pdos_down"]
                    )
                    for orbital in ["s", "p", "d", "f"]:
                        if orbital in atom_data:
                            if not np.array(
                                grouped_data[atom][orbital]["energies"]
                            ).any():
                                grouped_data[atom][orbital]["energies"] = atom_data[
                                    orbital
                                ]["energies"]
                            grouped_data[atom][orbital]["pdos_up"].append(
                                atom_data[orbital]["pdos_up"]
                            )
                            grouped_data[atom][orbital]["pdos_down"].append(
                                atom_data[orbital]["pdos_down"]
                            )
            # Convert lists to numpy arrays for easier manipulation
            for orbital in grouped_data[atom]:
                grouped_data[atom][orbital]["pdos_up"] = np.array(
                    grouped_data[atom][orbital]["pdos_up"]
                )
                grouped_data[atom][orbital]["pdos_down"] = np.array(
                    grouped_data[atom][orbital]["pdos_down"]
                )
                grouped_data[atom][orbital]["energies"] = np.array(
                    grouped_data[atom][orbital]["energies"]
                )

            # Sum the PDOS for the grouped atoms
            grouped_data[atom]["atom"]["pdos_up"] = np.sum(
                grouped_data[atom]["atom"]["pdos_up"], axis=0
            )
            grouped_data[atom]["atom"]["pdos_down"] = np.sum(
                grouped_data[atom]["atom"]["pdos_down"], axis=0
            )
            for orbital in ["s", "p", "d", "f"]:
                grouped_data[atom][orbital]["pdos_up"] = np.sum(
                    grouped_data[atom][orbital]["pdos_up"], axis=0
                )
                grouped_data[atom][orbital]["pdos_down"] = np.sum(
                    grouped_data[atom][orbital]["pdos_down"], axis=0
                )

        return grouped_data

    def plot_atomwise_pdos(self, figsize=(4, 4), energy_range=(-5, 5)):
        """Plot the PDOS for each atom."""

        for atom, data in self.pdos_data.items():
            plt.figure(figsize=figsize)
            plt.plot(
                data["atom"]["energies"],
                data["atom"]["pdos_up"],
                label="Up",
                color=colors[0],
            )
            plt.plot(
                data["atom"]["energies"],
                data["atom"]["pdos_down"],
                label="Down",
                color=colors[1],
            )
            plt.title(f"PDOS for {atom}")
            plt.xlabel("Energy (eV)")
            plt.ylabel("PDOS (states/eV)")
            plt.legend(loc="upper right")
            plt.grid(alpha=0.3, linestyle="--")
            plt.xlim(energy_range)
            plt.tight_layout()
            plt.savefig(self.pdos_folder / f"PDOS_{atom}.svg")
            plt.close()

    def plot_atomwise_orbital_pdos(self, figsize=(4, 4), energy_range=(-5, 5)):
        """Plot the PDOS for each atom and orbital."""
        for atom, data in self.grouped_orbitals.items():
            plt.figure(figsize=figsize)
            for orbital, orbital_data in data.items():
                if orbital == "atom":
                    continue
                # Choose a color based on the orbital type
                if "s" in orbital:
                    color = colors[0]
                elif "p" in orbital:
                    color = colors[1]
                elif "d" in orbital:
                    color = colors[2]
                else:
                    color = colors[3]  # Default color

                plt.plot(
                    orbital_data["energies"],
                    orbital_data["pdos_up"],
                    label=f"{orbital}",
                    color=color,
                )
                plt.plot(
                    orbital_data["energies"],
                    orbital_data["pdos_down"],
                    # label=f"{orbital}",
                    color=color,
                )
            plt.title(f"PDOS for {atom}")
            plt.xlabel("Energy (eV)")
            plt.ylabel("PDOS (states/eV)")
            plt.legend(loc="upper right")
            plt.grid(alpha=0.3, linestyle="--")
            plt.xlim(energy_range)
            plt.tight_layout()
            plt.savefig(self.pdos_folder / f"PDOS_{atom}_orbitals.svg")
            plt.close()

    def plot_species_orbital_pdos(self, figsize=(4, 4), energy_range=(-5, 5)):
        """Plot the PDOS for each species and orbital."""
        for species, data in self.grouped_species_orbitals.items():
            plt.figure(figsize=figsize)
            for orbital, orbital_data in data.items():
                if orbital == "atom":
                    continue
                # Choose a color based on the orbital type
                if "s" in orbital:
                    color = colors[0]
                elif "p" in orbital:
                    color = colors[1]
                elif "d" in orbital:
                    color = colors[2]
                else:
                    color = colors[3]  # Default color

                try:
                    plt.plot(
                        orbital_data["energies"],
                        orbital_data["pdos_up"],
                        label=f"{orbital}",
                        color=color,
                    )
                    plt.plot(
                        orbital_data["energies"],
                        orbital_data["pdos_down"],
                        color=color,
                    )
                except ValueError as e:
                    print(f"ValueError for {species} {orbital}. Skipping this orbital.")
                    continue
            plt.title(f"PDOS for {species}")
            plt.xlabel("Energy (eV)")
            plt.ylabel("PDOS (states/eV)")
            plt.legend(loc="upper right")
            plt.grid(alpha=0.3, linestyle="--")
            plt.xlim(energy_range)
            plt.tight_layout()
            plt.savefig(self.pdos_folder / f"PDOS_{species}_orbitals.svg")
            plt.close()


if __name__ == "__main__":
    # Example usage
    root_path = "./"  # Adjust the path as needed
    pdos_folder = "pdos_tetrahedron"  # Folder containing PDOS files

    plotter = PlotPDOS(
        pdos_file="f3gt.PDOS.Tetrahedron", root_path=root_path, pdos_folder=pdos_folder
    )
    plotter.plot_atomwise_pdos()  # Plot PDOS for each atom
    plotter.plot_atomwise_orbital_pdos()  # Plot PDOS for each atom and orbital
    plotter.plot_species_orbital_pdos()  # Plot PDOS for each species and orbital
