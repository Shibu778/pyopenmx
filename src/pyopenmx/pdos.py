# Script to plot partial DOS from collinear calculations
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib import style
import sys
import json

style_file = Path(__file__).parent / "pyopenmx.mplstyle"
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
        dostype="Tetrahedron",
    ):
        self.dostype = dostype
        self.root_path = Path(root_path)
        self.pdos_folder = self.root_path / pdos_folder
        self.pdos_filename, self.pdos_data = self.get_pdos_files()
        self.spinpolarization = spinpolarization
        self.pdos_data = self.extract_all_pdos()
        self.grouped_orbitals = self.group_orbitals()
        self.grouped_species_orbitals = self.group_species_and_orbitals()

    def get_pdos_files(self):
        """Get all PDOS files in the specified folder."""
        pdos_files = list(self.pdos_folder.glob(f"*.PDOS.{self.dostype}.*"))
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
        """
        Notes: In these files, the first and second columns are energy 
        in eV and DOS (eV$^{-1}$) or PDOS (eV$^{-1}$), and the third column 
        is the integrated DOS or PDOS. If a spin-polarized calculation using 
        'LSDA-CA', 'LSDA-PW', or 'GGA-PBE' is employed in the SCF calculation, 
        the second and third columns in these files correspond to DOS or PDOS 
        for up and down spin states, respectively, and the fourth and fifth 
        columns are the corresponding integrated values.
        Ref: https://www.openmx-square.org/openmx_man3.9/node70.html
        """
        data = np.loadtxt(pdos_file)
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
                self.pdos_data[atom]["atom"] = self.load_pdos(
                    self.pdos_folder / pdos_file
                )
            elif len(f) == 5:
                atom = f[3]
                orbital = f[4]
                self.pdos_data[atom][orbital] = self.load_pdos(
                    self.pdos_folder / pdos_file
                )
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

    def group_species(self, list_of_atoms={"Fe": [1, 2, 3], "Ge": [4], "Te": [5, 6]}):
        """Group PDOS data for a list of species."""
        grouped_data = {}
        for species, indices in list_of_atoms.items():
            grouped_data[species] = {
                "atom": {"energies": [], "pdos_up": [], "pdos_down": []},
                "s1": {"energies": [], "pdos_up": [], "pdos_down": []},
                "p1": {"energies": [], "pdos_up": [], "pdos_down": []},
                "p2": {"energies": [], "pdos_up": [], "pdos_down": []},
                "p3": {"energies": [], "pdos_up": [], "pdos_down": []},
                "d1": {"energies": [], "pdos_up": [], "pdos_down": []},
                "d2": {"energies": [], "pdos_up": [], "pdos_down": []},
                "d3": {"energies": [], "pdos_up": [], "pdos_down": []},
                "d4": {"energies": [], "pdos_up": [], "pdos_down": []},
                "d5": {"energies": [], "pdos_up": [], "pdos_down": []},
                "f1": {"energies": [], "pdos_up": [], "pdos_down": []},
                "f2": {"energies": [], "pdos_up": [], "pdos_down": []},
                "f3": {"energies": [], "pdos_up": [], "pdos_down": []},
                "f4": {"energies": [], "pdos_up": [], "pdos_down": []},
                "f5": {"energies": [], "pdos_up": [], "pdos_down": []},
                "f6": {"energies": [], "pdos_up": [], "pdos_down": []},
                "f7": {"energies": [], "pdos_up": [], "pdos_down": []},
            }
            for index in indices:
                atom_key = f"atom{index}"
                if atom_key in self.pdos_data:
                    atom_data = self.pdos_data[atom_key]
                    if not np.array(grouped_data[species]["atom"]["energies"]).any():
                        grouped_data[species]["atom"]["energies"] = atom_data["atom"][
                            "energies"
                        ]
                    grouped_data[species]["atom"]["pdos_up"].append(
                        atom_data["atom"]["pdos_up"]
                    )
                    grouped_data[species]["atom"]["pdos_down"].append(
                        atom_data["atom"]["pdos_down"]
                    )
            # Convert lists to numpy arrays for easier manipulation
            grouped_data[species]["atom"]["pdos_up"] = np.array(
                grouped_data[species]["atom"]["pdos_up"]
            )
            grouped_data[species]["atom"]["pdos_down"] = np.array(
                grouped_data[species]["atom"]["pdos_down"]
            )
            grouped_data[species]["atom"]["energies"] = np.array(
                grouped_data[species]["atom"]["energies"]
            )

            # Sum the PDOS for the grouped species
            grouped_data[species]["atom"]["pdos_up"] = np.sum(
                grouped_data[species]["atom"]["pdos_up"], axis=0
            )
            grouped_data[species]["atom"]["pdos_down"] = np.sum(
                grouped_data[species]["atom"]["pdos_down"], axis=0
            )

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

    def dump_data(self, filename="pdos_data.json", data=None):
        if data is None:
            data = self.pdos_data

        with open(filename, "w") as f:
            json.dump(data, f, indent=4, default=str)

        np_filename = filename.replace(".json", ".npy")
        np.save(np_filename, data)

    def load_data(self, filename="pdos_data.json"):
        with open(filename, "r") as f:
            data = json.load(f)
        self.pdos_data = data
        self.grouped_orbitals = self.group_orbitals()
        self.grouped_species_orbitals = self.group_species_and_orbitals()
        return data

    def plot_custom_pdos(
        self,
        atom_list,
        legend,
        orbital_list=None,
        total_dos=None,
        output_path="custom_pdos.svg",
        color_list=None,
        kwargs={},
    ):
        """Plot custom PDOS for a list of atoms and orbitals.

        Args:
            atom_list (list[list]): It sums data for each sublist and plots it.
                e.g. [['atom1'], ['atom2']] or [['atom1', 'atom2'], ['atom3']]
                First case will plot PDOS for atom1 and atom2 separately,
                second case will sum PDOS for atom1 and atom2 and plot it together.
            orbital_list (list[list], optional): List of orbitals to sum for each sublist.
                e.g. [['s'], ['p1', 'p2', 'p3']] or [['s'], ['s', 'p', 'd']]
                First case will plot s orbital for first sublist and p1+p2+p3 for second sublist.
                Second case will sum s+p+d for each sublist and plot it.
                If None, it will plot the total PDOS for each sublist.
            legend (list): Legend for each sublist.
            total_dos (path): Path to the total DOS file to plot alongside PDOS.

            Note: Length of legend must be equal to length of atom_list and orbital_list (if provided).
        """
        data_to_plot = []
        for i, atoms in enumerate(atom_list):
            summed_data = {"energies": None, "pdos_up": 0, "pdos_down": 0}
            for atom in atoms:
                if atom not in self.pdos_data:
                    print(f"Atom {atom} not found in PDOS data. Skipping.")
                    continue
                atom_data = self.pdos_data[atom]
                if orbital_list and i < len(orbital_list):
                    orbitals = orbital_list[i]
                    for orbital in orbitals:
                        if orbital in atom_data:
                            orbital_data = atom_data[orbital]
                            if summed_data["energies"] is None:
                                summed_data["energies"] = orbital_data["energies"]
                            summed_data["pdos_up"] += orbital_data["pdos_up"]
                            summed_data["pdos_down"] += orbital_data["pdos_down"]
                        else:
                            print(
                                f"Orbital {orbital} not found for atom {atom}. Skipping."
                            )
                else:
                    # Sum total PDOS from 'atom' key
                    total_data = atom_data["atom"]
                    if summed_data["energies"] is None:
                        summed_data["energies"] = total_data["energies"]
                    summed_data["pdos_up"] += total_data["pdos_up"]
                    summed_data["pdos_down"] += total_data["pdos_down"]
            data_to_plot.append(summed_data)
        if total_dos:
            total_dos_data = self.load_pdos(total_dos)
            data_to_plot.append(total_dos_data)
            legend.append("Total DOS")

        plt.figure(figsize=(6, 4))
        if color_list is None:
            color_list = colors
        for i, data in enumerate(data_to_plot):
            plt.plot(
                data["energies"],
                data["pdos_up"],
                label=f"{legend[i]}",
                color=color_list[i],
                **kwargs,
            )
            plt.plot(
                data["energies"],
                data["pdos_down"],
                color=color_list[i],
                **kwargs,
            )
        plt.xlabel("Energy (eV)")
        plt.ylabel("DOS (states/eV)")
        # plt.title("Custom PDOS Plot")
        plt.legend(loc="upper right")
        plt.grid(alpha=0.3, linestyle="--")
        plt.xlim(-5, 5)
        plt.tight_layout()
        plt.savefig(output_path)


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
    # Custom PDOS plot example
    atom_list = [["atom1", "atom2"], ["atom3"]]
    orbital_list = [
        ["d1", "d2", "d3", "d4", "d5"],
        ["d1", "d2", "d3", "d4", "d5"],
    ]
    labels = ["Fe-out (3d)", "Fe-in (3d)"]
    plotter.plot_custom_pdos(
        atom_list,
        labels,
        orbital_list,
        total_dos="f3gt.DOS.Tetrahedron",
        output_path="custom_pdos.svg",
    )
