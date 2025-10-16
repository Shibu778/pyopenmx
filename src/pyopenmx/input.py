# Class to read input data


class InputData:
    def __init__(self):
        self.system_dir = None
        self.system_name = None
        self.data_path = None
        self.species_number = None
        self.level_of_stdout = None
        self.level_of_fileout = None
        self.definition_of_atomic_species = {}
        self.atom_number = None
        self.species_coordinate_unit = None
        self.species_coordinate = {}
        self.unitvectors_unit = None
        self.unitvectors = []

        # SCF parameters (defaults from provided list)
        self.scf_mixing_every_pulay = 1  # scf.Mixing.EveryPulay
        self.scf_mixing_start_pulay = 30  # scf.Mixing.StartPulay
        self.scf_criterion = 1.0e-06  # scf.criterion (final value given)
        self.scf_kgrid = (15, 15, 1)  # scf.Kgrid
        self.scf_max_iter = 200  # scf.maxIter
        self.scf_mixing_history = 50  # scf.Mixing.History
        self.scf_electronic_temperature = 1000  # scf.ElectronicTemperature
        self.scf_energycutoff = 150  # scf.energycutoff
        self.scf_init_mixing_weight = 0.001  # scf.Init.Mixing.Weight
        self.scf_min_mixing_weight = 0.001  # scf.Min.Mixing.Weight
        self.scf_max_mixing_weight = 0.3  # scf.Max.Mixing.Weight
        self.scf_xctype = "GGA-PBE"  # scf.XcType
        self.scf_spin_polarization = True  # scf.SpinPolarization (On -> True)
        self.scf_eigenvalue_solver = "band"  # scf.EigenvalueSolver
        self.scf_mixing_type = "RMM-DIISK"  # scf.Mixing.Type

        # DOS parameters
        self.dos_fileout = True  # Dos.fileout (on -> True)
        self.dos_erange = (-25.0, 20.0)  # Dos.Erange
        self.dos_kgrid = (20, 20, 1)
