#!/usr/bin/env python3

#### >>>> Description of Code <<<< ####
__description__= """
Class for OpenMX DFT calculation
"""
__author__ = "Cha"
__email__ = "caca518@kaist.ac.kr"
__version__ = 0.0
#### >>>> Description of Code End <<<< ####

#### >>>> Modules and constant <<<< ####
from typing import Dict, List

#### >>>> Modules and constant End <<<< ####

#### >>>> Class <<<< ####
class Test:
    message = "Hello, World!"
    def print(self):
        print(self.message)
        
class OmxInput:
    """ Class for OpenMX input file"""
    def  __init__(self, datfile):
        self.sysname = None
        self.datfile = datfile
        self.options = self.get_options(self.datfile)
        self.structure = self.get_structure(self.datfile)
    def get_structure(self, datfile):
        """ Get structure from datfile"""
        structure = Structure()
        structure.read_datfile(datfile, spin_option = self.options.spin_polarization)
        return structure
    def get_options(self, datfile):
        """ Get options from datfile"""
        options = OmxOptions()
        options.read_datfile(datfile)
        return options
    
    def write_datfile(self, filepath = "./"):
        """ Write datfile"""
        with open(filepath + self.options.sysname + '.dat', 'w') as f:
            # Header
            f.write(f'System.CurrrentDirectory         ./    # default=./\n')
            f.write(f'DATA.PATH                        {OPENMX_DATA_PATH}\n')
            f.write(f'System.Name                      {self.sysname}\n')
            f.write(f'level.of.stdout                  {str(self.options.level_of_stdout)}\n')
            f.write(f'level.of.fileout                 {str(self.options.level_of_fileout)}\n')
            
            # Restart
            f.write(f'scf.restart                      {"on" if self.options.is_restart else "off"}\n') 
            # Hamiltonian
            f.write(f'HS.fileout                       {"on" if self.options.is_hamiltonian_fileout else "off"}\n')
            
            # Structure
            f.write(f'Species.Number                   {len(set(position[0] for position in self.structure.positions))}\n')
            f.write('<Definition.of.Atomic.Species\n')
            for element in set(position[0] for position in self.structure.positions):
                f.write(f'  {element}   {VPS_DIC[element][self.options.vpsoption]}    {VPS_DIC[element][0]}\n')
            f.write('Definition.of.Atomic.Species>\n')

            f.write(f'Atoms.Number                     {len(self.structure.positions)}\n')
            f.write(f'Atoms.SpeciesAndCoordinates.Unit FRAC\n')
            f.write('<Atoms.SpeciesAndCoordinates\n')
            for index, (position, spin) in enumerate(zip(self.structure.positions, self.structure.spin)):
                position_string = f'  {index+1} {position[0]}  {" ".join(str(i) for i in position[1])}'
                if self.options.spin_polarization in ['on', 'off']:
                    spin_string = f'   {format(float(VPS_DIC[position[0]][1])/2 + spin[0]/2, ".6f")} {format(float(VPS_DIC[position[0]][1])/2 - spin[0]/2, ".6f")}  off'
                elif self.options.spin_polarization == 'NC':
                    spin_string = f'   {format(float(VPS_DIC[position[0]][1])/2 + spin[0]/2, ".6f")} {format(float(VPS_DIC[position[0]][1])/2 - spin[0]/2, ".6f")}  {str(spin[1][0])} {str(spin[1][1])}  {str(spin[2][0])} {str(spin[2][1])}  off'
                f.write(position_string + spin_string + '\n')
            f.write('Atoms.SpeciesAndCoordinates>\n')

            f.write(f'Atoms.UnitVectors.Unit           {self.structure.unit}\n')
            f.write('<Atoms.UnitVectors\n')
            for vector in self.structure.unit_cell:
                f.write(f'  {" ".join(format(i, ".6f") for i in vector)}\n')
            f.write('Atoms.UnitVectors>\n')
            # pass
            
            # SCF
            f.write(f'scf.XcType                       {self.options.xctype}\n')
            f.write(f'scf.ElectronicTemperature        {self.options.electronic_temperature}\n')
            f.write(f'scf.maxIter                      {self.options.max_iteration}\n')
            f.write(f'scf.EigenvalueSolver             {self.options.eigen_solver}\n')
            f.write(f'scf.energycutoff                 {self.options.energy_cutoff}\n')
            f.write(f'scf.criterion                    {self.options.scf_criteria}\n')
            f.write(f'scf.Kgrid                        {" ".join(str(i) for i in self.options.kpoints)}\n')
            f.write(f'scf.Generation.Kpoint            {self.options.kpoints_generator}\n')

            # Spin
            f.write(f'scf.SpinPolarization             {self.options.spin_polarization}\n')
            f.write(f'scf.SpinOrbit.Coupling           {"on" if self.options.is_spin_orbit_coupling else "off"}\n')
            f.write(f'scf.Constraint.NC.Spin           {self.options.constraint_nc_spin}\n')
            f.write(f'scf.Constraint.NC.Spin.v         {self.options.constraint_nc_spin_value}\n')
            if self.options.so_factor is not None:
                f.write('<scf.SO.factor\n')
                for atom, orbitals in self.options.so_factor.items():
                    f.write(f'  {atom} {" ".join(str(i) for i in orbitals)}\n')
                f.write('scf.SO.factor>\n')

            # Hubbard U
            f.write(f'scf.Hubbard.U                    {"on" if self.options.is_dft_u_correction else "off"}\n')
            f.write(f'scf.DFTU.Type                    {self.options.dftu_type}\n')
            if self.options.hubbard_u_occupation is not None:
                f.write(f'scf.Hubbard.Occupation           {self.options.hubbard_u_occupation}\n')
            if self.options.hubbard_u_value is not None:
                f.write('<Hubbard.U.values\n')
                for atom, orbitals in self.options.hubbard_u_value.items():
                    string = f'  {atom}'
                    for orbital, U in orbitals.items():
                        string += f' {orbital} {U}'
                    f.write(string + '\n')
                f.write('Hubbard.U.values>\n')
            if self.options.hund_j_value is not None:
                f.write('<Hund.J.values\n')
                for atom, orbitals in self.options.hund_j_value.items():
                    string = f'  {atom}'
                    for orbital, J in orbitals.items():
                        string += f' {orbital} {J}'
                    f.write(string + '\n')
                f.write('Hund.J.values>\n')
            
            # Band
            f.write(f'Band.dispersion                 {"on" if self.options.is_draw_band else "off"}\n')
            if self.options.band_n_kpath is not None:
                f.write(f'Band.Nkpath                     {self.options.band_n_kpath}\n')
            if self.options.band_kpath is not None:
                f.write('<Band.kpath\n')
                for path in self.options.band_kpath:
                    f.write(f'  {path[0]}  {" ".join(str(i) for i in path[1])} {" ".join(str(i) for i in path[2])} {path[3]} {path[4]}\n')
                f.write('Band.kpath>\n')

            # DOS
            f.write(f'Dos.fileout                     {"on" if self.options.is_draw_dos else "off"}\n')
            if self.options.dos_energy_range is not None:
                f.write(f'Dos.Erange                      {" ".join(str(i) for i in self.options.dos_energy_range)}\n')
            if self.options.dos_kgrid is not None:
                f.write(f'Dos.Kgrid                       {" ".join(str(i) for i in self.options.dos_kgrid)}\n')

            # Mixing
            f.write(f'scf.Mixing.Type                 {self.options.mixing_type}\n')
            f.write(f'scf.Init.Mixing.Weight          {self.options.init_mixing_weight}\n')
            f.write(f'scf.Min.Mixing.Weight           {self.options.min_mixing_weight}\n')
            f.write(f'scf.Max.Mixing.Weight           {self.options.max_mixing_weight}\n')
            f.write(f'scf.Kerker.factor               {self.options.kerker_factor}\n')
            f.write(f'scf.Mixing.History              {self.options.mixing_history}\n')


class Structure:
    """ Class for structure"""
    def __init__(self):
        self.unit: str = None # Ang or AU
        self.unit_cell: List[List[float], List[float], List[float]] = None # [vec_a, vec_b, vec_c]
        self.positions: List[str, List[float, float, float]] = None # {atom: [x, y, z]}
        self.spin: List[float, List[float, float], List[float, float]] = None # [magmom, spin[theta, phi], orbital[theta, phi]]
    def read_datfile(self, datfile, spin_option: str):
        unit_vector_process_line = False
        positions_process_line = False
        with open(datfile, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if line.startswith('#'):
                    continue
                if 'Atoms.UnitVectors.Unit' in line:
                    self.unit = line.split()[1]
                if '<Atoms.UnitVectors' in line:
                    unit_vector_process_line = True
                    self.unit_cell = []
                    continue
                if 'Atoms.UnitVectors>' in line:
                    unit_vector_process_line = False
                    continue
                if unit_vector_process_line:
                    self.unit_cell.append([float(i) for i in line.split()])
                
                if '<Atoms.SpeciesAndCoordinates' in line:
                    positions_process_line = True
                    self.positions = []
                    self.spin = []
                    continue
                if 'Atoms.SpeciesAndCoordinates>' in line:
                    positions_process_line = False
                    continue
                if positions_process_line:
                    components = line.split()
                    self.positions.append([components[1], [float(i) for i in components[2:5]]])
                    if spin_option in ['on', 'off']:
                        self.spin.append([float(components[5]) - float(components[6]), [0, 0], [0, 0]])
                    elif spin_option == 'NC':
                        self.spin.append([float(components[5]) - float(components[6]), [float(components[7]), float(components[8])], [float(components[9]), float(components[10])]])

class OmxOptions:
    """ Class for OpenMX options"""
    def __init__(self, vpsoption=4):
        # Header
        self.sysname: str = None
        self.level_of_stdout: int = 1
        self.level_of_fileout: int = 1

        # Restart
        self.is_restart: bool = False

        # Hamiltonian
        self.is_hamiltonian_fileout: bool = False

        # Pseudo potential and basis set
        self.vpsoption: int = vpsoption  # 2: Quick, 3: Standard, 4: Precise

        # SCF
        self.xctype: str = None
        self.electronic_temperature: float = 300.0 # K
        self.max_iteration: int = 40
        self.eigen_solver: str = None
        self.energy_cutoff: float = 150.0 # Rydberg
        self.scf_criteria: float = 1.0e-6 # Hartree
        self.kpoints: List[int, int, int] = None
        self.kpoints_generator: str = 'regular'

        # Mixing
        self.mixing_type: str = None  # Simple, GR-Pulay, RMM-DIIS, Kerker, RMM-DIISK
        self.init_mixing_weight: float = 0.3 # 0-1, simple, GR-Pulay, RMM-DIIS, Kerker and RMM-DIISK
        self.min_mixing_weight: float = 0.01 # 0-1, simple, Kerker
        self.max_mixing_weight: float = 0.4  # 0-1, simple, Kerker
        self.kerker_factor: float = 1.0 # Kerker, RMM-DIISK
        self.mixing_history: int = 6  # GR-Pulay, RMM-DIIS, Kerker, RMM-DIISK
        self.mixing_start_pulay: int = 6 # GR-pulay, RMM-DIIS, Kerker, RMM-DIISK
        self.mixing_every_pulay: int = 5 # support only on mixing_type RMM-DIISK

        # Spin
        self.spin_polarization: str = 'off' # on, off, nc
        self.is_spin_orbit_coupling: bool = False
        self.constraint_nc_spin: str = 'off' # on, on2, off
        self.constraint_nc_spin_value: float = 0.0 # eV
        self.so_factor: Dict[str, Dict[str, float]] = None # {atom: {orbital: factor}

        # Hubbard U
        self.is_dft_u_correction: bool = False
        self.dftu_type: int = 1 # 1: simplified(Dudarev), 2: General
        self.hubbard_u_occupation: str = None # onsite, full, dual
        self.hubbard_u_value: Dict[str, Dict[str, float]] = None # {atom: {orbital: U}} 
        self.hund_j_value: Dict[str, Dict[str, float]] = None # {atom: {orbital: J}} 

        # Band
        self.is_draw_band: bool = False
        self.band_n_kpath: int = None
        self.band_kpath: List[List[int, List[float], List[float], str, str]] = None # [kgrid, start_point, end_point, label_start, label_end]

        # DOS
        self.is_draw_dos: bool = False
        self.dos_energy_range: List[float, float] = None # [min, max]
        self.dos_kgrid: List[int, int, int] = None 
    
    def read_datfile(self, datfile):
        so_factor_process_line = False
        hubbard_u_process_line = False
        hund_j_process_line = False
        band_kpath_process_line = False
        with open(datfile, 'r') as f:
            lines = f.readlines()
            for line in lines:
                # skip comment line
                if line.startswith('#'):
                    continue

                # read header
                if 'System.Name' in line:
                    self.sysname = line.split()[1]
                if 'level.of.stdout' in line:
                    self.level_of_stdout = int(line.split()[1])
                if 'level.of.fileout' in line:
                    self.level_of_fileout = int(line.split()[1])

                # read restart option and Hamiltonian fileout option
                if 'scf.restart' in line:
                    self.is_restart = True if line.split()[1] == 'on' else False
                if 'HS.fileout' in line:
                    self.is_hamiltonian_fileout = True if line.split()[1] == 'on' else False

                # read SCF options
                if 'scf.XcType' in line:
                    self.xctype = line.split()[1]
                if 'scf.ElectronicTemperature' in line:
                    self.electronic_temperature = float(line.split()[1])
                if 'scf.maxIter' in line:
                    self.max_iteration = int(line.split()[1])
                if 'scf.EigenvalueSolver' in line:
                    self.eigen_solver = line.split()[1]
                if 'scf.energycutoff' in line:
                    self.energy_cutoff = float(line.split()[1])
                if 'scf.criterion' in line: 
                    self.scf_criteria = float(line.split()[1])
                if 'scf.Kgrid' in line:
                    self.kpoints = [int(i) for i in line.split()[1:4]]
                if 'scf.Generation.Kpoint' in line:
                    self.kpoints_generator = line.split()[1]

                # read mixing options
                if 'scf.Mixing.Type' in line:
                    self.mixing_type = line.split()[1]
                if 'scf.Init.Mixing.Weight' in line:
                    self.init_mixing_weight = float(line.split()[1])
                if 'scf.Min.Mixing.Weight' in line:
                    self.min_mixing_weight = float(line.split()[1])
                if 'scf.Max.Mixing.Weight' in line:
                    self.max_mixing_weight = float(line.split()[1])
                if 'scf.Kerker.factor' in line:
                    self.kerker_factor = float(line.split()[1])
                if 'scf.Mixing.History' in line:
                    self.mixing_history = int(line.split()[1])
                if 'scf.Mixing.StartPulay' in line:
                    self.mixing_start_pulay = int(line.split()[1])
                if 'scf.Mixing.EveryPulay' in line:
                    self.mixing_every_pulay = int(line.split()[1])

                # read spin options
                if 'scf.SpinPolarization' in line:
                    self.spin_polarization = line.split()[1]
                # read so_factor
                if '<scf.SO.factor' in line:
                    so_factor_process_line = True
                    self.so_factor = {}
                    continue
                if 'scf.SO.factor>' in line:
                    so_factor_process_line = False
                    continue
                if so_factor_process_line:
                    components = line.split()
                    atom = components[0]
                    orbitals = {components[i]: float(components[i+1]) for i in range(1, len(components), 2)}
                    self.so_factor[atom] = orbitals
                if 'scf.SpinOrbit.Coupling' in line:
                    self.is_spin_orbit_coupling = True if line.split()[1] == 'on' else False
                if 'scf.Constraint.NC.Spin ' in line:
                    self.constraint_nc_spin = line.split()[1]
                if 'scf.Constraint.NC.Spin.v ' in line:
                    self.constraint_nc_spin_value = float(line.split()[1])
                if 'scf.Hubbard.U' in line:
                    self.is_dft_u_correction = True if line.split()[1] == 'on' else False
                if 'scf.DFTU.Type' in line:
                    self.dftu_type = int(line.split()[1])
                if 'scf.Hubbard.Occupation' in line:
                    self.hubbard_u_occupation = line.split()[1]
                # read Hubbard U value
                if '<Hubbard.U.values' in line:
                    hubbard_u_process_line = True
                    self.hubbard_u_value = {}
                    continue
                if 'Hubbard.U.values>' in line:
                    hubbard_u_process_line = False
                    continue
                if hubbard_u_process_line:
                    components = line.split()
                    atom = components[0]
                    orbitals = {components[i]: float(components[i+1]) for i in range(1, len(components), 2)}
                    self.hubbard_u_value[atom] = orbitals
                # read Hubbard J value
                if '<Hund.J.values' in line:
                    hund_j_process_line = True
                    self.hund_j_value = {}
                    continue
                if 'Hund.J.values>' in line:
                    hund_j_process_line = False
                    continue
                if hund_j_process_line:
                    components = line.split()
                    atom = components[0]
                    orbitals = {components[i]: float(components[i+1]) for i in range(1, len(components), 2)}
                    self.hund_j_value[atom] = orbitals
                                # read band options
                if 'Band.dispersion' in line:
                    self.is_draw_band = True if line.split()[1] == 'on' else False
                if 'Band.Nkpath' in line:
                    self.band_n_kpath = int(line.split()[1])
                # read band kpath
                if '<Band.kpath' in line:
                    band_kpath_process_line = True
                    self.band_kpath = []
                    continue
                if 'Band.kpath>' in line:
                    band_kpath_process_line = False
                    continue
                if band_kpath_process_line:
                    components = line.split()
                    kgrid = int(components[0])
                    start_point = [float(components[i]) for i in range(1, 4)]
                    end_point = [float(components[i]) for i in range(4, 7)]
                    label_start = components[7]
                    label_end = components[8]
                    self.band_kpath.append([kgrid, start_point, end_point, label_start, label_end])
                # read DOS options
                if 'Dos.fileout' in line:
                    self.is_draw_dos = True if line.split()[1] == 'on' else False
                if 'Dos.Erange' in line:
                    self.dos_energy_range = [float(line.split()[1]), float(line.split()[2])]
                if 'Dos.Kgrid' in line:
                    self.dos_kgrid = [int(line.split()[1]), int(line.split()[2]), int(line.split()[3])]



class OmxOutput:
    """ Class for OpenMX output file"""
    def __init__(self, outfile):
        with open(outfile, 'r') as f:
            _content = f.read()
        self.magmom = self.get_magmom(_content)
        self.energy = self.get_energy(_content)
        self.convergence = self.get_convergence(_content)

    def get_magmom(self, _content):
        """ Get magnetic moment"""
        lines = _content.split('\n')
        for line in lines:
            if 'Total spin moment (muB)' in line:
                _total_magmom = float(line.split()[-1])
                break
        else:
            raise ValueError('Magnetic moment not found, keyword was [Total spin moment (muB)]')
        
        return _total_magmom
    

    def get_energy(self, _content):
        """ Get energy"""
        # process content to get energy
        lines = _content.split('\n')
        for line in lines:
            if 'Uele.' in line:
                total_energy = float(line.split()[-1])
                break
        else:
            raise ValueError('Energy not found, keyword was [Uele.]')
        return total_energy

    def get_convergence(self, _content):
        """ Get convergence information: [normrd, duele]"""
        lines = _content.split('\n')
        last_two_steps = []
        for line in reversed(lines):
            if 'SCF=' in line:
                last_two_steps.append(line)
                if len(last_two_steps) == 2:
                    break
        else:
            raise ValueError('SCF steps not found, keyword was [SCF=]')
        normrd = float(last_two_steps[0].split()[3])
        duele = abs(float(last_two_steps[1].split()[5]) - float(last_two_steps[0].split()[5]))
        return [normrd, duele]

#### >>>> Class End <<<< ####

#### >>>> Constants <<<< ####
OPENMX_DATA_PATH = '/home/cha/program/openmx3.9/DFT_DATA19'
OPENMX_BIN = '/home/cha/program/openmx3.9/work/openmx'
VPS_DIC = {
    # element : [VPS, Valence electrons, Quick, Standard, Precise]
    'E' : ['E', '0.0 ', 'Kr10.0-s1p1   ', 'Kr10.0-s2p1d1 ', 'Kr10.0-s2p2d1f1'],
    'H' : ['H_PBE19 ', '1.0 ', 'H5.0-s2       ', 'H6.0-s2p1     ', 'H7.0-s2p2d1'],
    'He': ['He_PBE19 ', '2.0 ', 'He8.0-s1p1    ', 'He8.0-s2p1    ', 'He10.0-s2p2d1'],
    'Li': ['Li_PBE19 ', '3.0 ', 'Li8.0-s3p1    ', 'Li8.0-s3p2    ', 'Li8.0-s3p2d1'],
    'Be': ['Be_PBE19 ', '2.0 ', 'Be7.0-s2p1    ', 'Be7.0-s2p2    ', 'Be7.0-s3p2d1'],
    'B' : ['B_PBE19 ', '3.0 ', 'B7.0-s2p2     ', 'B7.0-s2p2d1   ', 'B7.0-s3p2d2'],
    'C' : ['C_PBE19 ', '4.0 ', 'C6.0-s2p2     ', 'C6.0-s2p2d1   ', 'C6.0-s3p2d2'],
    'N' : ['N_PBE19 ', '5.0 ', 'N6.0-s2p2     ', 'N6.0-s2p2d1   ', 'N6.0-s3p2d2'],
    'O' : ['O_PBE19 ', '6.0 ', 'O6.0-s2p2     ', 'O6.0-s2p2d1   ', 'O6.0-s3p2d2'],
    'F' : ['F_PBE19 ', '7.0 ', 'F6.0-s2p2     ', 'F6.0-s2p2d1   ', 'F6.0-s3p3d2f1'],
    'Ne': ['Ne_PBE19 ', '8.0 ', 'Ne9.0-s2p2    ', 'Ne9.0-s2p2d1  ', 'Ne9.0-s3p2d2'],
    'Na': ['Na_PBE19 ', '9.0 ', 'Na9.0-s3p2    ', 'Na9.0-s3p2d1  ', 'Na9.0-s3p2d2'],
    'Mg': ['Mg_PBE19 ', '8.0 ', 'Mg9.0-s2p2    ', 'Mg9.0-s3p2d1  ', 'Mg9.0-s3p2d2'],
    'Al': ['Al_PBE19 ', '3.0 ', 'Al7.0-s2p1d1  ', 'Al7.0-s2p2d1  ', 'Al7.0-s3p2d2'],
    'Si': ['Si_PBE19 ', '4.0 ', 'Si7.0-s2p1d1  ', 'Si7.0-s2p2d1  ', 'Si7.0-s3p3d2'],
    'P' : ['P_PBE19 ', '5.0 ', 'P7.0-s2p2d1   ', 'P7.0-s2p2d1f1 ', 'P7.0-s3p2d2f1'],
    'S' : ['S_PBE19 ', '6.0 ', 'S7.0-s2p2d1   ', 'S7.0-s2p2d1f1 ', 'S7.0-s3p2d2f1'],
    'Cl': ['Cl_PBE19 ', '7.0 ', 'Cl7.0-s2p2d1  ', 'Cl7.0-s2p2d1f1', 'Cl7.0-s3p2d2f1'],
    'Ar': ['Ar_PBE19 ', '8.0 ', 'Ar9.0-s2p2d1  ', 'Ar9.0-s2p2d1f1', 'Ar9.0-s3p2d2f1'],
    'K' : ['K_PBE19 ', '9.0 ', 'K10.0-s3p2    ', 'K10.0-s3p2d1  ', 'K10.0-s3p2d2'],
    'Ca': ['Ca_PBE19 ', '10.0', 'Ca9.0-s3p2    ', 'Ca9.0-s3p2d1  ', 'Ca9.0-s3p2d2'],
    'Sc': ['Sc_PBE19 ', '11.0', 'Sc9.0-s2p2d1  ', 'Sc9.0-s3p2d1  ', 'Sc9.0-s3p2d2'],
    'Ti': ['Ti_PBE19 ', '12.0', 'Ti7.0-s2p2d1  ', 'Ti7.0-s3p2d1  ', 'Ti7.0-s3p2d2f1'],
    'V' : ['V_PBE19 ', '13.0', 'V6.0-s2p2d1   ', 'V6.0-s3p2d1   ', 'V6.0-s3p2d2f1'],
    'Cr': ['Cr_PBE19 ', '14.0', 'Cr6.0-s2p2d1  ', 'Cr6.0-s3p2d1  ', 'Cr6.0-s3p2d2f1'],
    'Mn': ['Mn_PBE19 ', '15.0', 'Mn6.0-s2p2d1  ', 'Mn6.0-s3p2d1  ', 'Mn6.0-s3p2d2f1'],
#vps_dic['Fe_H'] = ['Fe_PBE19H', '16.0', 'Fe5.5H-s2p2d1 ', 'Fe5.5H-s3p2d1 ', 'Fe5.5H-s3p2d2f1']
#vps_dic['Fe_S'] = ['Fe_PBE19S', '14.0', 'Fe6.0S-s2p2d1 ', 'Fe6.0S-s3p2d1 ', 'Fe6.0S-s3p2d2f1']
#vps_dic['Co_H'] = ['Co_PBE19H', '17.0', 'Co6.0H-s2p2d1 ', 'Co6.0H-s3p2d1 ', 'Co6.0H-s3p2d2f1']
#vps_dic['Co_S'] = ['Co_PBE19S', '15.0', 'Co6.0S-s2p2d1 ', 'Co6.0S-s3p2d1 ', 'Co6.0S-s3p2d2f1']
#vps_dic['Ni_H'] = ['Ni_PBE19H', '18.0', 'Ni6.0H-s2p2d1 ', 'Ni6.0H-s3p2d1 ', 'Ni6.0H-s3p2d2f1']
#vps_dic['Ni_S'] = ['Ni_PBE19S', '16.0', 'Ni6.0S-s2p2d1 ', 'Ni6.0S-s3p2d1 ', 'Ni6.0S-s3p2d2f1']
#vps_dic['Cu_H'] = ['Cu_PBE19H', '19.0', 'Cu6.0H-s2p2d1 ', 'Cu6.0H-s3p2d1 ', 'Cu6.0H-s3p2d2f1']
#vps_dic['Cu_S'] = ['Cu_PBE19S', '11.0', 'Cu6.0S-s2p1d1 ', 'Cu6.0S-s3p2d1 ', 'Cu6.0S-s3p2d2f1']
#vps_dic['Zn_H'] = ['Zn_PBE19H', '20.0', 'Zn6.0H-s2p2d1 ', 'Zn6.0H-s3p2d1 ', 'Zn6.0H-s3p2d2f1']
#vps_dic['Zn_S'] = ['Zn_PBE19S', '12.0', 'Zn6.0S-s2p1d1 ', 'Zn6.0S-s3p2d1 ', 'Zn6.0S-s3p2d2f1']
    'Fe': ['Fe_PBE19H', '16.0', 'Fe5.5H-s2p2d1 ', 'Fe5.5H-s3p2d1 ', 'Fe5.5H-s3p2d2f1'],
    'Co': ['Co_PBE19H', '17.0', 'Co6.0H-s2p2d1 ', 'Co6.0H-s3p2d1 ', 'Co6.0H-s3p2d2f1'],
    'Ni': ['Ni_PBE19H', '18.0', 'Ni6.0H-s2p2d1 ', 'Ni6.0H-s3p2d1 ', 'Ni6.0H-s3p2d2f1'],
    'Cu': ['Cu_PBE19H', '19.0', 'Cu6.0H-s2p2d1 ', 'Cu6.0H-s3p2d1 ', 'Cu6.0H-s3p2d2f1'],
    'Zn': ['Zn_PBE19H', '20.0', 'Zn6.0H-s2p2d1 ', 'Zn6.0H-s3p2d1 ', 'Zn6.0H-s3p2d2f1'],
    'Ga': ['Ga_PBE19 ', '13.0', 'Ga7.0-s2p2d1  ', 'Ga7.0-s3p2d2  ', 'Ga7.0-s3p2d2f1'],
    'Ge': ['Ge_PBE19 ', '4.0 ', 'Ge7.0-s2p1d1  ', 'Ge7.0-s3p2d2  ', 'Ge7.0-s3p2d2f1'],
    'As': ['As_PBE19 ', '15.0', 'As7.0-s3p2d1  ', 'As7.0-s3p2d2  ', 'As7.0-s3p2d2f1'],
    'Se': ['Se_PBE19 ', '6.0 ', 'Se7.0-s3p2d1  ', 'Se7.0-s3p2d2  ', 'Se7.0-s3p2d2f1'],
    'Br': ['Br_PBE19 ', '7.0 ', 'Br7.0-s3p2d1  ', 'Br7.0-s3p2d2  ', 'Br7.0-s3p2d2f1'],
    'Kr': ['Kr_PBE19 ', '8.0 ', 'Kr10.0-s2p2d1 ', 'Kr10.0-s3p2d2 ', 'Kr10.0-s3p2d2f1'],
    'Rb': ['Rb_PBE19 ', '9.0 ', 'Rb11.0-s2p2d1 ', 'Rb11.0-s3p2d2 ', 'Rb11.0-s3p2d2f1'],
    'Sr': ['Sr_PBE19 ', '10.0', 'Sr10.0-s2p2d1 ', 'Sr10.0-s3p2d2 ', 'Sr10.0-s3p3d2f1'],
    'Y' : ['Y_PBE19 ', '11.0', 'Y10.0-s3p2d1  ', 'Y10.0-s3p2d2  ', 'Y10.0-s3p3d2f1'],
    'Zr': ['Zr_PBE19 ', '12.0', 'Zr7.0-s3p2d1  ', 'Zr7.0-s3p2d2  ', 'Zr7.0-s3p2d2f1'],
    'Nb': ['Nb_PBE19 ', '13.0', 'Nb7.0-s3p2d1  ', 'Nb7.0-s3p2d2  ', 'Nb7.0-s3p2d2f1'],
    'Mo': ['Mo_PBE19 ', '14.0', 'Mo7.0-s3p2d1  ', 'Mo7.0-s3p2d2  ', 'Mo7.0-s3p2d2f1'],
    'Tc': ['Tc_PBE19 ', '15.0', 'Tc7.0-s3p2d1  ', 'Tc7.0-s3p2d2  ', 'Tc7.0-s3p2d2f1'],
    'Ru': ['Ru_PBE19 ', '14.0', 'Ru7.0-s3p2d1  ', 'Ru7.0-s3p2d2  ', 'Ru7.0-s3p2d2f1'],
    'Rh': ['Rh_PBE19 ', '15.0', 'Rh7.0-s3p2d1  ', 'Rh7.0-s3p2d2  ', 'Rh7.0-s3p2d2f1'],
    'Pd': ['Pd_PBE19 ', '16.0', 'Pd7.0-s3p2d1  ', 'Pd7.0-s3p2d2  ', 'Pd7.0-s3p2d2f1'],
    'Ag': ['Ag_PBE19 ', '17.0', 'Ag7.0-s3p2d1  ', 'Ag7.0-s3p2d2  ', 'Ag7.0-s3p2d2f1'],
    'Cd': ['Cd_PBE19 ', '12.0', 'Cd7.0-s3p2d1  ', 'Cd7.0-s3p2d2  ', 'Cd7.0-s3p2d2f1'],
    'In': ['In_PBE19 ', '13.0', 'In7.0-s3p2d1  ', 'In7.0-s3p2d2  ', 'In7.0-s3p2d2f1'],
    'Sn': ['Sn_PBE19 ', '14.0', 'Sn7.0-s3p2d1  ', 'Sn7.0-s3p2d2  ', 'Sn7.0-s3p2d2f1'],
    'Sb': ['Sb_PBE19 ', '15.0', 'Sb7.0-s3p2d1  ', 'Sb7.0-s3p2d2  ', 'Sb7.0-s3p2d2f1'],
    'Te': ['Te_PBE19 ', '16.0', 'Te7.0-s3p2d2  ', 'Te7.0-s3p2d2f1', 'Te7.0-s3p3d2f1'],
    'I' : ['I_PBE19 ', '17.0', 'I7.0-s3p2d2   ', 'I7.0-s3p2d2f1 ', 'I7.0-s3p3d2f1'],
    'Xe': ['Xe_PBE19 ', '18.0', 'Xe9.0-s3p2d2  ', 'Xe9.0-s3p2d2f1', 'Xe9.0-s3p3d2f1'],
    'Cs': ['Cs_PBE19 ', '19.0', 'Cs11.0-s2p2d1 ', 'Cs11.0-s3p2d2 ', 'Cs11.0-s3p2d2f1'],
    'Ba': ['Ba_PBE19 ', '20.0', 'Ba10.0-s2p2d1 ', 'Ba10.0-s3p2d2 ', 'Ba10.0-s3p3d2f1'],
    'La': ['La_PBE19 ', '11.0', 'La10.0-s3p2d1 ', 'La10.0-s3p2d2 ', 'La10.0-s3p3d2f1'],
    'Ce': ['Ce_PBE19 ', '12.0', 'Ce10.0-s3p2d1 ', 'Ce10.0-s3p2d2 ', 'Ce10.0-s3p3d2f1'],
    'Pr': ['Pr_PBE19 ', '13.0', 'Pr10.0-s3p2d1 ', 'Pr10.0-s3p2d2 ', 'Pr10.0-s3p3d2f1'],
    'Nd': ['Nd_PBE19 ', '14.0', 'Nd10.0-s3p2d1 ', 'Nd10.0-s3p2d2 ', 'Nd10.0-s3p3d2f1'],
    'Pm': ['Pm_PBE19 ', '15.0', 'Pm10.0-s3p2d1 ', 'Pm10.0-s3p2d2 ', 'Pm10.0-s3p3d2f1'],
    'Sm': ['Sm_PBE19 ', '16.0', 'Sm10.0-s3p2d1 ', 'Sm10.0-s3p2d2 ', 'Sm10.0-s3p3d2f1'],
    'Eu': ['Eu_PBE19 ', '17.0', 'Eu10.0-s3p2d1 ', 'Eu10.0-s3p2d2 ', 'Eu10.0-s3p3d2f1'],
    'Gd': ['Gd_PBE19 ', '18.0', 'Gd10.0-s3p2d1 ', 'Gd10.0-s3p2d2 ', 'Gd10.0-s3p3d2f1'],
    'Tb': ['Tb_PBE19 ', '19.0', 'Tb10.0-s3p2d1 ', 'Tb10.0-s3p2d2 ', 'Tb10.0-s3p3d2f1'],
    'Dy': ['Dy_PBE19 ', '20.0', 'Dy10.0-s3p2d1 ', 'Dy10.0-s3p2d2 ', 'Dy10.0-s3p3d2f1'],
    'Ho': ['Ho_PBE19 ', '11.0', 'Ho10.0-s3p2d1 ', 'Ho10.0-s3p2d2 ', 'Ho10.0-s3p3d2f1'],
    'Er': ['Er_PBE19 ', '12.0', 'Er10.0-s3p2d1 ', 'Er10.0-s3p2d2 ', 'Er10.0-s3p3d2f1'],
    'Tm': ['Tm_PBE19 ', '13.0', 'Tm10.0-s3p2d1 ', 'Tm10.0-s3p2d2 ', 'Tm10.0-s3p3d2f1'],
    'Yb': ['Yb_PBE19 ', '14.0', 'Yb10.0-s3p2d1 ', 'Yb10.0-s3p2d2 ', 'Yb10.0-s3p3d2f1'],
    'Lu': ['Lu_PBE19 ', '15.0', 'Lu10.0-s3p2d1 ', 'Lu10.0-s3p2d2 ', 'Lu10.0-s3p3d2f1'],
    'Hf': ['Hf_PBE19 ', '16.0', 'Hf7.0-s3p2d1  ', 'Hf7.0-s3p2d2  ', 'Hf7.0-s3p2d2f1'],
    'Ta': ['Ta_PBE19 ', '17.0', 'Ta7.0-s3p2d1  ', 'Ta7.0-s3p2d2  ', 'Ta7.0-s3p2d2f1'],
    'W' : ['W_PBE19 ', '14.0', 'W7.0-s3p2d1   ', 'W7.0-s3p2d2   ', 'W7.0-s3p2d2f1'],
    'Re': ['Re_PBE19 ', '15.0', 'Re7.0-s3p2d1  ', 'Re7.0-s3p2d2  ', 'Re7.0-s3p2d2f1'],
    'Os': ['Os_PBE19 ', '16.0', 'Os7.0-s3p2d1  ', 'Os7.0-s3p2d2  ', 'Os7.0-s3p2d2f1'],
    'Ir': ['Ir_PBE19 ', '17.0', 'Ir7.0-s3p2d1  ', 'Ir7.0-s3p2d2  ', 'Ir7.0-s3p2d2f1'],
    'Pt': ['Pt_PBE19 ', '18.0', 'Pt7.0-s3p2d1  ', 'Pt7.0-s3p2d2  ', 'Pt7.0-s3p2d2f1'],
    'Au': ['Au_PBE19 ', '19.0', 'Au7.0-s3p2d1  ', 'Au7.0-s3p2d2  ', 'Au7.0-s3p2d2f1'],
    'Hg': ['Hg_PBE19 ', '12.0', 'Hg7.0-s3p2d1  ', 'Hg7.0-s3p2d2  ', 'Hg7.0-s3p2d2f1'],
    'Tl': ['Tl_PBE19 ', '13.0', 'Tl7.0-s3p2d1  ', 'Tl7.0-s3p2d2  ', 'Tl7.0-s3p2d2f1'],
    'Pb': ['Pb_PBE19 ', '14.0', 'Pb7.0-s3p2d1  ', 'Pb7.0-s3p2d2  ', 'Pb7.0-s3p2d2f1'],
    'Bi': ['Bi_PBE19 ', '15.0', 'Bi7.0-s3p2d1  ', 'Bi7.0-s3p2d2  ', 'Bi7.0-s3p2d2f1']
}

#### >>>> Constants End <<<< ####