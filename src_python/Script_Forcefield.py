"""
Script_Forcefield - Forcefield File Reading and Processing
Corresponds to Script_Forcefield.f90

Handles reading and processing complete forcefield files that contain:
- Atom type definitions
- Bond/angle/torsion type definitions  
- Molecule definitions with connectivity
- Forcefield parameters and types
- Unit specifications

This is the main orchestrator for forcefield file processing.
"""

import sys
import json
from typing import Tuple, Dict, List, Type, Optional, Any, Union
from .VarPrecision import dp

# Import our previously created Script modules
try:
    from .Script_FieldType import script_field_type
    from .Script_BondType import script_bond_type  
    from .Script_AngleType import script_angle_type
    from .Script_TorsionType import script_torsion_type
except ImportError:
    # Fallback functions if modules not available
    def script_field_type(line: str, ff_num: int):
        return None, 0
    def script_bond_type(line: str, bond_num: int):
        return None, 0
    def script_angle_type(line: str, angle_num: int):
        return None, 0
    def script_torsion_type(line: str, torsion_num: int):
        return None, 0

def GetXCommand(line: str, position: int) -> str:
    """Simple implementation of GetXCommand for parsing input"""
    parts = line.split()
    if position <= len(parts):
        return parts[position-1]
    return ""

def FindCommandBlock(lines: List[str], start_line: int, end_command: str) -> int:
    """Find the end of a command block"""
    for i in range(start_line + 1, len(lines)):
        line = lines[i].strip().lower()
        if line == end_command.lower():
            return i
    return len(lines)

def LowerCaseLine(text: str) -> str:
    """Convert line to lowercase"""
    return text.lower()

def LoadFile(filename: str) -> Tuple[List[str], List[int]]:
    """Load file and return lines and line numbers"""
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
        line_numbers = list(range(1, len(lines) + 1))
        return [line.strip() for line in lines], line_numbers
    except FileNotFoundError:
        print(f"ERROR: Could not find forcefield file: {filename}", file=sys.stderr)
        return [], []
    except Exception as e:
        print(f"ERROR: Could not read forcefield file {filename}: {e}", file=sys.stderr)
        return [], []

# =============================================================================
# Main Forcefield File Reading Functions
# =============================================================================

from .Molecule_Definition import Molecule_Type

def script_load_molecule(filename: str, atomtypes):
    json_data = {}
    try:
        with open(filename, 'r') as f:
            json_data = json.load(f)
    except FileNotFoundError:
        print(f"ERROR: Could not find molecule file: {filename}", file=sys.stderr)
    except json.JSONDecodeError:
        print(f"ERROR: Could not parse molecule file {filename}: {e}", file=sys.stderr)
    topofile = json_data.get("topology_file", "")
    if not topofile:
        print(f"ERROR: No topology file specified in {filename}", file=sys.stderr)
        return None
    
    newMolecule = Molecule_Type(topo_file=topofile)

def load_atomtypes(filename: str) -> List[str]:
    """Load atom types from a file"""
    atomtypes = []
    try:
        with open(filename, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    atomtypes.append(parts[0])  # Assuming first part is the atom type
    except FileNotFoundError:
        print(f"ERROR: Could not find atom types file: {filename}", file=sys.stderr)
    return atomtypes


def process_forcefield_type_section(line_store: List[str], i_line: int, line_buffer: int) -> int:
    """Process forcefield type definitions"""
    try:
        if _forcefield_data.energy_calculator:
            print("ERROR: Forcefield type has already been defined", file=sys.stderr)
            return -1
        
        n_items = line_buffer
        _forcefield_data.n_force_fields = n_items
        _forcefield_data.energy_calculator = []
        
        for i in range(1, n_items + 1):
            cur_line = i_line + i
            if cur_line < len(line_store):
                ff_instance, status = script_field_type(line_store[cur_line], i)
                if status < 0:
                    print("ERROR: Unknown forcefield type", file=sys.stderr)
                    return -1
                _forcefield_data.energy_calculator.append(ff_instance)
        
        return 0
        
    except Exception as e:
        print(f"ERROR: Failed to process forcefield type section: {e}", file=sys.stderr)
        return -1

def process_atom_def_section(line_store: List[str], i_line: int, line_buffer: int) -> int:
    """Process atom type definitions"""
    try:
        if _forcefield_data.atom_data:
            print("ERROR: Atom definitions have already been defined", file=sys.stderr)
            return -1
        
        n_items = line_buffer
        _forcefield_data.n_atom_types = n_items
        _forcefield_data.atom_data = []
        
        for i in range(1, n_items + 1):
            cur_line = i_line + i
            if cur_line < len(line_store):
                parts = line_store[cur_line].split()
                if len(parts) >= 2:
                    symbol = parts[0]
                    mass = float(parts[1])
                    atom = AtomData(symbol, mass)
                    _forcefield_data.atom_data.append(atom)
                    print(f"  Atom {i}: {atom}")
        
        return 0
        
    except (ValueError, IndexError) as e:
        print(f"ERROR: Failed to process atom definitions: {e}", file=sys.stderr)
        return -1

def process_bond_def_section(line_store: List[str], i_line: int, line_buffer: int) -> int:
    """Process bond type definitions"""
    try:
        if _forcefield_data.bond_data:
            print("ERROR: Bond definitions have already been defined", file=sys.stderr)
            return -1
        
        n_items = line_buffer
        _forcefield_data.n_bond_types = n_items
        _forcefield_data.bond_data = []
        
        for i in range(1, n_items + 1):
            cur_line = i_line + i
            if cur_line < len(line_store):
                bond_instance, status = script_bond_type(line_store[cur_line], i)
                if status < 0:
                    return -1
                
                bond_data = BondData()
                bond_data.bond_ff = bond_instance
                _forcefield_data.bond_data.append(bond_data)
                print(f"  Bond {i}: {bond_data}")
        
        return 0
        
    except Exception as e:
        print(f"ERROR: Failed to process bond definitions: {e}", file=sys.stderr)
        return -1

def process_angle_def_section(line_store: List[str], i_line: int, line_buffer: int) -> int:
    """Process angle type definitions"""
    try:
        if _forcefield_data.angle_data:
            print("ERROR: Angle definitions have already been defined", file=sys.stderr)
            return -1
        
        n_items = line_buffer
        _forcefield_data.n_angle_types = n_items
        _forcefield_data.angle_data = []
        
        for i in range(1, n_items + 1):
            cur_line = i_line + i
            if cur_line < len(line_store):
                angle_instance, status = script_angle_type(line_store[cur_line], i)
                if status < 0:
                    return -1
                
                angle_data = AngleData()
                angle_data.angle_ff = angle_instance
                _forcefield_data.angle_data.append(angle_data)
                print(f"  Angle {i}: {angle_data}")
        
        return 0
        
    except Exception as e:
        print(f"ERROR: Failed to process angle definitions: {e}", file=sys.stderr)
        return -1

def process_torsion_def_section(line_store: List[str], i_line: int, line_buffer: int) -> int:
    """Process torsion type definitions"""
    try:
        if _forcefield_data.torsion_data:
            print("ERROR: Torsion definitions have already been defined", file=sys.stderr)
            return -1
        
        n_items = line_buffer
        _forcefield_data.n_torsion_types = n_items
        _forcefield_data.torsion_data = []
        
        for i in range(1, n_items + 1):
            cur_line = i_line + i
            if cur_line < len(line_store):
                torsion_instance, status = script_torsion_type(line_store[cur_line], i)
                if status < 0:
                    return -1
                
                torsion_data = TorsionData()
                torsion_data.torsion_ff = torsion_instance
                _forcefield_data.torsion_data.append(torsion_data)
                print(f"  Torsion {i}: {torsion_data}")
        
        return 0
        
    except Exception as e:
        print(f"ERROR: Failed to process torsion definitions: {e}", file=sys.stderr)
        return -1

def process_molecule_section(line_store: List[str], i_line: int, line_buffer: int) -> int:
    """Process molecule definition section"""
    try:
        val = GetXCommand(line_store[i_line], 2)
        mol_type = int(val)
        
        if mol_type > _forcefield_data.n_mol_types:
            print(f"ERROR: Molecule index {mol_type} out of bounds", file=sys.stderr)
            return -1
        
        # Extract the molecule definition block
        mol_def_block = line_store[i_line+1:i_line+line_buffer]
        
        # Process molecule definition
        status = script_read_mol_def(mol_def_block, mol_type)
        if status < 0:
            return -1
        
        # Build molecular graph (placeholder)
        script_build_graph(mol_type)
        
        return 0
        
    except (ValueError, IndexError) as e:
        print(f"ERROR: Failed to process molecule section: {e}", file=sys.stderr)
        return -1

def process_molecule_types_section(line: str) -> int:
    """Process molecule types declaration"""
    try:
        val = GetXCommand(line, 2)
        n_mol_types = int(val)
        
        if _forcefield_data.mol_data:
            print("ERROR: Molecule types have already been defined", file=sys.stderr)
            return -1
        
        _forcefield_data.n_mol_types = n_mol_types
        _forcefield_data.mol_data = [MoleculeData() for _ in range(n_mol_types)]
        
        print(f"Allocated {n_mol_types} molecule types")
        return 0
        
    except (ValueError, IndexError) as e:
        print(f"ERROR: Failed to process molecule types: {e}", file=sys.stderr)
        return -1

def process_units_section(line: str) -> int:
    """Process units section"""
    try:
        status = script_set_units(line)
        return status
        
    except Exception as e:
        print(f"ERROR: Failed to process units: {e}", file=sys.stderr)
        return -1

# =============================================================================
# Molecule Definition Functions
# =============================================================================

def script_read_mol_def(cmd_block: List[str], mol_type: int) -> int:
    """
    Corresponds to Script_ReadMolDef in Fortran.
    Read molecule definition from command block.
    
    Args:
        cmd_block: Lines containing molecule definition
        mol_type: Molecule type index (1-based)
        
    Returns:
        int: Status code (0 = success, negative = error)
    """
    try:
        if mol_type <= 0 or mol_type > len(_forcefield_data.mol_data):
            print(f"ERROR: Invalid molecule type {mol_type}", file=sys.stderr)
            return -1
        
        mol_data = _forcefield_data.mol_data[mol_type - 1]
        n_lines = len(cmd_block)
        line_buffer = 0
        
        i_line = 0
        while i_line < n_lines:
            if line_buffer > 0:
                line_buffer -= 1
                i_line += 1
                continue
            
            line = cmd_block[i_line].strip()
            if not line or line.startswith('#'):
                i_line += 1
                continue
            
            command = GetXCommand(line, 1).lower()
            
            if command == "regrowthtype":
                status = script_regrow_type(line, mol_type)
                if status < 0:
                    return -1
                    
            elif command == "atoms":
                line_buffer = find_command_block_in_list(cmd_block, i_line, "end_atoms") - i_line - 1
                status = process_atoms_section(cmd_block, i_line, line_buffer, mol_data)
                if status < 0:
                    return -1
                    
            elif command == "bonds":
                line_buffer = find_command_block_in_list(cmd_block, i_line, "end_bonds") - i_line - 1
                status = process_bonds_section(cmd_block, i_line, line_buffer, mol_data)
                if status < 0:
                    return -1
                    
            elif command == "angles":
                line_buffer = find_command_block_in_list(cmd_block, i_line, "end_angles") - i_line - 1
                status = process_angles_section(cmd_block, i_line, line_buffer, mol_data)
                if status < 0:
                    return -1
                    
            elif command == "torsion":
                line_buffer = find_command_block_in_list(cmd_block, i_line, "end_torsion") - i_line - 1
                status = process_torsions_section(cmd_block, i_line, line_buffer, mol_data)
                if status < 0:
                    return -1
                    
            elif command == "misc":
                line_buffer = find_command_block_in_list(cmd_block, i_line, "end_misc") - i_line - 1
                status = process_misc_section(cmd_block, i_line, line_buffer, mol_data, mol_type)
                if status < 0:
                    return -1
                    
            else:
                print(f"WARNING: Unknown molecule command '{command}'", file=sys.stderr)
            
            i_line += 1
        
        print(f"Successfully processed molecule {mol_type}: {mol_data}")
        return 0
        
    except Exception as e:
        print(f"ERROR: Failed to process molecule definition: {e}", file=sys.stderr)
        return -1

def find_command_block_in_list(lines: List[str], start_line: int, end_command: str) -> int:
    """Find end of command block within a list of lines"""
    for i in range(start_line + 1, len(lines)):
        line = lines[i].strip().lower()
        if line == end_command.lower():
            return i
    return len(lines)

def process_atoms_section(cmd_block: List[str], i_line: int, line_buffer: int, mol_data: MoleculeData) -> int:
    """Process atoms section of molecule definition"""
    try:
        if mol_data.atom_type:
            print("ERROR: Atoms have already been defined for this molecule", file=sys.stderr)
            return -1
        
        n_items = line_buffer
        mol_data.n_atoms = n_items
        mol_data.atom_type = []
        
        for i in range(1, n_items + 1):
            cur_line = i_line + i
            if cur_line < len(cmd_block):
                parts = cmd_block[cur_line].split()
                if parts:
                    atom_type = int(parts[0])
                    mol_data.atom_type.append(atom_type)
        
        print(f"    Atoms: {mol_data.atom_type}")
        return 0
        
    except (ValueError, IndexError) as e:
        print(f"ERROR: Failed to process atoms section: {e}", file=sys.stderr)
        return -1

def process_bonds_section(cmd_block: List[str], i_line: int, line_buffer: int, mol_data: MoleculeData) -> int:
    """Process bonds section of molecule definition"""
    try:
        if mol_data.bond:
            print("ERROR: Bonds have already been defined for this molecule", file=sys.stderr)
            return -1
        
        n_items = line_buffer
        mol_data.n_bonds = n_items
        mol_data.bond = []
        
        # Initialize connectivity arrays
        mol_data.n_atm_bonds = [0] * mol_data.n_atoms
        mol_data.atm_bonds = [[0] * mol_data.n_atoms for _ in range(n_items)]
        
        for i in range(1, n_items + 1):
            cur_line = i_line + i
            if cur_line < len(cmd_block):
                parts = cmd_block[cur_line].split()
                if len(parts) >= 3:
                    bond_type = int(parts[0])
                    atm1 = int(parts[1]) - 1  # Convert to 0-based
                    atm2 = int(parts[2]) - 1  # Convert to 0-based
                    
                    bond = BondData()
                    bond.bond_type = bond_type
                    bond.mem1 = atm1
                    bond.mem2 = atm2
                    mol_data.bond.append(bond)
                    
                    # Update connectivity
                    mol_data.n_atm_bonds[atm1] += 1
                    mol_data.n_atm_bonds[atm2] += 1
        
        print(f"    Bonds: {len(mol_data.bond)} bonds defined")
        return 0
        
    except (ValueError, IndexError) as e:
        print(f"ERROR: Failed to process bonds section: {e}", file=sys.stderr)
        return -1

def process_angles_section(cmd_block: List[str], i_line: int, line_buffer: int, mol_data: MoleculeData) -> int:
    """Process angles section of molecule definition"""
    try:
        if mol_data.angle:
            print("ERROR: Angles have already been defined for this molecule", file=sys.stderr)
            return -1
        
        n_items = line_buffer
        mol_data.n_angles = n_items
        mol_data.angle = []
        
        # Initialize connectivity arrays
        mol_data.n_atm_angles = [0] * mol_data.n_atoms
        mol_data.atm_angles = [[0] * mol_data.n_atoms for _ in range(n_items)]
        
        for i in range(1, n_items + 1):
            cur_line = i_line + i
            if cur_line < len(cmd_block):
                parts = cmd_block[cur_line].split()
                if len(parts) >= 4:
                    angle_type = int(parts[0])
                    atm1 = int(parts[1]) - 1  # Convert to 0-based
                    atm2 = int(parts[2]) - 1
                    atm3 = int(parts[3]) - 1
                    
                    angle = AngleData()
                    angle.angle_type = angle_type
                    angle.mem1 = atm1
                    angle.mem2 = atm2
                    angle.mem3 = atm3
                    mol_data.angle.append(angle)
                    
                    # Update connectivity
                    mol_data.n_atm_angles[atm1] += 1
                    mol_data.n_atm_angles[atm2] += 1
                    mol_data.n_atm_angles[atm3] += 1
        
        print(f"    Angles: {len(mol_data.angle)} angles defined")
        return 0
        
    except (ValueError, IndexError) as e:
        print(f"ERROR: Failed to process angles section: {e}", file=sys.stderr)
        return -1

def process_torsions_section(cmd_block: List[str], i_line: int, line_buffer: int, mol_data: MoleculeData) -> int:
    """Process torsions section of molecule definition"""
    try:
        if mol_data.torsion:
            print("ERROR: Torsions have already been defined for this molecule", file=sys.stderr)
            return -1
        
        n_items = line_buffer
        mol_data.n_tors = n_items
        mol_data.torsion = []
        
        # Initialize connectivity arrays
        mol_data.n_atm_torsions = [0] * mol_data.n_atoms
        mol_data.atm_torsions = [[0] * mol_data.n_atoms for _ in range(n_items)]
        
        for i in range(1, n_items + 1):
            cur_line = i_line + i
            if cur_line < len(cmd_block):
                parts = cmd_block[cur_line].split()
                if len(parts) >= 5:
                    torsion_type = int(parts[0])
                    atm1 = int(parts[1]) - 1  # Convert to 0-based
                    atm2 = int(parts[2]) - 1
                    atm3 = int(parts[3]) - 1
                    atm4 = int(parts[4]) - 1
                    
                    torsion = TorsionData()
                    torsion.torsion_type = torsion_type
                    torsion.mem1 = atm1
                    torsion.mem2 = atm2
                    torsion.mem3 = atm3
                    torsion.mem4 = atm4
                    mol_data.torsion.append(torsion)
                    
                    # Update connectivity
                    mol_data.n_atm_torsions[atm1] += 1
                    mol_data.n_atm_torsions[atm2] += 1
                    mol_data.n_atm_torsions[atm3] += 1
                    mol_data.n_atm_torsions[atm4] += 1
        
        print(f"    Torsions: {len(mol_data.torsion)} torsions defined")
        return 0
        
    except (ValueError, IndexError) as e:
        print(f"ERROR: Failed to process torsions section: {e}", file=sys.stderr)
        return -1

def process_misc_section(cmd_block: List[str], i_line: int, line_buffer: int, mol_data: MoleculeData, mol_type: int) -> int:
    """Process miscellaneous interactions section"""
    try:
        if mol_data.misc_data:
            print("ERROR: Misc interactions have already been defined for this molecule", file=sys.stderr)
            return -1
        
        n_items = line_buffer
        mol_data.n_misc = n_items
        mol_data.misc_data = []
        
        for i in range(1, n_items + 1):
            cur_line = i_line + i
            if cur_line < len(cmd_block):
                misc_instance, status = script_misc_type(cmd_block[cur_line], mol_type, i)
                if status < 0:
                    return -1
                
                misc_data = MiscData()
                misc_data.misc_ff = misc_instance
                mol_data.misc_data.append(misc_data)
        
        print(f"    Misc: {len(mol_data.misc_data)} misc interactions defined")
        return 0
        
    except Exception as e:
        print(f"ERROR: Failed to process misc section: {e}", file=sys.stderr)
        return -1

# =============================================================================
# Utility Functions
# =============================================================================

def script_regrow_type(line: str, mol_type: int) -> int:
    """Process regrowth type specification"""
    # Placeholder - would implement regrowth type processing
    print(f"  Regrowth type for molecule {mol_type}: {line}")
    return 0

def script_misc_type(line: str, mol_type: int, misc_num: int) -> Tuple[object, int]:
    """Process miscellaneous interaction type"""
    # Placeholder - would implement misc type processing
    print(f"  Misc {misc_num} for molecule {mol_type}: {line}")
    return None, 0

def script_build_graph(mol_type: int):
    """Build molecular connectivity graph"""
    # Placeholder - would implement graph building
    print(f"  Building graph for molecule {mol_type}")

def script_set_units(line: str) -> int:
    """
    Process units specification.
    
    Args:
        line: Input line with units command
        
    Returns:
        int: Status code (0 = success, negative = error)
    """
    try:
        command = GetXCommand(line, 2).lower()
        unit_type = GetXCommand(line, 3)
        
        if command == "angle":
            _forcefield_data.ang_unit = find_angular_unit(unit_type)
        elif command == "energy":
            _forcefield_data.eng_unit = find_energy_unit(unit_type)
        elif command in ["distance", "length"]:
            _forcefield_data.len_unit = find_length_unit(unit_type)
        else:
            print(f"ERROR: Unknown unit type '{command}'", file=sys.stderr)
            return -1
        
        print(f"Set {command} unit to {unit_type}")
        return 0
        
    except Exception as e:
        print(f"ERROR: Failed to process units: {e}", file=sys.stderr)
        return -1

def find_angular_unit(unit_str: str) -> float:
    """Find angular unit conversion factor"""
    unit_map = {
        "degrees": 1.0,
        "radians": 57.2957795,  # 180/pi
        "deg": 1.0,
        "rad": 57.2957795
    }
    return unit_map.get(unit_str.lower(), 1.0)

def find_energy_unit(unit_str: str) -> float:
    """Find energy unit conversion factor"""
    unit_map = {
        "kcal/mol": 1.0,
        "kj/mol": 0.239006,
        "hartree": 627.509,
        "ev": 23.0609,
        "k": 0.0019872  # Boltzmann constant
    }
    return unit_map.get(unit_str.lower(), 1.0)

def find_length_unit(unit_str: str) -> float:
    """Find length unit conversion factor"""
    unit_map = {
        "angstrom": 1.0,
        "bohr": 0.529177,
        "nm": 10.0,
        "pm": 0.01,
        "a": 1.0,
        "au": 0.529177
    }
    return unit_map.get(unit_str.lower(), 1.0)

# =============================================================================
# Public Interface Functions
# =============================================================================

def get_forcefield_data() -> ForceFieldData:
    """Get the global forcefield data"""
    return _forcefield_data

def reset_forcefield_data():
    """Reset the global forcefield data"""
    global _forcefield_data
    _forcefield_data = ForceFieldData()

def validate_forcefield_data() -> bool:
    """Validate the loaded forcefield data"""
    try:
        if not _forcefield_data.atom_data:
            print("WARNING: No atom types defined")
        
        if not _forcefield_data.mol_data:
            print("WARNING: No molecule types defined")
        
        print(f"Forcefield validation:")
        print(f"  Atom types: {_forcefield_data.n_atom_types}")
        print(f"  Bond types: {_forcefield_data.n_bond_types}")
        print(f"  Angle types: {_forcefield_data.n_angle_types}")
        print(f"  Torsion types: {_forcefield_data.n_torsion_types}")
        print(f"  Molecule types: {_forcefield_data.n_mol_types}")
        print(f"  Force fields: {_forcefield_data.n_force_fields}")
        
        return True
        
    except Exception as e:
        print(f"ERROR: Forcefield validation failed: {e}", file=sys.stderr)
        return False

if __name__ == "__main__":
    # Test forcefield file reading
    import sys
    
    if len(sys.argv) > 1:
        filename = sys.argv[1]
        print(f"Testing forcefield file reading: {filename}")
        
        ff_data, status = script_read_field_file(filename)
        if status == 0:
            print("SUCCESS: Forcefield file processed successfully")
            validate_forcefield_data()
        else:
            print(f"FAILED: Status code {status}")
    else:
        print("Available forcefield processing functions:")
        print("  script_read_field_file(filename)")
        print("  get_forcefield_data()")
        print("  validate_forcefield_data()")
        print("  reset_forcefield_data()")
        
        # Test with sample data
        print("\nTesting with sample data:")
        reset_forcefield_data()
        _forcefield_data.n_atom_types = 2
        _forcefield_data.atom_data = [
            AtomData("C", 12.01),
            AtomData("H", 1.008)
        ]
        validate_forcefield_data()
