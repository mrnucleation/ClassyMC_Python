"""
Monte Carlo Moves Input Processing
Corresponds to Script_MCMoves.f90

Handles parsing input to instantiate the correct Monte Carlo move types.
"""

import sys
from typing import List, Dict, Optional, Tuple
from .VarPrecision import dp

# Import all available Monte Carlo moves
from .MC_Move_MolTranslation import MolTranslate
from .MC_Move_IsoVol import IsoVol
from .MC_Move_BasicSwap import BasicSwap
from .MC_Move_Delete import Delete
from .MC_Move_AVBMC import AVBMC

def GetXCommand(line: str, position: int) -> str:
    """Simple implementation to get command at position"""
    parts = line.split()
    if position <= len(parts):
        return parts[position-1]
    return ""

def script_mc_moves(line: str, move_num: int) -> Tuple[object, float, int]:
    """
    Corresponds to Script_MCMoves
    Parse Monte Carlo move type from input line and instantiate appropriate move.
    
    Args:
        line: Input line containing move type and probability
        move_num: Move number for indexing
        
    Returns:
        tuple: (move_instance, probability, status_code)
               status_code: 0 = success, negative = error
    """
    line_stat = 0
    
    try:
        # Parse move type and probability
        parts = line.strip().split()
        if len(parts) < 2:
            print(f"ERROR: Invalid move specification: {line}", file=sys.stderr)
            return None, 0.0, -1
        
        command = parts[0].lower()
        probability = float(parts[1])
        
        # Validate probability
        if probability < 0.0:
            print(f"ERROR: Move probability must be non-negative: {probability}", file=sys.stderr)
            return None, 0.0, -1
        
        # Instantiate the appropriate move type
        move_instance = None
        
        if command == "moltranslation" or command == "moltrans":
            move_instance = MolTranslate(None)  # BoxArray will be set later
            
        elif command == "isovol" or command == "volume":
            move_instance = IsoVol()
            
        elif command == "basicswap" or command == "swap":
            move_instance = BasicSwap()
            
        elif command == "delete" or command == "debugdelete":
            move_instance = Delete()
            
        elif command == "avbmc":
            move_instance = AVBMC()
            
        # Add more moves as they become available:
        # elif command == "atomtranslation" or command == "atomtrans":
        #     move_instance = AtomTranslate()
        # elif command == "atomexchange":
        #     move_instance = AtomExchange()
        # elif command == "cbmc":
        #     move_instance = CBMC()
        # elif command == "anisovol":
        #     move_instance = AnisoVol()
        # elif command == "particleexchange":
        #     move_instance = ParticleExchange()
        # elif command == "planerotate":
        #     move_instance = PlaneRotate()
        # elif command == "planetranslation":
        #     move_instance = PlaneTranslate()
        # elif command == "planeatomtranslation":
        #     move_instance = PlaneAtomTranslate()
        # elif command == "volexchange":
        #     move_instance = VolExchange()
        # elif command == "ubswap":
        #     move_instance = UBSwap()
        
        else:
            print(f"ERROR: Unknown Monte Carlo move type '{command}'", file=sys.stderr)
            return None, 0.0, -1
        
        # Initialize the move if instantiation was successful
        if move_instance is not None:
            move_instance.constructor()
            print(f"Monte Carlo move {move_num}: {command} (probability: {probability})")
        
        return move_instance, probability, line_stat
        
    except (ValueError, IndexError) as e:
        print(f"ERROR: Invalid move specification format: {line}", file=sys.stderr)
        print(f"Exception: {e}", file=sys.stderr)
        return None, 0.0, -1

def register_additional_moves():
    """
    Register additional move types as they become available.
    This function can be expanded as more moves are converted from Fortran.
    """
    # This is a placeholder for future move registrations
    # As more moves are converted, they can be added to the registry
    pass

def list_available_moves() -> List[str]:
    """
    List all available Monte Carlo move types.
    
    Returns:
        List[str]: Available move type names
    """
    available_moves = [
        "moltranslation",
        "moltrans", 
        "isovol",
        "volume",
        "basicswap",
        "swap",
        "delete",
        "debugdelete",
        "avbmc"
    ]
    
    # Future moves to be added:
    future_moves = [
        "atomtranslation",
        "atomtrans",
        "atomexchange", 
        "cbmc",
        "anisovol",
        "particleexchange",
        "planerotate",
        "planetranslation",
        "planeatomtranslation", 
        "volexchange",
        "ubswap"
    ]
    
    return available_moves

def process_moves_input(lines: List[str]) -> Tuple[List[object], List[float], int]:
    """
    Process multiple lines of Monte Carlo move input.
    
    Args:
        lines: List of input lines containing move specifications
        
    Returns:
        tuple: (list_of_moves, list_of_probabilities, status_code)
    """
    moves = []
    probabilities = []
    
    for i, line in enumerate(lines):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
            
        move_instance, prob, status = script_mc_moves(line, i + 1)
        if status < 0:
            return [], [], status
            
        if move_instance is not None:
            moves.append(move_instance)
            probabilities.append(prob)
    
    return moves, probabilities, 0

def validate_move_probabilities(probabilities: List[float]) -> bool:
    """
    Validate that move probabilities are reasonable.
    
    Args:
        probabilities: List of move probabilities
        
    Returns:
        bool: True if probabilities are valid
    """
    if not probabilities:
        return False
    
    total_prob = sum(probabilities)
    if total_prob <= 0.0:
        print("ERROR: Total move probability must be positive", file=sys.stderr)
        return False
    
    return True

def normalize_move_probabilities(probabilities: List[float]) -> List[float]:
    """
    Normalize move probabilities to sum to 1.0.
    
    Args:
        probabilities: List of move probabilities
        
    Returns:
        List[float]: Normalized probabilities
    """
    if not probabilities:
        return []
    
    total_prob = sum(probabilities)
    if total_prob <= 0.0:
        return probabilities
    
    return [p / total_prob for p in probabilities] 