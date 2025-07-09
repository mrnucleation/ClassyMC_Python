# Monte Carlo Move Translation Status

This document tracks the progress of translating Monte Carlo move classes from Fortran to Python.

## Completed Translations

### Base Classes
- ✅ `Template_MCMove.py` - Base class for all Monte Carlo moves
- ✅ `CoordinateTypes.py` - Displacement types used by moves

### Move Classes
- ✅ `MC_Move_MolTranslation.py` - Molecule translation moves (already existed)
- ✅ `MC_Move_IsoVol.py` - Isochoric volume change moves
- ✅ `MC_Move_BasicSwap.py` - Basic molecule swap moves (add/remove)
- ✅ `MC_Move_Delete.py` - Molecule deletion moves
- ⚠️ `MC_Move_AVBMC.py` - Aggregation Volume Bias Monte Carlo (partial, has dependencies)

## Fortran Move Classes Still to Translate

### Simple Moves (Recommended Priority)
- `Move_MC_AtomTranslation.f90` - Atom translation moves
- `Move_MC_AtomExchange.f90` - Atom exchange moves
- `Move_MC_PlaneAtomTranslate.f90` - Plane-constrained atom translation
- `Move_MC_PlaneRotate.f90` - Plane rotation moves
- `Move_MC_PlaneTranslate.f90` - Plane translation moves

### Complex Moves (Lower Priority)
- `Move_MC_AnisoVol.f90` - Anisotropic volume changes
- `Move_MC_CBMC.f90` - Configurational Bias Monte Carlo
- `Move_MC_EB_AVBMC.f90` - Energy-biased AVBMC
- `Move_MC_ParticleExchange.f90` - Particle exchange between boxes
- `Move_MC_ThermoLambda.f90` - Thermodynamic integration moves
- `Move_MC_UBSwap.f90` - Unbiased swap moves
- `Move_MC_VolExchange.f90` - Volume exchange moves

### Genetic Algorithm Moves
- `Move_GA_AtomTranslation.f90` - Genetic algorithm atom translation

## Key Dependencies Needed

### Missing Modules
- `Common_MolInfo.py` - Molecular information and data structures
- `RandomGen.py` - Random number generation utilities
- `BoxData.py` - Box array management
- `MolConstruct.py` - Molecule construction utilities

### Required Functions
- `Generate_UnitSphere()` - Generate random unit vector
- `FindMolecule()` - Find molecule from raw index
- `ListRNG()` - Weighted random selection
- Various molecule construction methods

## Translation Guidelines

### Structure
Each move class should:
1. Inherit from `MCMove`
2. Implement the `full_move()` method
3. Include proper error handling and validation
4. Follow the same interface as the Fortran version

### Key Methods to Implement
- `__init__()` - Initialize move parameters
- `constructor()` - Set up arrays and data structures
- `full_move()` - Perform the actual move
- `maintenance()` - Tune parameters if needed
- `prologue()` - Initialize before simulation
- `epilogue()` - Finalize after simulation
- `process_io()` - Handle input commands

### Coordinate Types
Use the appropriate displacement types from `CoordinateTypes.py`:
- `Displacement` - Basic atomic moves
- `Addition` - Molecule addition moves
- `Deletion` - Molecule deletion moves
- `OrthoVolChange` - Orthogonal volume changes
- `TriVolChange` - Triclinic volume changes
- `AtomExchange` - Atom exchange moves

## Testing Strategy

### Unit Tests Needed
1. Test each move class initialization
2. Test move execution with simple systems
3. Test acceptance/rejection logic
4. Test parameter tuning
5. Test input/output processing

### Integration Tests Needed
1. Test moves with different box types
2. Test moves with different force fields
3. Test moves with constraints
4. Test moves in multi-box simulations

## Next Steps

1. **Complete missing dependencies** - Create the required modules and functions
2. **Translate simple moves first** - Focus on atom translation and exchange moves
3. **Add proper error handling** - Ensure robust error checking
4. **Create unit tests** - Test each move class thoroughly
5. **Documentation** - Add detailed docstrings and examples

## Notes

- The existing `MC_Move_MolTranslation.py` needs to be updated to use the new base class
- Some moves have complex dependencies that need to be resolved
- The sampling system integration needs to be tested
- Performance optimization may be needed for production use 