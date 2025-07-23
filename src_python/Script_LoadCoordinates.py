from .Box_SimpleBox import SimpleBox
from .Box_CubeBox import CubeBox


#===============================================================                
def load_coords(coordinatefile, moldata):
    
    boxtypes = {
        'nobox': SimpleBox,
        'cube': CubeBox,
        "cubic": CubeBox
    }
    
    f = open(coordinatefile, 'r')
    lines = f.readlines()
    f.close()
    
    #Remove lines with comments
    lines = [line for line in lines if not line.startswith('#')]
    
    #Remove empty lines
    lines = [line for line in lines if line.strip() != '']
    
    # Pull the first non-empty line to determine the box type
    first_line = lines[0].strip().lower().split()
    

    box = boxtypes.get(first_line[0])
    if box is None:
        raise ValueError(f"Unknown box type: {first_line[0]}")
    box = box(moldata) if box else None
    
    #Extra columns after give the box dimensions
    if box is not None:
        if len(first_line) > 1:
            try:
                box.load_dimension([float(dim) for dim in first_line[1:]])
            except ValueError as e:
                raise ValueError(f"Invalid dimensions in the first line: {e}")
    box.load_coordinate(lines[1:], moldata)
    
    
    return box
#===============================================================
