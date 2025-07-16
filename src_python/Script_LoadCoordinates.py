from Box_SimpleBox import SimpleBox
from Script_BoxType import Script_BoxType


def load_coords(box_indx, coordinatefile, molTypes):
    
    box = None
    with open(coordinatefile, "r") as infile:
        for line in infile:
            if line.strip().startswith("#") or len(line.strip()):
                continue
            commands = line.split()
            if commands[0] == 'boxtype':
                box = Script_BoxType(commands[1])
            
            if box is None:
                raise ValueError("Box type must be defined before loading coordinates")
            
            if commands[0] == 'dimension':
                pass
            elif commands[0] == 'mol':
                pass
            elif commands[0] == 'molmin':
                pass
            elif commands[0] == 'molmax':
                pass
            else:
                print(f"Unknown command: {commands[0]}")
                
                
            
                
            
            
    