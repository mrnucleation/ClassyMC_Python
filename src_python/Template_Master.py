"""
Master Template Class
Corresponds to Template_Master.f90

Provides the base class with common lifecycle methods for all simulation components.
"""

class ClassyClass:
    """
    Base class for all simulation components.
    Provides common lifecycle methods that all derived classes can override.
    """
    
    def __init__(self):
        # Common properties that might be needed by all classes
        self.perMove = False
        self.IOUnit = -1
        self.UpdateFreq = -1
        
    def prologue(self):
        """Called at the beginning of simulation - override in subclasses"""
        pass
    
    def epilogue(self):
        """Called at the end of simulation - override in subclasses"""
        pass
    
    def maintenance(self):
        """Called periodically during simulation - override in subclasses"""
        pass
    
    def update(self):
        """Called to update internal state - override in subclasses"""
        pass
    
    def safety_check(self):
        """Called to verify system integrity - override in subclasses"""
        pass
    
    def process_io(self, line):
        """Process input/output commands - override in subclasses"""
        return 0 
    
    def screenout(self):
        pass