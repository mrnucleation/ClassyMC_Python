

class Simulation_Classy:
    def __init__(self, simulation_name, simulation_type, simulation_parameters. **kwargs):
        """
        Initializes the Simulation_Classy with a name, type, and parameters.
        """
        self.simulation_name = simulation_name
        self.simulation_type = simulation_type
        self.simulation_parameters = simulation_parameters
        
        self.constraints = [constr for constr in kwargs.get('constraints', [])]
        self.moves = [move for move in kwargs.get('moves', [])]
        self.traj_array = kwargs.get('traj_array', None)
        self.analysis_array = kwargs.get('analysis_array', None)
        
        assert len(self.moves) > 0, "At least one move must be defined for the simulation."
        
    def prologue(self):
        for constraint in self.constraints:
            constraint.prologue()
        for move in self.moves:
            move.prologue()
        if self.traj_array:
            for traj in self.traj_array:
                traj.prologue()
        if self.analysis_array:
            for analysis in self.analysis_array:
                analysis.prologue()
                
    def epilogue(self):
        for constraint in self.constraints:
            constraint.epilogue()
        for move in self.moves:
            move.epilogue()
        if self.traj_array:
            for traj in self.traj_array:
                traj.epilogue()
        if self.analysis_array:
            for analysis in self.analysis_array:
                analysis.epilogue()
                
   

    def run_simulation(self):
        print(f"Running {self.simulation_type} simulation: {self.simulation_name}")
        # Here you would add the logic to run the simulation based on the parameters
        # For now, we will just simulate a simple output
        return f"Simulation {self.simulation_name} completed successfully."