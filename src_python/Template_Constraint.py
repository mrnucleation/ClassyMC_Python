class Constraint():
    """
    A Python equivalent of the Fortran `constraint` type.
    """

    def __init__(self, box_id: int):
        """
        Constructor(self, boxID)
        """
        super().__init__()  # initialize base ClassyClass
        self.box_id = box_id

    def check_initial_constraint(self, trial_box) -> bool:
        """
        CheckInitialConstraint(self, trialBox, accept)
        Returns True if the constraint is initially satisfied.
        """
        return True

    def diff_check(self, trial_box, disp) -> bool:
        """
        DiffCheck(self, trialBox, disp, accept)
        Returns True if the differential constraint passes.
        """
        return True

    def post_energy(self,
                    trial_box: "SimBox",
                    disp: list["Perturbation"],
                    e_diff: float) -> bool:
        """
        PostEnergy(self, trialBox, disp, E_Diff, accept)
        Returns True if the post‐energy check passes.
        """
        # TODO: insert real logic here
        return True

    def process_io(self, line) -> int:
        """
        ProcessIO(self, line, lineStat)
        Parses a line of input and returns a status code.
        """
        # TODO: insert real logic here
        return 0


class ConstrainArray:
    """
    A container for multiple Constraint instances.
    """

    def __init__(self):
        # equivalent to `class(constraint), allocatable :: method`
        self.methods: list[Constraint] = []

    def add(self, constraint) -> None:
        """
        Add a new Constraint to the array.
        """
        self.methods.append(constraint)

    def __getitem__(self, idx: int) -> Constraint:
        return self.methods[idx]

    def __len__(self) -> int:
        return len(self.methods)
