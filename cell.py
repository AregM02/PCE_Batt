from models.composite import BatteryModel

class Cell:
    def __init__(self, initial_soc=0, temperature=25):
        self.model = BatteryModel()
        self.soc = initial_soc
        self.T = temperature

    def solve(self, current, time):
        """Invokes the solver method of the chosen model."""

        # TODO: add soc calculation from the initial value here
        
        soc = self.soc
        T = self.T
        voltage_mean, voltage_sigma = self.model.solve(current, time, soc, T)
        return voltage_mean, voltage_sigma