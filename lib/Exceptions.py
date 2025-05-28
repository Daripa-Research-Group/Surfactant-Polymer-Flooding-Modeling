class OutOfRangeError(Exception):
    """Exception raised for errors in the input if it is out of range."""

    def __init__(self, value, message="Value is out of the allowed range."):
        self.value = value
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f"{self.value} -> {self.message}"

class SimulationCalcInputException(Exception):
    """
    Exception Handling for required inputs within various stages of the simulation
    """
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return self.message

class UserInputException(Exception):
    """
    Exception Handling for user inputs retrieved from GUI
    """

