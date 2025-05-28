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
        """
        constructor exception object

        :param message: takes in the error message
        :type message: str
        """
        self.message = message

    def __str__(self):
        return self.message

class UserInputException(Exception):
    """
    Exception raised for invalid user inputs from the GUI.
    """

    def __init__(self, message: str, inputs: dict = None):
        """
        Initialize the exception with a message and optional input dictionary.

        :param message: Description of the validation error.
        :param inputs: Dictionary of user inputs (optional).
        """
        super().__init__(message)
        self.message = message
        self.user_inputs = inputs or {}

    def __str__(self):
        return f"UserInputException: {self.message}"

