class OutOfRangeError(Exception):
    """Exception raised for errors in the input if it is out of range."""

    def __init__(self, value, message="Value is out of the allowed range."):
        """
        Args:
        -----
            value (float): value that the user provided in GUI
            message (str): Error message
        """
        self._value = value
        self._message = message
        super().__init__(self.message)

    def __str__(self):
        """
        Returns: (str)
        --------------
            Returns the error message
        """
        return f"{self._value} -> {self._message}"


class SimulationCalcInputException(Exception):
    """
    Exception Handling for required inputs within various stages of the simulation
    """

    def __init__(self, message):
        """
        constructor exception object
        
        Args:
        -----
            message (str): takes in the error message
        """
        self._message = message

    def __str__(self):
        """
        Returns: (str)
        --------------
            Returns the error message
        """
        return self._message




class UserInputException(Exception):
    """
    Exception raised for invalid user inputs from the GUI.
    """

    def __init__(self, message: str, inputs: dict | None = None):
        """
        Initialize the exception with a message and optional input dictionary.

        Args:
        -----
            message (str): Description of the validation error.
            inputs (dict, None): Dictionary of user inputs (optional).
        """
        super().__init__(message)
        self._message = message
        self._user_inputs = inputs or {}

    def __str__(self):
        """
        Returns: (str)
        --------------
            Returns the error message
        """
        return f"UserInputException: {self._message}"
