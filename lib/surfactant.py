"""
This python script will contain the class definition for surfactants

"""

from Exceptions import OutOfRangeError

class Surfactant:
    def __init__(self, name, concentration, IFT_conc_equ):
        """
        Creates instance of Surfactant class

        :param name: Name of the surfactant
        :type name: str

        :param concentration: Concentration in wppm of surfactant
        :type concentration: float

        :param IFT_conc_equ: expression that relates surfactant concentration to interfacial tension b/t oil and water
        :type IFT_conc_equ: str
        """
        self.name = name
        self.concentration = concentration
        self.IFT_conc_equ = IFT_conc_equ


