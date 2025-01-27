"""
This python script will contain the class definition for surfactants

"""

class Surfactant:
    def __init__(self, name, initial_concentration, IFT_conc_equ, vec_concentration = None):
        """
        Creates instance of Surfactant class

        :param name: Name of the surfactant
        :type name: enum 'SurfactantList'

        :param concentration: Initial concentration in wppm of surfactant (scalar quantity)
        :type concentration: float

        :param IFT_conc_equ: expression that relates surfactant concentration to interfacial tension b/t oil and water
        :type IFT_conc_equ: lambda

        :param vec_concentration: vector representation of surfactant concentration in resevoir
        :type vec_concentration: np.array, None
        """
        self.name = name
        self.concentration = initial_concentration
        self.vec_concentration = vec_concentration
        self.IFT_conc_equ = IFT_conc_equ


