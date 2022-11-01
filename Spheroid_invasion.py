
from cc3d import CompuCellSetup


from Spheroid_invasionSteppables import CalculateARandFrontLeaders
CompuCellSetup.register_steppable(steppable=CalculateARandFrontLeaders(frequency=1))

from Spheroid_invasionSteppables import Calculate_P
CompuCellSetup.register_steppable(steppable=Calculate_P(frequency=1))


from Spheroid_invasionSteppables import Spheroid_invasionSteppable
CompuCellSetup.register_steppable(steppable=Spheroid_invasionSteppable(frequency=1))



CompuCellSetup.run()
