
from cc3d import CompuCellSetup
        

from Invitro_NTC_2DSteppables import CreateCellClusters
CompuCellSetup.register_steppable(steppable=CreateCellClusters(frequency=1))

from Invitro_NTC_2DSteppables import NeuraltubeModelCellBehaviors
CompuCellSetup.register_steppable(steppable=NeuraltubeModelCellBehaviors(frequency=1))

# from Invitro_NTC_2DSteppables import MitosisSteppableClusters
# CompuCellSetup.register_steppable(steppable=MitosisSteppableClusters(frequency=10))

from Invitro_NTC_2DSteppables import MitosisSteppable
CompuCellSetup.register_steppable(steppable=MitosisSteppable(frequency=10))


# from Invitro_NTC_2DSteppables import ConstraintInitializerSteppable
# CompuCellSetup.register_steppable(steppable=ConstraintInitializerSteppable(frequency=1))

# from Invitro_NTC_2DSteppables import MitosisSteppable
# CompuCellSetup.register_steppable(steppable=MitosisSteppable(frequency=10))


CompuCellSetup.run()
