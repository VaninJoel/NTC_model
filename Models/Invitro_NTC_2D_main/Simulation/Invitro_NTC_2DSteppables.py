# bATCH
import sys
import numpy as np
from os.path import join
import os
# Import project libraries and classes
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))))

# Import project libraries and classes
sys.path.append(os.path.dirname(__file__))
from BatchRun import BatchRunLib

# Import toolkit
sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from BatchInputs import *

##########################################################

from cc3d.core.PySteppables import *
import math
import random
import parameters as p

#################################################
## Define Constants for referencing Cell types ##
MEDIUM = 0  ##
ENDODERM = 1  ##
NOTOCHORD = 2  ##
MESODERM_LEFT = 3  ##
MESODERM_RIGHT = 4  ##
MHP_BASAL = 5  ##
MHP_LATERAL = 6  ##
MHP_APICAL = 7  ##
NE_BASAL_LEFT = 8  ##
NE_LATERAL_LEFT = 9  ##
NE_APICAL_LEFT = 10  ##
NE_BASAL_RIGHT = 11  ##
NE_LATERAL_RIGHT = 12  ##
NE_APICAL_RIGHT = 13  ##
ECTODERM_LEFT = 14  ##
ECTODERM_RIGHT = 15  ##
WALL = 16  ##
ECM = 17
fluid = 18  ##
#################################################
## Constants to translate cell ID to cell name ##
type_list = [  ##
    "Medium",  ##
    "Endoderm",  ##
    "Notochord",  ##
    "Mesoderm_Left",  ##
    "Mesoderm_Right",  ##
    "MHP_Basal",
    "MHP_Lateral",
    "MHP_Apical",
    "NE_Basal_Left",
    "NE_Lateral_Left",
    "NE_Apical_Left",
    "NE_Basal_Right",
    "NE_Lateral_Right",
    "NE_Apical_Right",
    "Ectoderm_Left",
    "Ectoderm_Right",
    "Wall",
    "ECM",
    "fluid"]

#################################################
n_cell_types = 19  # includes medium
mitosis_count_list = [0 for i in range(n_cell_types)]

# global variables to store the plot window and distance data
plot_window = None
distance_data = []


class CreateCellClusters(SteppableBasePy):
    def __init__(self, frequency=1):

        SteppableBasePy.__init__(self, frequency)
        import parameters
        BatchRunLib.apply_external_multipliers(__name__, parameters)

    def start(self):
        self.build_wall(self.WALL)

        for cell in self.cell_list_by_type(self.MHP_APICAL):
            for cell2 in self.cell_list_by_type(self.MHP_LATERAL):
                if cell.xCOM == cell2.xCOM:
                    self.reassign_cluster_id(cell2, cell.clusterId)
                    break

            for cell2 in self.cell_list_by_type(self.MHP_BASAL):
                if cell.xCOM == cell2.xCOM:
                    self.reassign_cluster_id(cell2, cell.clusterId)
                    break

        for cell in self.cell_list_by_type(self.NP_APICAL_L):
            for cell2 in self.cell_list_by_type(self.NP_LATERAL_L):
                if cell.xCOM == cell2.xCOM:
                    self.reassign_cluster_id(cell2, cell.clusterId)
                    break

            for cell2 in self.cell_list_by_type(self.NP_BASAL_L):
                if cell.xCOM == cell2.xCOM:
                    self.reassign_cluster_id(cell2, cell.clusterId)
                    break

        for cell in self.cell_list_by_type(self.NP_APICAL_R):
            for cell2 in self.cell_list_by_type(self.NP_LATERAL_R):
                if cell.xCOM == cell2.xCOM:
                    self.reassign_cluster_id(cell2, cell.clusterId)
                    break

            for cell2 in self.cell_list_by_type(self.NP_BASAL_R):
                if cell.xCOM == cell2.xCOM:
                    self.reassign_cluster_id(cell2, cell.clusterId)
                    break

        for cell in self.cell_list_by_type(self.MHP_LATERAL, self.NP_LATERAL_L, self.NP_LATERAL_R):
            for cell2 in self.cell_list_by_type(self.MHP_APICAL, self.MHP_BASAL, self.NP_BASAL_L, self.NP_BASAL_R,
                                                self.NP_APICAL_L, self.NP_APICAL_R):
                if cell.clusterId == cell2.clusterId:
                    # access/modification of a dictionary attached to cell - make sure to declare in main script that
                    # you will use such attribute
                    # access/modification of a dictionary attached to cell - make sure to declare in main script that
                    # you will use such attribute
                    cell.dict['contact'] = 0

                    link = self.new_fpp_internal_link(cell, cell2, 40, 3, 40)
                    # To select which types of links to remove, set any of the following keyword arguments to True
                    #   Remove all links          : links
                    #   Remove all internal links : internal_links
                    #   Remove all anchors        : anchors
                    # If none of these are specified, then all links of any type are removed
                    # self.remove_all_cell_fpp_links(cell)

    def step(self, mcs):
        if not mcs % 1000:
            self.save_sim_as_numpy(mcs=mcs)
        if mcs == 10:

            # for cell in self.cell_list_by_type(self.MHP_APICAL):            
            # self.remove_all_cell_fpp_links(cell)
            # for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):

            # if neighbor and neighbor.type==self.MHP_APICAL:
            # lambda_distance = 30
            # target_distance = 1
            # max_distance = 20
            # link = self.new_fpp_link(cell, neighbor, lambda_distance, target_distance, max_distance)

            for cell in self.cell_list_by_type(self.MHP_LATERAL):
                for cell2 in self.cell_list_by_type(self.MHP_APICAL, self.MHP_BASAL):
                    if cell.clusterId == cell2.clusterId:
                        link = self.new_fpp_internal_link(cell, cell2, 40, 3, 40)

            for cell in self.cell_list_by_type(self.MHP_LATERAL, self.NP_LATERAL_L, self.NP_LATERAL_R):
                for cell2 in self.cell_list_by_type(self.MHP_APICAL, self.MHP_BASAL, self.NP_BASAL_L, self.NP_BASAL_R,
                                                    self.NP_APICAL_L, self.NP_APICAL_R):
                    if cell.clusterId == cell2.clusterId:
                        link = self.new_fpp_internal_link(cell, cell2, 40, 3, 40)

            for cell in self.cell_list_by_type(self.NP_APICAL_L, self.NP_APICAL_R):
                for cell2 in self.cell_list_by_type(self.NP_BASAL_L, self.NP_BASAL_R):
                    if cell.clusterId == cell2.clusterId:
                        link = self.new_fpp_internal_link(cell, cell2, 40, 12, 40)

            for cell in self.cell_list_by_type(self.MHP_APICAL):
                for cell2 in self.cell_list_by_type(self.MHP_BASAL):
                    if cell.clusterId == cell2.clusterId:
                        link = self.new_fpp_internal_link(cell, cell2, 40, 12, 40)

            for cell in self.cell_list_by_type(self.NP_APICAL_R):
                for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                    if neighbor:
                        if neighbor.type in [NE_APICAL_RIGHT]:
                            apical_link = self.new_fpp_link(cell, neighbor, 20, 5, 20)
                            # apical_link = self.get_fpp_link_by_cells(apical_parent, neighbor)
                            apical_link.setMaxNumberOfJunctions(1)

    def save_sim_as_numpy(self, mcs=12000):

        mega_image_cheat_sheet = "Channel list: \n " \
                                 "- Channel 0 = cell id \n" \
                                 "- Channel 1 = cluster id \n" \
                                 "- Channel 2 = type \n" \
                                 "- Channel 3 = BMP4 \n" \
                                 "- Channel 4 = SHH \n" \
                                 "- Channel 5 = COM x \n" \
                                 "- Channel 6 = COM y \n" \
                                 "- Channel 7 = pressure"

        number_channels = len(mega_image_cheat_sheet.split("\n")) - 1

        mega_image = np.zeros([self.dim.x, self.dim.y, number_channels], dtype=float)
        BMP4 = self.field.BMP4
        SHH = self.field.SHH

        for x, y, z in self.every_pixel():
            cell = self.cell_field[x, y, z]
            if cell:
                mega_image[x, y, 0] = cell.id
                mega_image[x, y, 1] = cell.clusterId
                mega_image[x, y, 2] = cell.type
                mega_image[x, y, 5] = cell.xCOM
                mega_image[x, y, 6] = cell.yCOM
                mega_image[x, y, 7] = cell.pressure
            else:
                mega_image[x, y, 0] = np.nan
                mega_image[x, y, 1] = np.nan
                mega_image[x, y, 2] = np.nan
                mega_image[x, y, 5] = np.nan
                mega_image[x, y, 6] = np.nan
                mega_image[x, y, 7] = np.nan
            mega_image[x, y, 3] = BMP4[x, y, 0]
            mega_image[x, y, 4] = SHH[x, y, 0]


        with open(join(self.output_dir, "image_cheat_sheet.txt"), "w+") as fout:
            fout.write(mega_image_cheat_sheet+"\n")

        np.save(join(self.output_dir, f"sim_field_{mcs}.npy"), mega_image)

    def finish(self):
        self.save_sim_as_numpy()
        return

    def on_stop(self):
        self.finish()
        return


class NeuraltubeModelCellBehaviors(SteppableBasePy):
    def __init__(self, frequency):
        SteppableBasePy.__init__(self, frequency)

        # self.track_cell_level_scalar_attribute(field_name='cell_id', attribute_name='id')
        # self.track_cell_level_scalar_attribute(field_name='target_volume', attribute_name='tvolume')
        # self.track_cell_level_scalar_attribute(field_name='target_surface', attribute_name='tsurface')
        # self.track_cell_level_scalar_attribute(field_name='lambda_volume', attribute_name='lvolume')
        # self.track_cell_level_scalar_attribute(field_name='lambda_surface', attribute_name='lsurface')
        # self.track_cell_level_scalar_attribute(field_name='x_force_push', attribute_name="X_PUSH_DIRECTION")
        # self.track_cell_level_scalar_attribute(field_name='cell_pressure', attribute_name='pressure')
        # self.track_cell_level_scalar_attribute(field_name='cell_cycle', attribute_name="CELL_CYCLE")

        ## general cell parameters ##
        self.standard_target_volume = 49
        self.standard_lambda_volume = 2.0
        self.target_surface_modifier = 0
        self.lambda_surface_modifier = 0.3
        self.standard_target_surface = 25  # self.target_surface_modifier * (math.pi * (2*math.sqrt(self.standard_target_volume / math.pi)))
        self.standard_lambda_surface = 2.0  # self.standard_lambda_volume * self.lambda_surface_modifier
        self.time_to_stop_mesodermal_growth = 4000  # 12000 #mcs default value that stops cell growth
        self.time_to_stop_ectodermal_growth = 4000
        self.cell_cycle_time = 4000
        # self.mesoderm_growth_rate = 0.01 #was 0.015 #was 0.0030
        # self.mesoderm_max_growth_rate = 1.0
        self.ectoderm_growth_rate = 0.0010
        self.somite_growth_rate = 0.0
        # self.ECM_growth_rate = 0.01

        self.mesodermal_differentiation_probability = 5 / 1000
        # self.apoptosis_chance = 0.001

        ## gene parameters ## Include genes from Kazburn et al, 2021
        # ATRA --> Time based parameter, starts Low, ends high
        # FGF --> Time based variable, starts high, ends low.
        # NOG --> Inhibits BMP
        # N/E Cadherin switch --> for fusion, regulated by Snail2
        # Snail2 is regulated by WNT
        self.gene_status_dict = {
            "BMP4": 1.0,
            "SHH": 1.0
        }

        ## BMP4 ##
        # self.BMP4_secretion_constant = 0.020
        self.AC50_BMP4 = 50
        self.BMP4_gradient_low = 0.1
        self.BMP4_gradient_high = 0.8

        ## SHH ##
        # self.SHH_secretion_constant = 0.010
        self.SHH_gradient_low = 0.1
        self.SHH_gradient_high = 0.6

        ## Toxicological perturbation parameters ##
        self.dose_toxic_substance = False

        ##Apical constriction parameters##
        self.enable_apical_constriction_behavior = True
        self.enable_mhp_formation = True
        # self.ne_apical_target_volume = p.ne_apical_target_volume#21 # The target volume of apical cell under apical constriction
        # self.ne_basal_target_volume = 77 # The target volume of basal cell under apical constriction
        # self.mhp_apical_target_volume = 21
        # self.mhp_basal_target_volume = 77

        self.mhp_constrict_frequency = 1
        self.apical_constrict_frequency = 1  # The number of MCS for one step in apical constriction
        self.delta_target_volume = 0.04  # the amount the target volume changes per step frequency
        self.delta_lambda_volume = 0.008  # makes the small cells not dissapear # the amount the lambda volume changes per step frequency
        self.delta_target_volume_apical = 0.04
        self.apical_constriction_elongation = 18  # the target distance of fpp links between the basal and apical cells under apical constriction
        self.delta_elongation = 0.01  # the step size of elongation difference
        self.standard_elongation = 12  # the standard target distance of fpp links between the basal and apical cells

        # self.freeze_apical_constriction_time = 5500

        ## Apical constriction FPP links parameters
        self.initial_lambda_distance = 20.0  # The initial lambda distance of FPP links
        self.initial_target_distance = 5.0  # The initial target distance of FPP links
        self.initial_max_distance = 20.0  # The initial max distance of FPP links
        self.apical_constrict_lambda_distance = 40.0  # The lambda distance of FPP links under Apical constriction
        self.apical_constrict_target_distance = 0.0  # The lambda distance of FPP links under Apical constriction
        self.basal_constrict_target_distance = 9.0
        self.apical_constrict_max_distance = 15.0  # The max distance of FPP links under Apical constriction
        self.delta_distance = 0.005

        ## Cell layout parameters
        self.use_ECM = True

    def set_standard_target_lambda_volume(self, cell):
        if cell.type in [MESODERM_LEFT, MESODERM_RIGHT, ECTODERM_LEFT, ECTODERM_RIGHT]:
            cell.targetVolume = random.uniform(self.standard_target_volume - 10, self.standard_target_volume + 10)
            cell.lambdaVolume = self.standard_lambda_volume
        if cell.type in [ECM]:
            cell.targetVolume = 25
            cell.lambdaVolume = 0.25
        else:
            cell.targetVolume = self.standard_target_volume
            cell.lambdaVolume = self.standard_lambda_volume
        # cell.targetSurface = self.standard_target_surface
        # cell.lambdaSurface = self.standard_target_surface

    def hill_function(self, _fieldValue, _maxRate, _AC50, _n=1):
        """
        field_value = The concentration of the ligand
        max_rate = the maximum growth rate
        AC50 = estimates the concentration at which a chemical produces the half-maximal response along a sigmoidal curve
        n = the Hill coefficient 
        """
        return _maxRate * ((_fieldValue ** _n) / (_AC50 ** _n + _fieldValue ** _n))

    # fieldSHH.get(pt), _maxRate=maxSecreteRate, _AC50=20, _n=1
    def make_mhp_fpp_links(self, cell):
        if cell.type in [MHP_APICAL]:
            for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                if neighbor and neighbor.type in [MHP_APICAL]:
                    fpp_link = self.get_fpp_link_by_cells(cell, neighbor)
                    fpp_link.setLambdaDistance(self.apical_constrict_lambda_distance)
                    fpp_link.setTargetDistance(self.apical_constrict_target_distance)
                    fpp_link.setMaxDistance(self.apical_constrict_max_distance)

    def mhp_apical_constriction(self, cell, apical_tv, basal_tv, delta_tv, delta_lv, constrict_frequency):
        if self.mcs % constrict_frequency == 0:
            if cell.type in [MHP_APICAL]:
                if cell.targetVolume >= apical_tv:  # slightly modify target volume
                    cell.targetVolume -= delta_tv
                    cell.lambdaVolume += delta_lv

                if self.mcs > 10:
                    cluster_cell_list = self.get_cluster_cells(cell.clusterId)
                    basal_cell = cluster_cell_list[0]
                    lateral_cell = cluster_cell_list[1]
                    apical_cell = cluster_cell_list[2]

                    fpp_link = self.get_fpp_internal_link_by_cells(basal_cell, apical_cell)
                    target_distance = fpp_link.getTargetDistance()
                    if target_distance < self.apical_constriction_elongation:
                        fpp_link.setTargetDistance(target_distance + self.delta_elongation)

                    # Adjust FPP links between apical NE cells
                    for neighbor, common_surface_area in self.get_cell_neighbor_data_list(apical_cell):
                        if neighbor and neighbor.type in [MHP_APICAL]:
                            fpp_link = self.get_fpp_link_by_cells(apical_cell, neighbor)
                            if fpp_link:  # possibly at a check to see whether neighbor is True/False compartment
                                target_distance = fpp_link.getTargetDistance()
                                if target_distance > self.apical_constrict_target_distance:
                                    fpp_link.setTargetDistance(target_distance - self.delta_distance)

                    # Adjust FPP links between basal NE cells
                    for neighbor, common_surface_area in self.get_cell_neighbor_data_list(basal_cell):
                        if neighbor and neighbor.type in [MHP_BASAL]:
                            fpp_link = self.get_fpp_link_by_cells(basal_cell, neighbor)
                            if fpp_link:
                                target_distance = fpp_link.getTargetDistance()
                                if target_distance < self.basal_constrict_target_distance:
                                    fpp_link.setTargetDistance(target_distance + self.delta_distance)

            if cell.type in [MHP_BASAL]:
                if cell.targetVolume <= basal_tv:  # slightly modify target volume
                    cell.targetVolume += delta_tv
                    cell.lambdaVolume += delta_lv

                elif cell.targetVolume >= basal_tv:
                    for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                        if neighbor:
                            if neighbor.type == NOTOCHORD:
                                link = self.get_fpp_link_by_cells(cell, neighbor)
                                if not link:
                                    apical_link = self.new_fpp_link(cell, neighbor, 20, 6, 20)
                                    apical_link = self.get_fpp_link_by_cells(cell, neighbor)
                                    apical_link.setMaxNumberOfJunctions(1)

    def neuralectoderm_apical_constriction(self, cell, apical_tv, basal_tv, delta_tv, delta_lv, constrict_frequency):
        if cell.type in [NE_LATERAL_RIGHT]:  # use the lateral_cell to determine the behaviour of cells in the cluster
            constrict_flag = cell.dict["ENABLE_CONSTRICT"]
            if self.mcs % constrict_frequency == 0:
                # get cells from the cluster
                cluster_cell_list = self.get_cluster_cells(cell.clusterId)

                basal_cell = cluster_cell_list[0]
                lateral_cell = cluster_cell_list[1]
                apical_cell = cluster_cell_list[2]

                if constrict_flag:  # start apical constriction behavior --> change cell compartment sizes
                    # Change Apical_cell target volume                
                    if apical_cell.targetVolume > apical_tv:  # slightly modify target volume
                        apical_cell.targetVolume -= self.delta_target_volume_apical
                        apical_cell.lambdaVolume += delta_lv

                    # Change basal cell
                    if basal_cell.targetVolume < basal_tv:  # slightly modify target volume
                        basal_cell.targetVolume += delta_tv
                        basal_cell.lambdaVolume += delta_lv

                    # Adjust FPP links between apical and basal cell
                    fpp_link = self.get_fpp_internal_link_by_cells(basal_cell, apical_cell)
                    target_distance = fpp_link.getTargetDistance()
                    if target_distance < self.apical_constriction_elongation:
                        fpp_link.setTargetDistance(target_distance + self.delta_elongation)

                    # Adjust FPP links between apical NE cells
                    for neighbor, common_surface_area in self.get_cell_neighbor_data_list(apical_cell):
                        if neighbor and neighbor.type in [NE_APICAL_RIGHT]:
                            fpp_link = self.get_fpp_link_by_cells(apical_cell, neighbor)
                            if fpp_link:  # possibly at a check to see whether neighbor is True/False compartment
                                target_distance = fpp_link.getTargetDistance()
                                if target_distance > self.apical_constrict_target_distance:
                                    fpp_link.setTargetDistance(target_distance - self.delta_distance)

                    # Adjust FPP links between basal NE cells
                    for neighbor, common_surface_area in self.get_cell_neighbor_data_list(basal_cell):
                        if neighbor and neighbor.type in [NE_BASAL_RIGHT]:
                            fpp_link = self.get_fpp_link_by_cells(basal_cell, neighbor)
                            if fpp_link:
                                target_distance = fpp_link.getTargetDistance()
                                if target_distance < self.basal_constrict_target_distance:
                                    fpp_link.setTargetDistance(target_distance + self.delta_distance)

                if not constrict_flag:  # go back to normal cell compartment sizes, same steps as constriction rate
                    # Change Apical_cell target volume                
                    if apical_cell.targetVolume < self.standard_target_volume:  # slightly modify target volume
                        apical_cell.targetVolume += self.delta_target_volume_apical
                        apical_cell.lambdaVolume -= delta_lv

                    # Change basal cell
                    if basal_cell.targetVolume > self.standard_target_volume:  # slightly modify target volume
                        basal_cell.targetVolume -= delta_tv
                        basal_cell.lambdaVolume -= delta_lv

                    # Adjust FPP links between apical and basal cell
                    fpp_link = self.get_fpp_internal_link_by_cells(basal_cell, apical_cell)
                    target_distance = fpp_link.getTargetDistance()
                    if target_distance > self.standard_elongation:
                        fpp_link.setTargetDistance(target_distance - self.delta_elongation)

                    # Adjust FPP links between apical NE cells
                    for neighbor, common_surface_area in self.get_cell_neighbor_data_list(basal_cell):
                        if neighbor and neighbor.type in [NE_APICAL_RIGHT]:
                            fpp_link = self.get_fpp_link_by_cells(basal_cell, neighbor)
                            if fpp_link:
                                target_distance = fpp_link.getTargetDistance()
                                if target_distance < self.initial_target_distance:
                                    fpp_link.setTargetDistance(target_distance + self.delta_distance)

                    # Adjust FPP links between basal NE cells
                    for neighbor, common_surface_area in self.get_cell_neighbor_data_list(basal_cell):
                        if neighbor and neighbor.type in [NE_BASAL_RIGHT]:
                            fpp_link = self.get_fpp_link_by_cells(basal_cell, neighbor)
                            if fpp_link:
                                target_distance = fpp_link.getTargetDistance()
                                if target_distance > self.initial_target_distance:
                                    fpp_link.setTargetDistance(target_distance - self.delta_distance)

    def set_behaviour_flag_apical_constriction_BMP4(self, cell, field, field_low, field_high):

        if cell.type in [NE_LATERAL_RIGHT]:  # [NE_APICAL_LEFT, NE_APICAL_RIGHT, NE_BASAL_LEFT, NE_BASAL_RIGHT]:
            field_concentration_at_cell = field[int(cell.xCOM), int(cell.yCOM), int(cell.zCOM)]

            if field_low <= field_concentration_at_cell <= field_high:
                cell.dict["ENABLE_CONSTRICT"] = True
            if not field_low <= field_concentration_at_cell <= field_high:
                cell.dict["ENABLE_CONSTRICT"] = False

    def set_behaviour_flag_apical_constriction_SHH(self, cell, field, field_high):

        if cell.type in [NE_LATERAL_RIGHT]:  # [NE_APICAL_LEFT, NE_APICAL_RIGHT, NE_BASAL_LEFT, NE_BASAL_RIGHT]:
            field_concentration_at_cell = field[int(cell.xCOM), int(cell.yCOM), int(cell.zCOM)]
            if not field_concentration_at_cell <= field_high:
                cell.dict["ENABLE_CONSTRICT"] = False

    def calc_apical_constriction_fpplink_change(self, cell, fpp_low, fpp_high):
        if cell.type in [NE_LATERAL_RIGHT]:  # [NE_APICAL_LEFT, NE_APICAL_RIGHT]:
            if cell.dict["ENABLE_CONSTRICT"]:
                if fpp_low <= cell.dict["ENABLE_FPPLINK_CHANGE"] <= fpp_high:
                    cell.dict["ENABLE_FPPLINK_CHANGE"] += 1

            if not cell.dict["ENABLE_CONSTRICT"]:
                if fpp_low <= cell.dict["ENABLE_FPPLINK_CHANGE"] <= fpp_high:
                    cell.dict["ENABLE_FPPLINK_CHANGE"] -= 1

            # ensure that the FPP_link_change is not out of bounds    
            if cell.dict["ENABLE_FPPLINK_CHANGE"] > fpp_high:
                cell.dict["ENABLE_FPPLINK_CHANGE"] = fpp_high
            elif cell.dict["ENABLE_FPPLINK_CHANGE"] < fpp_low:
                cell.dict["ENABLE_FPPLINK_CHANGE"] = fpp_low

    def constrict_all_links(self, cell, constrict_none=True):
        if constrict_none:
            pass
        else:

            if cell.type in [NE_APICAL_LEFT, NE_APICAL_RIGHT]:
                for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                    if neighbor and neighbor.type in [NE_APICAL_LEFT, NE_APICAL_RIGHT]:
                        fpp_link = self.get_fpp_link_by_cells(cell, neighbor)
                        fpp_link.setLambdaDistance(self.apical_constrict_lambda_distance)
                        fpp_link.setTargetDistance(self.apical_constrict_target_distance)
                        fpp_link.setMaxDistance(self.apical_constrict_max_distance)

    def track_cell_cycle(self, cell):
        chance = random.randint(0, 3)
        cell.dict["CELL_CYCLE"] += chance

        if cell.dict["CELL_CYCLE"] >= self.cell_cycle_time:
            cell.dict["ALLOW_DIVISION"] = True

    def force_push(self, cell, mcs, push_time, type_list, x, y, z):
        if mcs == push_time:
            if cell.type in type_list:
                cell.lambdaVecX = x  # was -0.5# force component pointing along X axis - towards positive X's
                cell.lambdaVecY = y  # was -0.5 force component pointing along Y axis - towards negative Y's
                cell.lambdaVecZ = z

        # if cell in [ECTODERM_LEFT]:
        # cell.lambdaVecX = random.uniform(-4.0, -0.5)
        # cell.lambdaVecY = random.uniform(-1.0, 0.0)

        # if cell in [ECTODERM_RIGHT]:
        # cell.lambdaVecX = random.uniform(0.5, 4.0)
        # cell.lambdaVecY = random.uniform(-1.0, 0.0)

    def add_mesoderm_matrix_cell(self, cell):
        if cell.type in [MESODERM_LEFT, MESODERM_RIGHT]:
            # if _cell.volume>2 and _cell.targetVolume>2:
            pixelTrackerDataList = CellBoundaryPixelList(self.boundaryPixelTrackerPlugin,
                                                         cell)  # careful, not really a python list
            pixelList = []
            for boundaryPixelTrackerData in pixelTrackerDataList:
                point3dpixel = boundaryPixelTrackerData.pixel
                np_pixel = self.point_3d_to_numpy(point3dpixel)
                l_pixel = np_pixel.tolist()
                pixelList.append(l_pixel)
                # print(l_pixel)
            pt = random.choice(pixelList)
            # print("Random choosen pixel " + str(pt))
            # f pt:

            # newMatrix = self.cell_field[pt[0], pt[1], pt[2]]
            x = pt[0]
            y = pt[1]
            size = 2
            # chance = random.random()
            # if chance < random.random():
            ncell = self.new_cell(self.ECM)
            # size of cell will be SIZExSIZEx1
            self.cell_field[x:x + size - 1, y:y + size - 1, 0] = ncell

    def condense_somite_cells(self):
        # iterating over cells of type 1
        # list of  cell types (capitalized)
        # mesoderm_cell_list_L = []
        if [x for x in self.cell_list_by_type(MESODERM_RIGHT)] == []:
            for condense_cell in self.cell_list_by_type(MESODERM_LEFT):
                #       mesoderm_cell_list_L.append(cell)
                # mesoderm_cell_list = [cell for cell in self.cell_list_by_type(self.MESODERM_LEFT, self.MESODERM_RIGHT)]
                # random.shuffle(mesoderm_cell_list_L)

                # condense_cell = mesoderm_cell_list_L[0]

                for neighbor, common_surface_area in self.get_cell_neighbor_data_list(condense_cell):

                    if neighbor:
                        if neighbor.type in [MESODERM_LEFT]:
                            vec = self.invariant_distance_between_cells(condense_cell, neighbor)
                            # print(vec)
                            if vec <= 8.0:
                                link = self.new_fpp_link(
                                    condense_cell,
                                    neighbor,
                                    self.apical_constrict_lambda_distance,
                                    0.0,
                                    5)

    def condense_notochord_cells(self):
        # iterating over cells of type 1
        # list of  cell types (capitalized)
        notochord_cell_list = []
        for cell in self.cell_list_by_type(NOTOCHORD):
            notochord_cell_list.append(cell)
        # mesoderm_cell_list = [cell for cell in self.cell_list_by_type(self.MESODERM_LEFT, self.MESODERM_RIGHT)]
        random.shuffle(notochord_cell_list)

        condense_cell = notochord_cell_list[0]

        for neighbor, common_surface_area in self.get_cell_neighbor_data_list(condense_cell):

            if neighbor:
                if neighbor.type in [NOTOCHORD]:
                    link = self.new_fpp_link(
                        condense_cell,
                        neighbor,
                        self.apical_constrict_lambda_distance,
                        -1.0,
                        8)

        # mesoderm_cell_list_R = []
        # for cell in self.cell_list_by_type(MESODERM_RIGHT):
        # mesoderm_cell_list_R.append(cell)
        # #mesoderm_cell_list = [cell for cell in self.cell_list_by_type(self.MESODERM_LEFT, self.MESODERM_RIGHT)]
        # random.shuffle(mesoderm_cell_list_R)

        # condense_cell = mesoderm_cell_list_R[0]
        # for neighbor, common_surface_area in self.get_cell_neighbor_data_list(condense_cell):
        # if neighbor:
        # if neighbor.type in [MESODERM_RIGHT]:
        # link = self.new_fpp_link(
        # condense_cell,
        # neighbor,
        # self.apical_constrict_lambda_distance,
        # -1.0,
        # 10)

    def differentiate_mesoderm_to_somite(self, cell):
        prob_to_diff = self.mesodermal_differentiation_probability
        if cell.type in [MESODERM_RIGHT]:
            # for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
            #    if neighbor.type in [MESODERM_LEFT]:
            #        prob_to_diff = prob_to_diff + 0.01

            if random.random() < prob_to_diff:
                cell.type = MESODERM_LEFT

    def freeze_apical_constriction_status(self):
        self.enable_apical_constriction_behavior = False
        if mcs < self.freeze_apical_constriction_time:
            return True

        if mcs > self.freeze_apical_constriction_time:
            return False

    def get_outermost_cells(self):
        # Get the neural ectoderm cell group
        # ne_group = self.getCellIdsByType("np_lateral_R")
        # iterating over cells of type 1
        # list of  cell types (capitalized)
        ne_group = [cell for cell in self.cell_list_by_type(self.NP_LATERAL_R)]

        # Get the leftmost and rightmost cells in the group
        leftmost_cell = None
        rightmost_cell = None
        for cell in ne_group:
            # cell = self.cellList[cell_id]
            if leftmost_cell is None or cell.xCOM < leftmost_cell.xCOM:
                leftmost_cell = cell
            if rightmost_cell is None or cell.xCOM > rightmost_cell.xCOM:
                rightmost_cell = cell

        for cell in ne_group:
            if cell.id == leftmost_cell.id:
                cell.dict["OUTERMOST"] = "LEFT"
            elif cell.id == rightmost_cell.id:
                cell.dict["OUTERMOST"] = "RIGHT"
            else:
                cell.dict["OUTERMOST"] = False

    def calculate_outermost_cell_distance(self):

        leftmost_cell_l = [cell for cell in self.cell_list_by_type(self.NP_LATERAL_R) if
                           cell.dict["OUTERMOST"] == "LEFT"]
        rightmost_cell_l = [cell for cell in self.cell_list_by_type(self.NP_LATERAL_R) if
                            cell.dict["OUTERMOST"] == "RIGHT"]

        leftmost_cell = leftmost_cell_l[0]
        rightmost_cell = rightmost_cell_l[0]

        # Compute the Euclidean distance between the leftmost and rightmost cells
        dx = rightmost_cell.xCOM - leftmost_cell.xCOM
        dy = rightmost_cell.yCOM - leftmost_cell.yCOM
        dz = rightmost_cell.zCOM - leftmost_cell.zCOM
        distance = sqrt(dx * dx + dy * dy + dz * dz)

        return distance

    def start(self):
        self.get_outermost_cells()
        for cell in self.cellList:
            self.set_standard_target_lambda_volume(cell)
            self.orientedGrowthPlugin.setElongationEnabled(cell,
                                                           False)  # Make sure to enable or disable elongation in ALL cells

            ## Attach dictionary value to each cell ##
            cell.dict["CONSTRAINT_WIDTH"] = 5
            cell.dict["ENABLE_CONSTRICT"] = False
            cell.dict["ENABLE_FPPLINK_CHANGE"] = random.randint(0, 200)  # Pick a random starting value
            cell.dict["CELL_CYCLE"] = random.randint(0, self.cell_cycle_time / 2)
            cell.dict["ALLOW_DIVISION"] = False
            cell.dict["X_PUSH_DIRECTION"] = 0

            cell.dict["id"] = float(cell.id)
            cell.dict["tvolume"] = cell.targetVolume
            cell.dict["tsurface"] = cell.targetSurface
            cell.dict["lvolume"] = cell.lambdaVolume
            cell.dict["lsurface"] = cell.lambdaSurface

            if cell.type in [ECTODERM_LEFT]:
                if cell.xCOM > 7:
                    cell.dict["X_PUSH_DIRECTION"] = random.uniform(-2.0, -1.25)

            if cell.type in [ECTODERM_RIGHT]:
                if cell.xCOM < 203:
                    cell.dict["X_PUSH_DIRECTION"] = random.uniform(1.25, 2.0)

            if self.use_ECM:
                self.add_mesoderm_matrix_cell(cell)

                # with open("output.txt", "w+") as f:
        #     for cell in self.cell_list:
        #             f.write("CLUSTER CELL ID=" + str(cell.clusterId) + "\t" + " type=" + type_list[cell.type] + "\n")

        # self.plot_win_distance = self.add_new_plot_window(
        #     title='Distance',
        #     x_axis_title='MonteCarlo Step (MCS)',
        #     y_axis_title='Distance',
        #     x_scale_type='linear',
        #     y_scale_type='linear',
        #     grid=True # only in 3.7.6 or higher
        # )

        # self.plot_win_distance.add_plot("Distance", style='Lines', color='red', size=5)

    def step(self, mcs):

        distance = self.calculate_outermost_cell_distance()
        # self.plot_win_distance.add_data_point("Distance", mcs, distance)
        ## Get all fields ###############################
        BMP4_field = self.field.BMP4
        SHH_field = self.field.SHH
        # TOXIC_FLUID = self.field.TOXIC_FLUID

        ## Get all field secretors ######################
        BMP4_secretor = self.get_field_secretor("BMP4")
        SHH_secretor = self.get_field_secretor("SHH")
        TOXIC_FLUID_secretor = self.get_field_secretor("TOXIC_FLUID")
        ## Add a cell dictionary with a number to it (FPP link change or something, make it a new function that is modifyable.
        ## when the number => than the set number, change the FPP links to allow apical constriction (base the change on apical_constriction TRUE)
        ## If below set number, change apical constriction to false

        # if 5000 < mcs < 5500:
        # self.condense_somite_cells()
        # if 5000 < mcs < 6000:
        # self.condense_notochord_cells()

        ## mechanistic changes that need to get biological meaning ##
        # if mcs == 7000:
        # self.freeze_apical_constriction()

        for cell in self.cellList:  ## This is the main loop through all the cells ##
            cell.dict["id"] = float(cell.id)
            cell.dict["tvolume"] = cell.targetVolume
            cell.dict["tsurface"] = cell.targetSurface
            cell.dict["lvolume"] = cell.lambdaVolume
            cell.dict["lsurface"] = cell.lambdaSurface

            self.track_cell_cycle(cell)
            if mcs > 3250:  # was 3250
                self.differentiate_mesoderm_to_somite(cell)
            # e = self.isBoundaryCell(cell, [ENDODERM])

            if mcs == 10:
                self.constrict_all_links(cell)
                # if self.enable_mhp_formation:
                # self.make_mhp_fpp_links(cell)

            # Uptake of BMP4 in each of the cells (to compensate for too high diffusion constants)
            # if mcs == 100:
            # if cell.type in [ECTODERM_LEFT, ECTODERM_RIGHT]:
            # secretor.secreteOutsideCellAtBoundary(cell, 0.1)
            # Make sure Secretion plugin is loaded
            # make sure this field is defined in one of the PDE solvers
            # you may reuse secretor for many cells. Simply define it outside the loop
            # secretor = self.get_field_secretor("FIELDNAME")
            # secretor.secreteOutsideCellAtBoundary(cell, 300)
            # tot_amount = secretor.secreteOutsideCellAtBoundaryTotalCount(cell, 300).tot_amount

            # BMP4_concentration_at_cell = BMP4_field[int(cell.xCOM), int(cell.yCOM), int(cell.zCOM)]
            if cell.type in [ECTODERM_LEFT, ECTODERM_RIGHT]:
                BMP4_secretor.secreteInsideCell(cell, p.BMP4_secretion_constant * self.gene_status_dict["BMP4"])

            if cell.type in [NOTOCHORD]:
                SHH_secretor.secreteInsideCell(cell, p.SHH_secretion_constant * self.gene_status_dict["SHH"])

            if mcs > 3000:
                if cell.type in [MHP_BASAL, MHP_APICAL, MHP_LATERAL]:
                    SHH_secretor.secreteInsideCell(cell, p.SHH_secretion_constant * self.gene_status_dict["SHH"])

            # if self.dose_toxic_substance:
            # if 5000 <= mcs <= 5010:
            # TOXIC_FLUID[75, 198, 0] = TOXIC_FLUID[0, 0, 0] + 200
            # TOXIC_FLUID[2, 198, 0] = TOXIC_FLUID[0, 0, 0] + 200
            # TOXIC_FLUID[148, 198, 0] = TOXIC_FLUID[0, 0, 0] + 200
            # TOXIC_FLUID[2, 100, 0] = TOXIC_FLUID[0, 0, 0] + 20
            # TOXIC_FLUID[148, 100, 0] = TOXIC_FLUID[0, 0, 0] + 20

            ## cell death when touching toxic substance, write in a function ##
            # if mcs > 6000:
            # if cell.type in [NE_APICAL_LEFT, NE_APICAL_RIGHT, NE_LATERAL_LEFT, NE_LATERAL_RIGHT, NE_BASAL_LEFT, NE_BASAL_RIGHT]:
            # chance = random.random()
            # fluid_concentration_at_cell = TOXIC_FLUID[int(cell.xCOM), int(cell.yCOM), int(cell.zCOM)]
            # if chance < fluid_concentration_at_cell:
            # cell.targetVolume = 0

            ## Start cell growth loop ##
            if mcs <= self.time_to_stop_mesodermal_growth:
                if cell.type in [MESODERM_RIGHT]:
                    growth_rate = p.mesoderm_growth_rate
                    mesoderm_max_growth_rate = p.mesoderm_max_growth_rate
                    BMP4_concentration_at_cell = BMP4_field[int(cell.xCOM), int(cell.yCOM), int(cell.zCOM)]
                    growth_rate += (
                        self.hill_function(BMP4_concentration_at_cell, mesoderm_max_growth_rate, self.AC50_BMP4))
                    if cell.targetVolume < 80:
                        cell.targetVolume += growth_rate

                if cell.type in [ECM]:
                    cell.targetVolume += p.ECM_growth_rate

                if cell.type in [MESODERM_LEFT]:
                    cell.targetVolume += self.somite_growth_rate

            if mcs <= self.time_to_stop_ectodermal_growth:
                if cell.type in [ECTODERM_LEFT, ECTODERM_RIGHT]:
                    growth_rate = self.ectoderm_growth_rate + random.uniform(-self.ectoderm_growth_rate,
                                                                             self.ectoderm_growth_rate)
                    cell.targetVolume += growth_rate

            # if cell.type in [NOTOCHORD]:
            # growth_rate = self.notochord_growth_rate + random.uniform (-self.notochord_growth_rate, self.notochord_growth_rate)
            # cell.targetVolume += growth_rate

            ## Median hinge point apical constriction behavior ##
            if self.enable_mhp_formation:
                self.mhp_apical_constriction(  # call MHP apical constriction function
                    cell,
                    p.apical_constrict_target_volume,
                    p.basal_constrict_target_volume,
                    self.delta_target_volume,
                    self.delta_lambda_volume,
                    self.mhp_constrict_frequency
                )

                ## Set neuroectoderm apical constriction behavior ##
            if self.enable_apical_constriction_behavior:
                if mcs >= 3000:  # start apical constriction later in time

                    self.set_behaviour_flag_apical_constriction_BMP4(
                        # set ENABLE_CONSTRICT dict flag to True / False depending on BMP4 concentration
                        cell,
                        field=BMP4_field,
                        field_low=self.BMP4_gradient_low,
                        field_high=self.BMP4_gradient_high
                    )

                    self.set_behaviour_flag_apical_constriction_SHH(
                        # set ENABLE_CONSTRICT dict flag to True / False depending on BMP4 concentration
                        cell,
                        field=SHH_field,
                        field_high=self.SHH_gradient_high
                    )

                    self.neuralectoderm_apical_constriction(  # call NE apical constriction function
                        cell,
                        p.apical_constrict_target_volume,
                        p.basal_constrict_target_volume,
                        self.delta_target_volume,
                        self.delta_lambda_volume,
                        self.apical_constrict_frequency
                    )

                # Calculate neuroectoderm apical constriction behaviour FPP link changes
                self.calc_apical_constriction_fpplink_change(
                    cell,
                    fpp_low=0,
                    fpp_high=200
                )

            ## Forced cell migration caused by physical forces ##
            self.force_push(  # call force_push function
                cell,
                mcs=mcs,  # current mcs
                push_time=500,  # push time, start at mcs 10
                type_list=[ECTODERM_RIGHT],
                x=cell.dict["X_PUSH_DIRECTION"],  # x
                y=0,  # -0.25,#-0.25,# was -0.5
                z=0)  # z

            self.force_push(  # call force_push function
                cell,
                mcs=mcs,  # current mcs
                push_time=500,  # push time, start at mcs 10
                type_list=[ECTODERM_LEFT],
                x=cell.dict["X_PUSH_DIRECTION"],  # x
                y=0,  # -0.25,#-0.25,# was -0.5
                z=0)

            ## apoptosis_chance ##
            # if self.mcs % 10 == 0:
            # if cell.type in [ECTODERM_LEFT, ECTODERM_RIGHT, MESODERM_LEFT, MESODERM_RIGHT]:
            # coinflip = random.random()
            # if self.apoptosis_chance > coinflip:
            # cell.targetVolume = 0

        # for cell in self.cell_list:

        # # you can use here any fcn of concentrationAtCOM
        # cell.targetVolume += 0.01 * concentrationAtCOM

        # for cell in self.cell_list_by_type(MESODERM_LEFT, MESODERM_RIGHT):
        # concentrationAtCOM = field[int(cell.xCOM), int(cell.yCOM), int(cell.zCOM)]
        # if cell.volume < 200:
        # cell.targetVolume += (concentrationAtCOM/300)
        # # Make sure Secretion plugin is loaded
        # # make sure this field is defined in one of the PDE solvers
        # # you may reuse secretor for many cells. Simply define it outside the loop
        # secretor = self.get_field_secretor("shh")

        # # arguments are: cell, max uptake, relative uptake
        # secretor.uptakeInsideCell(cell, 2.0, (concentrationAtCOM/100))
        # tot_amount = secretor.uptakeInsideCellTotalCount(cell, 2.0, 0.2).tot_amount

        # if cell.dict["ENABLE_GROWTH"] == 1:
        # cell.targetVolume += 0.01
        # elif cell.dict["ENABLE_GROWTH"] == 2:
        # cell.targetVolume += 0.005
        # Make sure Secretion plugin is loaded
        # make sure this field is defined in one of the PDE solvers


class MitosisSteppable(MitosisSteppableBase):

    def __init__(self, frequency=10):
        MitosisSteppableClustersBase.__init__(self, frequency)
        # self.cells_to_divide = []
        self.do_not_divide_type_list = [MEDIUM, WALL]
        self.divide_cell_type_list = [MESODERM_LEFT, MESODERM_RIGHT, ECTODERM_LEFT, ECTODERM_RIGHT, NOTOCHORD]

    def step(self, mcs):
        cells_to_divide = []
        for cell in self.cellList:

            if cell.type in self.divide_cell_type_list:
                if cell.dict["ALLOW_DIVISION"] == True and cell.volume > 60:
                    cells_to_divide.append(cell)
                    cell.dict["ALLOW_DIVISION"] = False
                    cell.dict["CELL_CYCLE"] = random.randint(0, 1500)

            if cell.type in [ECM]:
                if cell.volume > 50:
                    cells_to_divide.append(cell)

        # mitosis_cluster_id_list = []
        # mitosis_cell_type_list = []
        # for compartment_list in self.clusterList:
        # # print( "cluster has size=",compartment_list.size())
        # cluster_volume = 0
        # cluster_id = 0
        # for cell in CompartmentList(compartment_list):
        # if cell.type in self.cluster_base_type_list:
        # if cell.dict["ALLOW_DIVISION"] == True:
        # cluster_volume += cell.volume
        # if cluster_volume >= 40:
        # cluster_id = cell.clusterId
        # mitosis_cluster_id_list.append(cluster_id)
        # mitosis_cell_type_list.append(cell.type)
        # cell.dict["ALLOW_DIVISION"] = False
        # cell.dict["CELL_CYCLE"] = random.randint(0,1500)

        for cell in cells_to_divide:
            if cell.type in [MESODERM_LEFT, MESODERM_RIGHT, NOTOCHORD]:
                self.divideCellRandomOrientation(cell)
            elif cell.type in [ECTODERM_LEFT]:
                self.divideCellRandomOrientation(cell)
                # self.divideCellOrientationVectorBased(cell,-1,0,0)
            elif cell.type in [ECTODERM_RIGHT]:
                self.divideCellRandomOrientation(cell)
                # self.divideCellOrientationVectorBased(cell,1,0,0)

        # for i in range(len(mitosis_cluster_id_list)):
        # if mitosis_cell_type_list[i] in [MESODERM_LEFT, MESODERM_RIGHT]:
        # self.divide_cluster_random_orientation(mitosis_cluster_id_list[i])
        # elif mitosis_cell_type_list[i] in [ECTODERM_LEFT]:
        # self.divide_cluster_orientation_vector_based(mitosis_cluster_id_list[i], -1, 0, 0)
        # elif mitosis_cell_type_list[i] in [ECTODERM_RIGHT]:
        # self.divide_cluster_orientation_vector_based(mitosis_cluster_id_list[i], 1, 0, 0)

        # # other valid options - to change mitosis mode leave one of the below lines uncommented
        # self.divide_cluster_orientation_vector_based(cluster_id, 1, 0, 0)
        # self.divide_cluster_along_major_axis(cluster_id)
        # self.divide_cluster_along_minor_axis(cluster_id)

    def update_attributes(self):
        # compartments in the parent and child clusters are
        # listed in the same order so attribute changes require simple iteration through compartment list
        # compartment_list_parent = self.get_cluster_cells(self.parent_cell.clusterId)
        parentCell = self.mitosisSteppable.parentCell
        childCell = self.mitosisSteppable.childCell
        save_volume = parentCell.targetVolume
        self.clone_parent_2_child()
        # self.clone_parent_cluster_2_child_cluster()#copies chemotaxis attributes too, don't want that in differentiation below
        parentCell.targetVolume = save_volume / 2  # random.randint(20,30) #save_volume / 2 # #distribute sizes so mitosis time is also random
        childCell.targetVolume = save_volume / 2  # max(20,save_volume-parentCell.targetVolume) #save_volume /2 #conservation of mass
        childCell.dict["CELL_CYCLE"] = random.randint(0, 1500)
        if childCell.type in [ECTODERM_LEFT]:
            childCell.dict["X_PUSH_DIRECTION"] = random.uniform(-1.0, 0.0)
            childCell.lambdaVecX = childCell.dict["X_PUSH_DIRECTION"]

        if childCell.type in [ECTODERM_RIGHT]:
            childCell.dict["X_PUSH_DIRECTION"] = random.uniform(0.0, 1.0)
            childCell.lambdaVecX = childCell.dict["X_PUSH_DIRECTION"]
        # childCell.lambdaVolume=2.0
        # parentCell.targetSurface =
        # childCell.lambdaSurface = 2.0
        # parentCell.targetSurface = 1.0 * (math.pi * (2*math.sqrt(parentCell.targetVolume / math.pi)))
        # parentCell.targetSurface=1.50*math.pow((36*math.pi*pow(parentCell.targetVolume, 2)), 1.0/3.0)
        # childCell.targetSurface=1.0 * (math.pi * (2*math.sqrt(childCell.targetVolume / math.pi)))

        # for i in range(len(compartment_list_parent)):
        # compartment_list_parent[i].targetVolume /= 2

        # 
# class MitosisSteppableClusters(MitosisSteppableClustersBase):

#     def __init__(self, frequency=10):
#         MitosisSteppableClustersBase.__init__(self, frequency)

#     def step(self, mcs):

#         parent_vector_list = []
#         mitosis_cluster_id_list = []
#         for cell in self.cell_list_by_type(NE_LATERAL_RIGHT):
#             cluster_cell_list = self.get_cluster_cells(cell.clusterId)
#             if len(cluster_cell_list) < 3:
#                 return "bug"


#             basal_cell = cluster_cell_list[0]
#             lateral_cell = cluster_cell_list[1]
#             apical_cell = cluster_cell_list[2]

#             # Make sure CenterOfMass plugin is loaded
#             # READ ONLY ACCESS
#             basal_xCOM = basal_cell.xCOM
#             basal_yCOM = basal_cell.yCOM
#             basal_zCOM = basal_cell.zCOM 

#             apical_xCOM = apical_cell.xCOM
#             apical_yCOM = apical_cell.yCOM
#             apical_zCOM = apical_cell.zCOM

#             print("basal_cells",[basal_xCOM, basal_yCOM, basal_zCOM])
#             print("apical_cells",[apical_xCOM, apical_yCOM, apical_zCOM])

#             parent_x = apical_xCOM - basal_xCOM
#             parent_y = apical_yCOM - basal_yCOM

#             #normalize vector
#             # combasal - com_apical

#             #divide by the square of the sum of the squares.

#             # deltaX DelataY delataZ


#             if mcs == 4900:
#                 lateral_cell.dict["CELL_CYCLE"] = random.randint(0,2000)


#         #random.shuffle(cluster_id_list)
#         #print(cluster_id_list)


#             if mcs > 500:
#                 if cell.dict["ENABLE_CONSTRICT"] == False:
#                     basal_cell.targetVolume += 0.2# = cluster_cell_list[0]
#                     lateral_cell.targetVolume += 0.2
#                     apical_cell.targetVolume += 0.2

#                     divide_cell = lateral_cell

#                     if divide_cell.dict["ALLOW_DIVISION"] == True and divide_cell.volume > 60:
#                 #divide_cell = cluster_id_list[0]
#                         mitosis_cluster_id_list.append(divide_cell)


#                         parent_vector_list.append((parent_x, parent_y))
#                         divide_cell.dict["ALLOW_DIVISION"] = False


#             #print("cluster_volume=", cluster_volume)

#             # condition under which cluster mitosis takes place
#             #if cluster_volume > 250:
#                 # instead of doing mitosis right away we store ids for clusters which should be divide.
#                 # This avoids modifying cluster list while we iterate through it
#                 #mitosis_cluster_id_list.append(cluster_id)

#         for i in range(len(mitosis_cluster_id_list)):
#             #if self.mcs % 1000 == 0:

#             self.divide_cluster_orientation_vector_based(
#             mitosis_cluster_id_list[i], 
#             parent_vector_list[i][0], 
#             parent_vector_list[i][1], 
#             0)


#     def update_attributes(self): #Cell attributes become mixed up?

#         # compartments in the parent and child clusters are
#         # listed in the same order so attribute changes require simple iteration through compartment list
#         compartment_list_parent = self.get_cluster_cells(self.parent_cell.clusterId)
#         compartment_list_child = self.get_cluster_cells(self.child_cell.clusterId)
#         basal_parent = compartment_list_parent[0]
#         lateral_parent = compartment_list_parent[1] 
#         apical_parent = compartment_list_parent[2]

#         basal_child = compartment_list_child[0]
#         lateral_child = compartment_list_child[1] 
#         apical_child = compartment_list_child[2]          


#         ## create FPP links for child based on parent links ##
#         # Basal - Apical link
#         p_b_a_link = self.get_fpp_internal_link_by_cells(basal_parent, apical_parent)
#         p_b_a_ld = p_b_a_link.getLambdaDistance()
#         p_b_a_td = p_b_a_link.getTargetDistance()  
#         p_b_a_md = p_b_a_link.getMaxDistance()   
#         link = self.new_fpp_internal_link(basal_child, apical_child, p_b_a_ld, p_b_a_td, p_b_a_md)

#         # Lateral - Apical link
#         p_l_a_link = self.get_fpp_internal_link_by_cells(lateral_parent, apical_parent)
#         p_l_a_ld = p_l_a_link.getLambdaDistance()
#         p_l_a_td = p_l_a_link.getTargetDistance()  
#         p_l_a_md = p_l_a_link.getMaxDistance()
#         link = self.new_fpp_internal_link(lateral_child, apical_child, p_l_a_ld, p_l_a_td, p_l_a_md)

#         # Lateral - Basal link
#         p_l_b_link = self.get_fpp_internal_link_by_cells(lateral_parent, basal_parent)
#         p_l_b_ld = p_l_b_link.getLambdaDistance()
#         p_l_b_td = p_l_b_link.getTargetDistance()  
#         p_l_b_md = p_l_b_link.getMaxDistance()
#         link = self.new_fpp_internal_link(lateral_child, basal_child, p_l_b_ld, p_l_b_td, p_l_b_md)

#         # remove all parent links and create new internal links
#         self.remove_all_cell_fpp_links(apical_parent)
#         link = self.new_fpp_internal_link(basal_parent, apical_parent, p_b_a_ld, p_b_a_td, p_b_a_md)
#         link = self.new_fpp_internal_link(lateral_parent, apical_parent, p_l_a_ld, p_l_a_td, p_l_a_md)
#         link = self.new_fpp_internal_link(lateral_parent, basal_parent, p_l_b_ld, p_l_b_td, p_l_b_md)


#         # apical parent - apical child link
#         apical_link = self.new_fpp_link(apical_child, apical_parent, 20, 4, 20)
#         apical_link = self.get_fpp_internal_link_by_cells(lateral_parent, basal_parent)
#         apical_link.setMaxNumberOfJunctions(2)

#         # Other neighbor links
#         for neighbor, common_surface_area in self.get_cell_neighbor_data_list(apical_parent):
#             if neighbor:
#                 if neighbor.type == apical_parent.type:
#                     link = self.get_fpp_link_by_cells(apical_parent, neighbor)
#                     if not link:
#                         apical_link = self.new_fpp_link(apical_parent, neighbor, 20, 4, 20)
#                         apical_link = self.get_fpp_link_by_cells(apical_parent, neighbor)
#                         apical_link.setMaxNumberOfJunctions(2)

#         for neighbor, common_surface_area in self.get_cell_neighbor_data_list(apical_child):
#             if neighbor:
#                 if neighbor.type == apical_child.type:
#                     link = self.get_fpp_link_by_cells(apical_child, neighbor)
#                     if not link:
#                         apical_link = self.new_fpp_link(apical_child, neighbor, 20, 4, 20)
#                         apical_link = self.get_fpp_link_by_cells(apical_child, neighbor)
#                         apical_link.setMaxNumberOfJunctions(2)


#         print(compartment_list_parent[0])
#         print(len(compartment_list_child))

#                             # fpp_link = self.get_fpp_internal_link_by_cells(basal_cell, apical_cell)
#                     # target_distance = fpp_link.getTargetDistance()
#                     # if target_distance < self.apical_constriction_elongation:
#                         # fpp_link.setTargetDistance(target_distance + self.delta_elongation)

#         print([x.type for x in compartment_list_child])
#         print([x for x in compartment_list_child])

#         for i in range(len(compartment_list_parent)):
#             compartment_list_parent[i].targetVolume /= 2.0
#         self.clone_parent_cluster_2_child_cluster()

#         lateral_parent.dict["CELL_CYCLE"] = random.randint(0,1000)
#         lateral_child.dict["CELL_CYCLE"] = random.randint(0,1000)
