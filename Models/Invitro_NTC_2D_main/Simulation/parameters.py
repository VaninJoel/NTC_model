## general cell parameters        
standard_target_volume = 49  # initial volume of the cells based on the 7x7 cell layout
standard_lambda_volume = 2.0  # lambda volume
time_to_stop_mesodermal_growth = 4000  # after x mcs
time_to_stop_ectodermal_growth = 4000  # after x mcs
cell_cycle_time = 4000  # number of mcs for once cell cycle to complete
mesoderm_growth_rate = 0.01  # delta volume that is added per mcs in mesoderm
mesoderm_max_growth_rate = 1.0  # maximum growth rate, input for hill function
ectoderm_growth_rate = 0.0010  # delta volume that is added per mcs in ectoderm
somite_growth_rate = 0.0  # delta volume that is added per mcs in the somite
ECM_growth_rate = 0.01  # delta volume that is added per mcs in the ECM
mesodermal_differentiation_probability = 5 / 1000  # probability for a mesodermal cell to differentiate into somite
apoptosis_chance = 1 / 10000  # probability for a cell to undergo apoptosis

## gene parameters ## Include genes from Kazburn et al, 2021
# ATRA --> Time based parameter, starts Low, ends high
# FGF --> Time based variable, starts high, ends low.
# NOG --> Inhibits BMP
# N/E Cadherin switch --> for fusion, regulated by Snail2
# Snail2 is regulated by WNT
gene_status_dict = {  # genes and their expression (attach to each cell?)
    "BMP4": 1.0,
    "SHH": 1.0
}

## BMP4  parameters ##
BMP4_secretion_constant = 0.020  # secretion constant of BMP4, (can be modified by changing the status dict)
AC50_BMP4 = 50  # AC50 BMP (arbitraty value)
BMP4_gradient_low = 0.1  # lowerbound for apical constriction
BMP4_gradient_high = 0.8  # upperbound for apical constriction

## SHH ##
SHH_secretion_constant = 0.010  # secretion constant of SHH, (can be modified by changing the status dict)
SHH_gradient_low = 0.1  # lowerbound for apical constriction
SHH_gradient_high = 0.6  # upperbound for apical constriction

## Toxicological perturbation parameters ##
dose_toxic_substance = False  # set to true to induce toxic field that kills NE

##Apical constriction parameters##
enable_apical_constriction_behavior = True
enable_mhp_formation = True
ne_apical_target_volume = 21  # The target volume of apical cell under apical constriction
ne_basal_target_volume = 77  # The target volume of basal cell under apical constriction
mhp_apical_target_volume = 21
mhp_basal_target_volume = 77
apical_constrict_target_volume = 21
basal_constrict_target_volume = 77

mhp_constrict_frequency = 1  # per how many mcs a apical constriction step is performed
apical_constrict_frequency = 1  # The number of MCS for one step in apical constriction
delta_target_volume = 0.04  # the amount the target volume changes per step frequency
delta_lambda_volume = 0.008  # makes the small cells not dissapear # the amount the lambda volume changes per step frequency
delta_target_volume_apical = 0.04
apical_constriction_elongation = 18  # the target distance of fpp links between the basal and apical cells under apical constriction
delta_elongation = 0.01  # the step size of elongation difference
standard_elongation = 12  # the standard target distance of fpp links between the basal and apical cells

## Apical constriction FPP links parameters
initial_lambda_distance = 20.0  # The initial lambda distance of FPP links
initial_target_distance = 5.0  # The initial target distance of FPP links
initial_max_distance = 20.0  # The initial max distance of FPP links
apical_constrict_lambda_distance = 40.0  # The lambda distance of FPP links under Apical constriction
apical_constrict_target_distance = 0.0  # The lambda distance of FPP links under Apical constriction
basal_constrict_target_distance = 9.0  # FPP link distance on the basal side
apical_constrict_max_distance = 15.0  # The max distance of FPP links under Apical constriction
delta_distance = 0.005  # delta change of the FPP link distance

## Cell layout parameters
use_ECM = True  # initialize extra ECM around the mesodermal cells

##### for the cluster computer #####

BMP4_secretion_constant = 0.020
SHH_secretion_constant = 0.010
mesoderm_growth_rate = 0.010

#####
