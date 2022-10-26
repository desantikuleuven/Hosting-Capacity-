using PowerModels
using PowerModelsAnalytics
using PowerPlots
using ColorSchemes
using Setfield
using Plots
using JuMP, Ipopt
using StatsBase
using Random
import InfrastructureModels
const _IM = InfrastructureModels
using DataFrames
using XLSX

include("C:/Workdir/Develop/PF_simulations/My_functions.jl")
include("C:/Workdir/Develop/PF_simulations/My_ref/My_ref_lazy_flex.jl")
include("HC_functions.jl")

#Parameters

flex = 20                              # flexibility offered %
congestion_limit = 100  	             # congestion limit %
threshold = 100                        # value used to identify branches with current rating higher than threshold

seed = 99                              # seed for random DG buses choice
gen_number_per_feeder = 5              # number of random DGs per feeder
power_target_per_feeder = 5            # total capacity installed in each feeder

step = 1.5                             # Incremental size of generation in MW

# Input file
file_name = "Official_urban.m"
file_path = "C://Users//u0152683//Desktop//Networks//Experiments//Official_urban.m"
net_data = parse_file(file_path)

#Add flexibility % that each load can offer
net_data["flex"] = flex
update_data!(net_data, add_load_flexibility(net_data, net_data["flex"]/100))

# Add congestion capacity for each line
cong_cap = add_congestion_capacity(net_data, congestion_limit/100)
update_data!(net_data, cong_cap)

# Get feeder info 
feeder_ID_1, mv_busbar = get_feeder_data(net_data, file_name)

# Random choice of buses
random_generators = get_random_generators(feeder_ID_1, gen_number_per_feeder, seed)  

# Add new generators
size_std = power_target_per_feeder/gen_number_per_feeder  #size of each DG
add_generators(net_data, random_generators, size_std)
gen_ID = get_gen_info(net_data, feeder_ID_1)

HC = power_target_per_feeder*length(feeder_ID_1) # initial HC value 

# Find HC
feeder_HC, gen_ID = calc_hosting_capacity_w_flex(net_data, build_pf_all_flex_lazy, HC, step ; feeder_ID_1, gen_ID, gen_number_per_feeder )

# Print results

[println(id," HC: ", feeder["HC"], " MW") for (id,feeder) in sort(feeder_HC)]

# Compute voltage profiles
feeder_ID, path_volt, mv_busbar = calc_voltage_profile(net_data, result, file_name)
calc_branch_loading(net_data, feeder_ID, gen_ID, threshold)

# Evaluate types of violations occurred in feeders
add_violation_type_flex(net_data, feeder_ID, feeder_HC)

# Claculate flexibility offered 
flex_loads, p_load, q_load = calc_flexibility_offered_new(net_data, result)

# Compute branch flows and losses

flows = calc_branch_flow_ac(net_data)
update_data!(net_data, flows)

losses, p_loss, q_loss = calc_power_losses(net_data)
update_data!(net_data, losses)

# Printing statements
printing_statements_HC_flex(feeder_HC, feeder_ID, flex_loads, p_load, q_load)

create_df_HC_flex(feeder_HC)

#Look at function description in My_functions to see which argument you can pass
#bus, gen, branch
plot_grid(net_data, "p_flex","basic", "loading"; display_flow = true)


