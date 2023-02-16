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
include("C:/Workdir/Develop/PF_simulations/My_ref/My_ref_lazy_DG_curtail.jl")
include("HC_function_curt.jl")

#Parameters

flex = 20                              # flexibility offered %
congestion_limit = 100  	             # congestion limit %
threshold = 100                        # value used to identify branches with current rating higher than threshold
curtailment = 20                       # DG curtailment %

seed = 99                              # seed for random DG buses choice
gen_number_per_feeder = 5              # number of random DGs per feeder
power_target_per_feeder = 5            # total capacity installed in each feeder

gen_step = 1.5                             # Incremental size of generation in MW

# Input file
file_name = "Official_urban.m"
file_path = "C:/Workdir//Develop//"*file_name
net_data = parse_file(file_path)

#Add flexibility % that each load can offer
net_data["flex"] = flex
net_data["curtail"] = curtailment
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
add_generators(net_data, random_generators, size_std, curtailment/100)

# Create dict for DGs (CALL ONLY BEFORE PF)
gen_ID = get_gen_info(net_data, feeder_ID_1)

HC = power_target_per_feeder*length(feeder_ID_1) # initial HC value 

# Find HC
feeder_HC = Dict()

healthy_feeders = []

feeder_names = Vector{String}()
[push!(feeder_names, feeder["Name"]) for (i,feeder) in feeder_ID_1]
sort!(feeder_names)
healthy_feeders = copy(feeder_names)  # store only feeder names with no violations

iter = 0
stable = true  #if stable = false, it means ALL feeders have reached network violations

while stable
    stable = false
    iter +=1
    println("\n\n Starting iteration $iter: Capacity installed = $HC MW\n" )

    println("Running PF with: \n")
    [println(f_name," ", round(sum([gen["p_nominal"] for (id, gen) in gen_ID if gen["feeder"] == f_name]), digits = 4), " MW") for f_name in feeder_names]

    global result = HC_evaluation_flex_curt(net_data, ACPPowerModel, build_pf_dg_curtail ,gen_ID, feeder_ID_1, feeder_HC, gen_step, healthy_feeders, feeder_names)
    
    if !isempty(healthy_feeders) #healthy_feeders gets updated as soon as a feeder violates constraints
        stable = true
        feeders_iterator = copy(healthy_feeders)

        for f_name in feeders_iterator

            gen_increase_size(net_data, gen_ID, f_name, gen_step)
            println("\n Increasing capacity generation at $f_name of ",gen_step * gen_number_per_feeder," MW")
            HC += gen_step * gen_number_per_feeder
            
        end

    else
        HC = sum([feeder["HC"] for (id,feeder) in feeder_HC])
        print("\n Hosting Capacity evaluation completed!\n Maximum Network HC: $HC MW \n")
    end
    
end

update_data!(net_data, result["solution"])

# Print results

[println(id," HC: ", feeder["HC"], " MW") for (id,feeder) in sort(feeder_HC)]

# Voltage profiles: specify path if you want to save
save_path =  ""
feeder_ID, path_volt, mv_busbar = calc_voltage_profile(net_data, result, file_name, save_path, true, true, false)

# Add branch loadings to net_data
calc_branch_loading(net_data, feeder_ID, gen_ID, threshold)

# Compute branch flows and losses

flows = calc_branch_flow_ac(net_data)
update_data!(net_data, flows)

losses, p_loss, q_loss = calc_power_losses(net_data)
update_data!(net_data, losses)

# Evaluate types of violations occurred in feeders
add_violation_type_flex(net_data, feeder_ID, feeder_HC)

# Claculate flexibility offered 
flex_loads, p_load, q_load = calc_flexibility_offered_new(net_data, result)

# Printing statements
printing_statements_HC_flex_curt(feeder_HC, feeder_ID, net_data, result)

# Add curtailment info to gen_ID, net_data and create new dict 
feeder_curtailment = feeder_dg_curtailment(net_data, result, gen_ID)

create_df_HC_flex_curt(feeder_HC, feeder_curtailment)

#Look at function description in My_functions to see which argument you can pass
#bus, gen, branch
plot_grid(net_data, "p_flex","curtailment","loading"; zoom =false, display_flow = false, save_fig = false)


