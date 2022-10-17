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

#=
COMPUTE HOSTING CAPACITY WITHOUT FELXIBILITY 
=#

include("C:/Users/u0152683/Desktop/Networks/PF simulation/My_functions.jl")
#include("C:/Users/u0152683/Desktop/Networks/PF simulation/My_ref/My_ref_lazy_ALL.jl")
include("HC_functions.jl")
include("Ref_no_flex.jl")

# Input file
file_name = "Official_rural.m"
file_path = "C://Users//u0152683//Desktop//Networks//Experiments//Official_rural.m"
net_data = parse_file(file_path)

# Add congestion capacity for each line
congestion_limit = 100
threshold = 100
cong_cap = add_congestion_capacity(net_data, congestion_limit/100)
update_data!(net_data, cong_cap)

# Get feeder info 

feeder_ID_1, mv_busbar = get_feeder_data(net_data, file_name)

# Random choice of buses
gen_number_per_feeder = 5 #number of random DGs per feeder
random_generators = get_random_generators(feeder_ID_1, gen_number_per_feeder, 99)  

# Add new generators
power_target_per_feeder = 5
size_std = power_target_per_feeder/gen_number_per_feeder  #size of each DG
add_generators(net_data, random_generators, size_std)
gen_ID = get_gen_info(net_data, feeder_ID_1)

HC = power_target_per_feeder*length(feeder_ID_1) # initial HC value 

feeder_HC, gen_ID = calc_hosting_capacity(net_data, build_ref_no_flex, HC, 0.5 , threshold; feeder_ID_1, gen_ID, file_name, gen_number_per_feeder )

[println(id," HC: ", feeder["HC"], " MW") for (id,feeder) in sort(feeder_HC)]

# Update net_data with branch loadings 
#calc_branch_loading e calc_voltage_profile already performed inside the HC function

#plot_grid_new(net_data, "basic")

#= Compute flexibility offered by each load 
flex_loads, p_load, q_load = calc_flexibility_offered_new(net_data, result)
=#
# Compute branch flows and losses

flows = calc_branch_flow_ac(net_data)
update_data!(net_data, flows)

losses, p_loss, q_loss = calc_power_losses(net_data)
update_data!(net_data, losses)


for (id,feeder) in sort(feeder_HC)
    println(" $id HC: ",feeder["HC"], " MW")
    if Bool(feeder["Volt_violation"]) && Bool(feeder["Congestion_violation"])
        println(" Limiting factor for $id: voltage bounds and congestion constraint \n")
    elseif !Bool(feeder["Volt_violation"])
        println(" Limiting factor for $id: thermal limits \n")
    else
        println(" Limiting factor for $id: voltage bounds \n") 
    end
end   

    

print("\n Slack bus power generated: ", result["solution"]["gen"]["1"]["pg"], "MW \n")
print("\n Immaginary power generated: ", result["solution"]["gen"]["1"]["qg"], "MVar \n")

print("\n MV busbar voltage : ", round(result["solution"]["bus"]["$mv_busbar"]["vm"],digits =5), "p.u. \n")

#print("\n Real power requested: ", p_load, "MW \n")
#print("\n Immaginary power requested: ", q_load, "MVar \n")

print("\n Real power losses: ", sum(p_loss), "MW \n")
print("\n Immaginary power losses: ", sum(q_loss), "MVar \n")

#[println("\n Vmin: ", round(feeder_ID[gen["ref"]]["vmin"], digits = 5)," Vmax: ",round(feeder_ID[gen["ref"]]["vmax"],digits = 5)," for ",gen["feeder"]," (where generator $i is connected)") for (i,gen) in gen_ID]
#[println("\n Feeder loading: ", round(gen["max_branch_load"],digits = 3),"% in branch ",gen["Critical_branch"]) for (i,gen) in gen_ID]



for (id, gen) in sort(gen_ID)
    println(" Generator $id is connected at bus ", gen["bus"], " in ", gen["feeder"], " with ",gen["pg"],"MW installed")
end

println()
#=
for f_name in feeder_names

    for (ref,data) in feeder_ID
        if data["Name"] == f_name

            println(" ",f_name," --> Vmin: ", round(data["vmin"], digits = 5), " Vmax: ", round(data["vmin"], digits = 5))
            
        end
    end
end
=#
println()


vmin = Vector()
vmax = Vector()
[push!(vmin,feeder["vmin"]) for (id, feeder) in feeder_ID ]
println("\n Lowest voltage magnitude in the grid: ", minimum(vmin))
[push!(vmin,feeder["vmax"]) for (id, feeder) in feeder_ID ]
println("\n Highest voltage magnitude in the grid: ", maximum(vmin))

plot_grid_new(net_data, "vm", true)
