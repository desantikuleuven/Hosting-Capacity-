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

#=
HOSTING CAPACITY SENSITIVITY with FLEXIBILITY 
=#

include("C:/Users/u0152683/Desktop/Networks/PF simulation/My_functions.jl")
include("C:/Users/u0152683/Desktop/Networks/PF simulation/My_ref/My_ref_lazy_ALL.jl")
include("HC_functions.jl")
include("Ref_no_flex.jl")

# Input file
file_name = "Official_rural.m"
file_path = "C://Users//u0152683//Desktop//Networks//Experiments//Official_rural.m"
net_data = parse_file(file_path)

##### PARAMETERS #####

# Flexibility 
flx = [10,15,20,25,30]

# Generators number per feeder
gen_number = [5, 8, 12]


# Initial value = 5MW
# Step = 0.1 MW
# seed = 99

#######################

# Add congestion capacity for each line
congestion_limit = 100
threshold = 100
cong_cap = add_congestion_capacity(net_data, congestion_limit/100)
update_data!(net_data, cong_cap)

# Get feeder info 
feeder_ID_1, mv_busbar = get_feeder_data(net_data, file_name)

for flex_value in flx

    #Add flexibility % that each load can offer
    net_data["flex"] = flex_value
    flex = add_load_flexibility(net_data, net_data["flex"]/100)
    update_data!(net_data, flex)
    df = DataFrame()

    println("\n INITIALIZING ANALYSIS WITH $flex_value % OF FLEXIBILITY AVAILABLE ")

    for gnf in gen_number

        # Random choice of buses
        gen_number_per_feeder = gnf #number of random DGs per feeder
        random_generators = get_random_generators(feeder_ID_1, gen_number_per_feeder, 99)  

        # Add new generators
        power_target_per_feeder = 5
        size_std = power_target_per_feeder/gen_number_per_feeder  #size of each DG
        add_generators(net_data, random_generators, size_std)
        gen_ID = get_gen_info(net_data, feeder_ID_1)

        HC = power_target_per_feeder*length(feeder_ID_1) # initial HC value 

        println("\n INTALLED GENERATORS PER FEEDER: $gnf ")

        # Find HC
        feeder_HC, gen_ID = calc_hosting_capacity_w_flex(net_data, build_pf_all_flex_lazy, HC, 0.1 ; feeder_ID_1, gen_ID, gen_number_per_feeder )

        # Evaluate types of violations occurred in feeders
        add_violation_type_flex(net_data, feeder_ID_1, feeder_HC)

        # Claculate flexibility offered 
        flex_loads, p_load, q_load = calc_flexibility_offered_new(net_data, result)

        for (x, feeder) in feeder_ID_1
            f_name = feeder["Name"]
        
            feeder_HC[f_name]["up_p_flex"] = 0
            feeder_HC[f_name]["up_q_flex"] = 0
            feeder_HC[f_name]["down_p_flex"] = 0
            feeder_HC[f_name]["down_q_flex"] = 0 
            
            for bus in feeder["Buses"]
                if flex_loads["$bus"]["diff_real"] >= 0
                    feeder_HC[f_name]["up_p_flex"] += flex_loads["$bus"]["diff_real"]
                else 
                    feeder_HC[f_name]["down_p_flex"] += flex_loads["$bus"]["diff_real"]
                end
                if flex_loads["$bus"]["diff_imm"] >= 0
                    feeder_HC[f_name]["up_q_flex"] += flex_loads["$bus"]["diff_imm"]
                else
                    feeder_HC[f_name]["down_q_flex"] += flex_loads["$bus"]["diff_imm"]
                end
        
            end
        end

        #create df
     
        
        df = hcat(df, create_df_HC_flex(feeder_HC), makeunique = true)
        
        for (id,gen) in gen_ID
            if id != "1"
                delete!(gen_ID,id)
                delete!(net_data["gen"],id)
            end
        end

    end

    @eval $(Symbol(:df, flex_value)) = df

    XLSX.openxlsx("C:\\Users\\u0152683\\Desktop\\Networks\\Hosting Capacity\\Sensitivity.xlsx", mode="rw") do xf
        #XLSX.addsheet!(xf,"$flex_value")
        sheet = xf["$flex_value"]
        for r in 1:size(df,1), c in 1:size(df,2)
            sheet[XLSX.CellRef(r+2,c)] = df[r,c]
        end
    end

end




#= Print results

[println(id," HC: ", feeder["HC"], " MW") for (id,feeder) in sort(feeder_HC)]

# Compute voltage profiles
feeder_ID, path_volt, mv_busbar = calc_voltage_profile(net_data, result, file_name)
calc_branch_loading(net_data, feeder_ID, gen_ID, threshold)

# Claculate flexibility offered 
flex_loads, p_load, q_load = calc_flexibility_offered_new(net_data, result)

# Compute branch flows and losses

flows = calc_branch_flow_ac(net_data)
update_data!(net_data, flows)

losses, p_loss, q_loss = calc_power_losses(net_data)
update_data!(net_data, losses)

# Printing statements
printing_statements_HC_flex(feeder_HC, feeder_ID, flex_loads, p_load, q_load)


#plot_grid_new(net_data, "p_flex")
#plot_grid_new(net_data, "vm" )

create_df_HC(feeder_HC)
=#