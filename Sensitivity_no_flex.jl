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
HOSTING CAPACITY SENSITIVITY WITHOUT FLEXIBILITY 
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

df = DataFrame()

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

    # Find HC
    feeder_HC, gen_ID = calc_hosting_capacity(net_data, build_ref_no_flex, HC, 0.1 , threshold; feeder_ID_1, gen_ID, file_name, gen_number_per_feeder )

    #create df
    df = hcat(df, create_df_HC(feeder_HC), makeunique = true)

    for (id,gen) in gen_ID
        if id != "1"
            delete!(gen_ID,id)
            delete!(net_data["gen"],id)
        end
    end

end

XLSX.openxlsx("C:\\Users\\u0152683\\Desktop\\Networks\\Hosting Capacity\\Sensitivity.xlsx", mode="rw") do xf
    3#XLSX.addsheet!(xf,"No_flex")
    sheet = xf["No_flex"]
    for r in 1:size(df,1), c in 1:size(df,2)
        sheet[XLSX.CellRef(r+1,c)] = df[r,c]
    end
end

