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
using CSV


include("C:/Workdir/Develop/PF_simulations/My_functions.jl")
include("C:/Workdir/Develop/PF_simulations/My_ref/My_ref_lazy_DG_curtail.jl")
include("HC_function_curt.jl")

#Parameters

flex = 20                              # flexibility offered %
congestion_limit = 100  	             # congestion limit %
threshold = 100                        # value used to identify branches with current rating higher than threshold
curtailment = 0                       # DG curtailment %

seed = 99                              # seed for random DG buses choice
gen_number_per_feeder = 5              # number of random DGs per feeder
power_target_per_feeder = 5            # total capacity installed in each feeder

gen_step = 0.1/5                            # Incremental size of generation in MW

# Input file
file_name = "Official_semiurban.m"
file_path = "C://Users//u0152683//Desktop//Networks//Experiments//Official_semiurban.m"
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
gen_ID = get_gen_info(net_data, feeder_ID_1)

HC = power_target_per_feeder*length(feeder_ID_1) # initial HC value 

# Find HC
feeder_HC = Dict()
healthy_feeders = []

feeder_names = Vector{String}()
[push!(feeder_names, feeder["Name"]) for (i,feeder) in feeder_ID_1]
sort!(feeder_names)
healthy_feeders = copy(feeder_names)  # store only feeder names with no violations

plot_info = Dict{Any, Dict{Any, Any}}()
iter = 0
stable = true  #if stable = false, it means ALL feeders have reached network violations

for (k, feeder) in feeder_ID_1  # Aggiungi appartenenza ad ogni feeder per bus, lines e load
    for (id, branch) in net_data["branch"]
        if id in feeder["Branches"]
            net_data["branch"][id]["feeder"] = feeder["Name"]
        end
    end
    for (idx, bus) in net_data["bus"]
        if parse(Int64,idx) in feeder["Buses"][2:end]
            net_data["bus"][idx]["feeder"] = feeder["Name"]
        end
    end
    for (idx, load) in net_data["load"]
        if load["load_bus"] in feeder["Buses"][2:end]
            net_data["load"][idx]["feeder"] = feeder["Name"]
        end
    end
end

[net_data["branch"][id]["feeder"] = "HV/MV" for (id,branch) in net_data["branch"] if !haskey(branch,"feeder")]
[net_data["bus"][id]["feeder"] = "HV/MV" for (id,bus) in net_data["bus"] if !haskey(bus,"feeder")]

while stable 
    stable = false
    iter +=1
    println("\n\n Starting iteration $iter: Capacity installed = $HC MW\n" ) 
    
    push!(plot_info, HC =>Dict(f_name => Dict("p_installed" => 0.0, "p_generated"=> 0.0, 
    "curtailment"=> 0.0, "p_load_nominal" => 0.0, "p_load_final"=> 0.0, "q_load_nominal" => 0.0, 
    "q_load_final"=> 0.0, "p_flex_up" => 0.0,"p_flex_dwn" => 0.0, "q_flex_up" => 0.0,"q_flex_dwn" => 0.0, "n_volt_vio" =>0, "n_cong_vio" => 0 ) for f_name in healthy_feeders))

    [plot_info[HC][gen["feeder"]]["p_installed"] += gen["p_nominal"] for (i,gen) in gen_ID if gen["feeder"] in healthy_feeders]

    [plot_info[HC][load["feeder"]]["p_load_nominal"] += load["pd"] for (i,load) in net_data["load"] if load["feeder"] in healthy_feeders]
    [plot_info[HC][load["feeder"]]["q_load_nominal"] += load["qd"] for (i,load) in net_data["load"] if load["feeder"] in healthy_feeders]


    println("Running PF with: \n")
    [println(f_name," ", round(sum([gen["p_nominal"] for (id, gen) in gen_ID if gen["feeder"] == f_name]), digits = 4), " MW") for f_name in feeder_names]

    global result = test(net_data, ACPPowerModel, build_pf_dg_curtail ,gen_ID, feeder_ID_1, feeder_HC, gen_step, healthy_feeders, feeder_names, plot_info, HC)
    
    if !isempty(healthy_feeders) #healthy_feeders gets updated as soon as a feeder violates constraints
        stable = true
        feeders_iterator = copy(healthy_feeders)

        for f_name in feeders_iterator

            gen_increase_size(net_data, gen_ID, f_name, gen_step)
            println("\n Increasing capacity generation at $f_name of ", gen_step*gen_number_per_feeder," MW")
            HC += gen_step*gen_number_per_feeder   
            
        end

    else
        HC = sum([feeder["HC"] for (id,feeder) in feeder_HC])
        print("\n Hosting Capacity evaluation completed!\n Maximum Network HC: $HC MW \n")
    end
    
end


update_data!(net_data, result["solution"])

plot_info = sort(plot_info)

[println(id," HC: ", feeder["HC"], " MW") for (id,feeder) in sort(feeder_HC)]

#plot_grid_new(net_data, "p_flex", "curtailment")

for f_name in feeder_names

    df = DataFrame()

    q_load_final = Vector()
    p_load_final = Vector()
    curtailment = Vector()
    p_installed = Vector()
    p_load_nominal = Vector()
    q_load_nominal = Vector()
    p_generated = Vector()
    n_volt_vio = Vector()
    n_cong_vio = Vector()
    p_flex_up = Vector()
    q_flex_up = Vector()
    p_flex_dwn = Vector()
    q_flex_dwn = Vector()

    for (HC,x) in plot_info
        if haskey(x,f_name)
            dati = x[f_name]

            push!(q_load_final, dati["q_load_final"])
            push!(p_load_final, dati["p_load_final"])
            push!(curtailment, dati["curtailment"])
            push!(p_installed, dati["p_installed"])
            push!(q_load_nominal, dati["q_load_nominal"])
            push!(p_load_nominal, dati["p_load_nominal"])
            push!(p_generated, dati["p_generated"])
            push!(n_volt_vio, dati["n_volt_vio"])
            push!(n_cong_vio, dati["n_cong_vio"])
            push!(p_flex_up, dati["p_flex_up"])
            push!(q_flex_up, dati["q_flex_up"])
            push!(p_flex_dwn, dati["p_flex_dwn"])
            push!(q_flex_dwn, dati["q_flex_dwn"])
        end
    end

    df.Q_load_final = round.(q_load_final, digits = 4)
    df.P_load_final = round.(p_load_final, digits = 4)
    df.Curtailment = round.(curtailment, digits = 4)
    df.P_installed = round.(p_installed, digits = 4)
    df.P_load_nominal = round.(p_load_nominal, digits = 4)
    df.Q_load_nominal = round.(q_load_nominal, digits = 4)
    df.P_generated = round.(p_generated, digits = 4)
    df.N_volt_vio = n_volt_vio
    df.N_cong_vio = n_cong_vio
    df.P_flex_up = round.(p_flex_up, digits = 4)
    df.Q_flex_up = round.(q_flex_up, digits = 4)
    df.P_flex_dwn = round.(p_flex_dwn, digits = 4)
    df.Q_flex_dwn = round.(q_flex_dwn, digits = 4)

    file_name = replace(file_name, "Official_" => "")
    file_name = replace(file_name, ".m" => " HC")
    file_name = uppercasefirst(file_name)

    #XLSX.openxlsx("test3.xlsx", mode="rw") do xf
    XLSX.openxlsx(file_name*".xlsx", mode="rw") do xf
    
        try
            XLSX.addsheet!(xf, f_name)
        catch
        end
        
        sheet = xf[f_name]

        column_names = names(df)

        for c in 1:length(column_names)
            sheet[XLSX.CellRef(1,c)] = column_names[c]
        end

        for r in 1:size(df,1), c in 1:size(df,2)
            sheet[XLSX.CellRef(r+1,c)] = df[r,c]
        end
    end
end

        
