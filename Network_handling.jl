import Pkg
#Pkg.add("FileIO")
#Pkg.add("PowerModelsAnalytics")
#Pkg.add("ColorSchemes")
#Pkg.add("Setfield")
#Pkg.add("Ipopt")
using PowerModels
using PowerModelsAnalytics
using PowerPlots
using DataFrames
using ColorSchemes
using Setfield
using JuMP, Ipopt
using FileIO
#=
1) lump_load, Methodology:
    - From the list of MV and LV, retrives which LV nodes are connected to the MV nodes by creating a directed graph
    - The load of each LV feeder is added to the MV node to which is connected 

    !! If eliminate_LV_branches is not applied afterwards, then LV nodes are still present in the network, hence
    we have just increased the load at th MV network. !!

2) map_nodes, inverse_map: needed in case the index of the bus does not match with the bus number 

3) eliminate_LV_branches, Methodology: 
    - Distinguish between MV and LV nodes
    - All proprieties linked to LV nodes (i.e. load, bus, etc) are eliminated from the data

=#

function lump_load(net_data::Dict)

    list_LV = []
    list_MV = []
    
    # GET LIST OF MV AND LV BUSES
    
    for (ind, bus) in net_data["bus"]
    
        if bus["base_kv"] < 10  # save buses whose voltage is higher than 10 (i.e. MV)
            push!(list_LV,bus["bus_i"])
        else
            push!(list_MV,bus["bus_i"])
        end
    end
    for i in keys(net_data["branch"])
       
        if net_data["branch"][i]["br_status"] == 0
            pop!(net_data["branch"],i)
        end
    end
    
    list_LV=sort(list_LV)
    list_MV = sort(list_MV)
    
    # Get list of secondary substations: MV bus --> LV bus
    
    transformers = sort(Dict( bus["f_bus"] => bus["t_bus"] for (i,bus) in net_data["branch"] if bus["transformer"] && bus["t_bus"] in list_LV && bus["f_bus"] in list_MV))
    
    MV_transformers = Int64[]
    [push!(MV_transformers, key) for key in keys(transformers)] # Save list of MV buses that are secondary substations
    
    # CREATE DIRECTED GRAPH
    
    ref = PowerModels.build_ref(net_data)[:it][:pm][:nw][0]
    ref[:load] = sort(ref[:load])
    arcs_from = sort(ref[:arcs_from])
    n_busses = length(keys(ref[:bus]))
    map = map_nodes(net_data)
    Gdir = LightGraphs.SimpleDiGraph(n_busses)
    
    if is_cyclic(Gdir) 
        print("ERROR cyclic network")
    end
    
    # Add arcs in directed Graph
    
    for (l, i, j) in arcs_from
    
        if i in list_MV && j in list_MV  #Consider only branches between MV-LV buses and LV-LV buses
            continue
        else
            i = map[i]
            j = map[j]
            LightGraphs.add_edge!(Gdir, i, j) 
        end
    end 
    
    # down_MV gives me the LV node on the secondary side of the trafo
    # down_LV gives me all the LV nodes below the LV node of the trafo 
    
    down_MV = sort(Dict(b => [LightGraphs.dst(e) for e in collect(LightGraphs.edges(LightGraphs.dfs_tree(Gdir, b; dir=:out)))] for b in MV_transformers))
    
    down_LV = sort(Dict(LV[1] => [LightGraphs.dst(e) for e in collect(LightGraphs.edges(LightGraphs.dfs_tree(Gdir, LV[1]; dir=:in)))] for (MV, LV) in down_MV))
    #=
    # SECURITY CHECKS
    
    all_lv_nodes = sort(Dict( i => 0 for i in keys(ref[:bus]) if i in list_LV)) # used to check if LV nodes are double counted
    
    for (lv_bus, nodes) in down_LV
    
        all_lv_nodes[lv_bus] += 1
    
        if any(in(list_MV).(nodes))  # Don't consider eventual presence of MV nodes
            list = nodes[in(list_MV).(nodes)]
            [all_lv_nodes[i] += 1 for i in nodes if i ∉ list]
        else
            [all_lv_nodes[i] += 1 for i in nodes]
        end
    
    end
    
    [println("Missing LV bus: ",j) for (j,i) in all_lv_nodes if i == 0]
    [print("ERROR bus :", i, "should not be present") for (i,j) in all_lv_nodes if i in list_MV && j>0]
    [print("ERROR bus :", i, "should not be present") for (i,j) in all_lv_nodes if i in list_MV]
    
    s = []
    [push!(s,i) for (j,i) in all_lv_nodes]
    s = sum(s)
    
    # If the two numbers are equal, it means that down_LV contains all the LV nodes with no double counting
    
    println("Number of LV buses: ", length(list_LV))
    println("Counted LV buses: ", s)
    =#
    
    # COMPUTE LOAD FOR EACH SECONDARY SUBSTATION
    
    loads = Dict(i =>[] for i in keys(down_MV))
    
    for (bus,lv_node) in down_MV
        
        lv_node = lv_node[1]
        Ptot = 0
        Qtot = 0
        
        # Add the power data of the LV node of the secondary substation
        for (ind, data) in net_data["load"]
            if data["load_bus"] == lv_node
                Ptot = data["pd"]
                Qtot = data["qd"]
                break
            end
        end
        
        # Add the remaining power profiles of the LV nodes connected to the same secondary substation
        for i in down_LV[lv_node]
            for (ind,data) in net_data["load"]
                if data["load_bus"] == i
                    Ptot += data["pd"]
                    Qtot += data["qd"]
                    break
                end
            end
        end
        
        loads[bus] = [Ptot,Qtot]
    end
    
    # Check if power coincide
    p = []
    [push!(p,i[1]) for (j,i) in loads]
    p = round(sum(p), digits = 5)
    
    println("\nP demand computed on LV= ",p, " kW")
    
    s = []
    [push!(s,data["pd"]) for (ind,data) in ref[:load] if data["load_bus"] ∉ list_MV]
    s = round(sum(s), digits = 5)
    
    println("P demand on LV side from dataset= ",s, " kW")

    t = []
    [push!(t,data["pd"]) for (ind,data) in ref[:load] if data["load_bus"] in list_MV]
    t = round(sum(t), digits = 5)
    
    println("P demand on MV side = ",t," kW")
    println("\nTotal power demand in the grid: ", t+s, " kW")
    
    # ADD DEMAND TO MV NODE 
    println("\nAdding demand profiles")

    for (bus, power) in loads
        
        index = length(keys(net_data["load"]))+1
        flag = true
        
        # Update load of MV node in secondary substation if node already present in the list of loads in the dataset
        for (ind, data) in net_data["load"]
    
            if data["load_bus"] == bus
                flag = false
                net_data["load"][ind]["pd"] = data["pd"] + power[1]
                net_data["load"][ind]["qd"] = data["qd"] + power[2]
                print("CIAOOOOOOOOOO \n\n")
                break

            end
        end
        

        if flag
            
            power[1] = round(power[1], digits = 5)
            power[2] = round(power[2], digits = 5)
            net_data["load"]["$index"] = Dict{String, Any}("source_id"=> ["bus",index+1], "load_bus" => bus, "status" => 1,"qd" => power[2],"pd" => power[1],"index" => index)
            
        end

    end
    return net_data,loads
end

function map_nodes(data::Dict)

    list_bus = []

    for (ind, bus) in data["bus"]
        push!(list_bus, bus["bus_i"])
    end
    
    list_bus = sort(list_bus)

    map = sort(Dict( true_ind => index for (index, true_ind) in enumerate(list_bus)))

    return map
end

function inverse_map(data::Dict)
    
    list_bus = []

    for (ind, bus) in data["bus"]
        push!(list_bus, bus["bus_i"])
    end
    
    list_bus = sort(list_bus)

    inv_map = sort(Dict( index => true_ind for (index, true_ind) in enumerate(list_bus)))

    return inv_map
end

function eliminate_LV_branches(case::Dict)
    
    dummy = deepcopy(case)

    list_LV = []

    for (ind, bus) in case["bus"]

        if bus["base_kv"] < 10  # save buses whose voltage is higher than 10 (i.e. MV)
            pop!(dummy["bus"], ind)
            push!(list_LV,parse(Int64,ind))
        end
    end

    for (ind,load) in case["load"]

        if load["load_bus"] in list_LV
            pop!(dummy["load"], ind)
        end
    end


    list_MV = []
    inactive = []
    [push!(list_MV,parse(Int64,i)) for i in keys(dummy["bus"])]   # save as a list of integers the Mv busses

    for i in keys(case["branch"])

        # I check if each edge of a branch has a MV bus 

        if case["branch"][i]["f_bus"] in list_MV && case["branch"][i]["t_bus"] in list_MV
            continue
        else
            pop!(dummy["branch"],i)
            pop!(dummy["distance"],i)
        end
    end

    # Don't consider inactive branches
   
    for i in keys(case["branch"])
        if case["branch"][i]["br_status"] == 0
            push!(inactive,i)
            pop!(dummy["branch"],i)
            pop!(dummy["distance"],i)
        end
    end

    [print("The inactive branches are indexed as: $i \n\n") for i in inactive]

    return dummy

end



# EXPORT MODIFIED NETOWRK

#export_matpower(replace(file, "matpower" => "Test"),net_data)

#powerplot(net_data, width = 800, height = 800, branch_size = 2, bus_size = 50, gen_size = 100)
#=
result = solve_ac_pf(net_data, Ipopt.Optimizer)
update_data!(net_data, result["solution"])

flows = calc_branch_flow_ac(net_data)
update_data!(net_data, result["solution"])

print_summary(result["solution"])
=#
#=
#print_summary(result["solution"]["flows"])
display(net_data)
print_summary(result["solution"])
=#
