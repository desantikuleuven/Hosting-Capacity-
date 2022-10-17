#####################################################
################## COMMON FUNCTIONS #################
#####################################################


# Same function as calc_voltage_profile but without showing voltage profiles, needed just to get feeders info
function get_feeder_data(net_data::Dict, file_name)
    
    # INPUT FOR PATH CALC

    down_node, g, map = downstreamcalcs(net_data)

    extremes = []
    [push!(extremes, i) for (i,j) in down_node if j ==[]]

    tot_paths = [Tuple[]]
    dist = AbstractFloat[]

    len = Dict((map[j["from_b"]],map[j["to_b"]]) => j["length"] for (ind,j) in net_data["distance"])

    inv_map = inverse_map(net_data)
    node_dist = []

    tot_paths, dist, node_dist = paths_calc(g, extremes, len, dist, tot_paths, inv_map, file_name, node_dist)

    paths = []
    feeders = []

    mv_busbar = 0

    #COUNT NUMBEER OF FEEDERS
    for (path,line) in enumerate(tot_paths)
    
        if line[1][2]!=1 && line[1][2] ∉ feeders
            push!(feeders, line[1][2])
        end
    end

    # CREATE DICT CONTAINING ALL DATA RELEVANT FOR FEEDERS
    alphabet = 'A':'Z'
    feeder_ID = Dict( j => Dict("Name" => "Feeder $i", "Paths" => [], "Paths_ID" => [], "Paths_distance"=> []) for (i,j) in zip('A':alphabet[length(feeders)], feeders) )

    # Write tot_paths in a vectorial way. 
    for vec in tot_paths

        feed = []
        for j in 1:length(vec)

            push!(feed,vec[j][1])
            
        end
    
        push!(paths,feed)  #Rewrite how paths are made

    end

    for (idx,path) in enumerate(paths)

        if path[end]!=1 && path[2] in keys(feeder_ID)
            push!(feeder_ID[path[2]]["Paths"], path)
            push!(feeder_ID[path[2]]["Paths_ID"], idx)
            push!(feeder_ID[path[2]]["Paths_distance"], node_dist[idx])
        end

    end
    
    # Add ID of branches belonging to the same feeder

    for (ref, data) in feeder_ID

        branches = []
    
        for path in feeder_ID[ref]["Paths"]
            
            for (ind,j) in enumerate(path)
        
                if ind<length(path)
        
                    t_bus = path[ind]
                    f_bus = path[ind+1]
        
                    for (idx,branch) in net_data["branch"]
        
                        if f_bus == branch["f_bus"] && t_bus == branch["t_bus"]
                            if idx ∉ branches
                                push!(branches,idx)
                            end
                        end
                    end
                end
            end
        end
    
        feeder_ID[ref]["Branches"] = branches
    
    end
    
    # Add ID of bus belonging to the same feeder
    for (ref, data) in feeder_ID

        nodes = []

        for path in feeder_ID[ref]["Paths"]
            [push!(nodes, bus) for bus in path if bus ∉ nodes]
        end

        feeder_ID[ref]["Buses"] = nodes
        mv_busbar = nodes[1]
    end

    return feeder_ID, mv_busbar
end

#=
gen_number_per_feeder --> define how many DGs you want in each feeder
This number will define how many buses will be picked to become active in the grid, chosen randomly. 
Gives back dictionary with the buses selected for DGs installation    
=#
function get_random_generators(feeder_ID_1, gen_number_per_feeder, seed_value::Int, seed = true)

    generators = Dict()
    if seed
        Random.seed!(seed_value)
    end
    for (id, feeder) in feeder_ID_1

        if length(feeder["Buses"])>gen_number_per_feeder
            active_nodes = sample(feeder["Buses"][2:end], gen_number_per_feeder, replace = false)  #Choose n random buses in each feeder, each one with same size determined by p_target/n
            generators[feeder["Name"]] = active_nodes
        else
            active_nodes = sample(feeder["Buses"][2:end], length(feeder["Buses"])-2, replace = false)  # Set number of generators equal to number of nodes - 2
            generators[feeder["Name"]] = active_nodes
        end
    end

    return generators
end

#=
Each generator will get installed a power equal to size_std. 
These generators are added to the model. 
The list of generators are passed as a Dict, where you give a list of buses for each feeders. This is done by calling before get_random_generators
=#
function add_generators(net_data::Dict, generators::Dict{Any, Any}, size_std)
    
    x = []
    [append!(x,v) for v in values(generators)]

    for (i,gen) in enumerate(x)
        i+=1
        net_data["gen"]["$i"] = Dict("pg" =>size_std, "qg" =>0, "pmin" => size_std , "pmax"=> size_std , "qmin" =>0, "qmax"=>0, "gen_bus" => gen, "gen_status"=>1, "index" => i, "source_id" => ["gen", i])
        net_data["bus"]["$gen"]["bus_type"] = 2
    end
end

function get_gen_info(net_data::Dict, feeder_ID::Dict)

    generator_ID = Dict()

    for (j, feeder) in feeder_ID

        for (i,gen) in net_data["gen"]
            
            present = false

            for k in feeder["Paths"]
                if gen["gen_bus"] in k
                    present = true
                end
            end

            if present && gen["gen_bus"]!="1"  #exclude slack generator
                generator_ID[i] = Dict("ref" => j, "feeder"=> feeder["Name"], "bus"=>gen["gen_bus"],"pmin" => gen["pmin"], "p_nominal" => gen["pg"], "q_nominal"=> gen["qg"])
            end
        end
    end


    return generator_ID

end

# Function used to increase exisiting capacity of generators in a specific feeder. Each generator gets increased by "additiona_power" value
function gen_increase_size(net_data, gen_ID, feeder, additional_power)

    for feeder_name in [feeder]
        for (id,gen) in gen_ID

            if gen["feeder"] == feeder_name
                net_data["gen"][id]["pg"] += additional_power
                net_data["gen"][id]["pmax"] += additional_power
                net_data["gen"][id]["pmin"] += additional_power
                gen_ID[id]["p_nominal"] += additional_power
            end
        end
    end

end

# Update net_data generation values based on data stored on gen_ID
function update_gens(net_data, gen_ID)

    for (id, gen) in gen_ID
        net_data["gen"][id]["pg"] = gen["p_nominal"]
    end

end


##################################################################
################## FUNCTIONS FOR HC WITHOUT FLEX #################
##################################################################



# Function to check if voltage violation occurrred in the feeder
function check_voltage_constraints(feeder_ID, healthy_feeders)
    feeder_voltage_constraint = Dict()
    
   
    if !isempty(healthy_feeders)

        [feeder_voltage_constraint[i] = false for i in healthy_feeders]
        
        for (id, data) in feeder_ID
            if data["Name"] in healthy_feeders
                v_max = data["vmax"]
                v_min = data["vmin"]

                if v_min < 0.9
                    
                    feeder_voltage_constraint[data["Name"]] = true 

                end

                if v_max > 1.1

                    feeder_voltage_constraint[data["Name"]] = true 
                    
                end
            end

        end
    end

    return feeder_voltage_constraint
end

function check_congestion_constraints(feeder_ID, healthy_feeders)
    feeder_congestion_constraint = Dict()

    [feeder_congestion_constraint[i] = false for i in healthy_feeders]
    
    if !isempty(healthy_feeders)
        for (id, data) in feeder_ID
            if data["Name"] in healthy_feeders
                if !isempty(data["Congested_lines"])
                    feeder_congestion_constraint[data["Name"]] = true
                end
            end
        end
    end

    return feeder_congestion_constraint
end

# Adds to each feeder which kind of violation occurred
function add_violation_type(feeder_HC, volt_violation, congestion_violation)

    if volt_violation && congestion_violation
        feeder_HC["Volt_violation"] = feeder_HC["Congestion_violation"] = true
        print(" Limited by both Thermal and Voltage constraints\n")
    elseif volt_violation
        feeder_HC["Volt_violation"]= true
        feeder_HC["Congestion_violation"] = false
        print(" Limited by Voltage bounds\n")
    elseif congestion_violation
        feeder_HC["Congestion_violation"] = true
        feeder_HC["Volt_violation"] = false
        print(" Limited by Thermal constraints\n")
    end

end

#=
    This function computes the HC for each feeder. An ACPF is performed with DGs present in each feeder.
    For each PF, a network constraint analysis is done to each feeder to check any voltage or congestion violation. 
    If no violation occurs, the generation capacity of the feeder is increase by "power_target_per_feeder".
    Return the HC of each feeder and the new generation values. 
    
    Ref: build_ref_no_flex does not have any voltage or congestion constraint. DGs defined as PQ buses. 

    ! No reactive injection taken into acocunt !

=#
function calc_hosting_capacity(net_data, ref, HC, power_target_per_feeder, threshold; feeder_ID_1, gen_ID, file_name, gen_number_per_feeder )
    feeder_HC = Dict()
    gen_step = power_target_per_feeder/gen_number_per_feeder
    
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

        pm = instantiate_model(net_data, ACPPowerModel, ref)
        global result = optimize_model!(pm, optimizer=JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0))
        update_data!(net_data, result["solution"])
        println("\n Termination status: ", result["termination_status"])

        # get info on voltage and loading profiles for each feeder
        global feeder_ID, path_volt, mv_busbar = calc_voltage_profile(net_data, result, file_name)
        calc_branch_loading(net_data, feeder_ID, gen_ID, threshold)

        # check network violations
        feeder_voltage_constraint = check_voltage_constraints(feeder_ID, healthy_feeders)
        feeder_congestion_constraint = check_congestion_constraints(feeder_ID, healthy_feeders)

        if !isempty(healthy_feeders) #healthy_feeders gets updated as soon as a feeder violates constraints
            stable = true
            feeders_iterator = copy(healthy_feeders)

            for f_name in feeders_iterator

                volt_violation = feeder_voltage_constraint[f_name]
                congestion_violation = feeder_congestion_constraint[f_name]

                if !volt_violation && !congestion_violation

                    gen_increase_size(net_data, gen_ID, f_name, gen_step)
                    println("\n Increasing capacity generation at $f_name of $(power_target_per_feeder) MW")
                    HC += power_target_per_feeder
                else
                    #Reduce Pg of each gen to the previous safe generation value
                    feeder_HC[f_name] = Dict("HC"=> sum([gen["p_nominal"] - gen_step for (id, gen) in gen_ID if gen["feeder"] == f_name]))
                    println("\n Network constraint violated at $f_name for ",feeder_HC[f_name]["HC"], " MW installed")

                    for (id, gen) in gen_ID
                        if gen["feeder"] == f_name
                            gen_ID[id]["p_nominal"] -= gen_step
                        end
                    end
                    add_violation_type(feeder_HC[f_name], volt_violation, congestion_violation)
                    update_gens(net_data, gen_ID)
                    deleteat!(healthy_feeders, findall(x->x==f_name, healthy_feeders))
                end
            end

        else
            HC = sum([feeder["HC"] for (id,feeder) in feeder_HC])
            print("\n Hosting Capacity evaluation completed!\n Maximum Network HC: $HC MW \n")
        end


    end
    
    return feeder_HC, gen_ID

end


##################################################################
################## FUNCTIONS FOR HC WITH FLEX ####################
##################################################################



#=
    This function computes the HC for each feeder. An ACPF is performed with DGs present in each feeder.
    For each PF, a network constraint analysis is done to each feeder to check any voltage or congestion violation. 
    If no violation occurs, the generation capacity of the feeder is increased by "power_target_per_feeder".
    Return the HC of each feeder and the new generation values. 

    If network constraints are violated, a new PF is performed to check if problems can be solved by exploiting flexibility.
    DGs defined as PQ buses. 

    ! No reactive injection taken into acocunt !

=#

function calc_hosting_capacity_w_flex(net_data, ref, HC, gen_step; feeder_ID_1, gen_ID, gen_number_per_feeder )
    feeder_HC = Dict()
    gen_step
    
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

        global result = solve_branch_voltage_cuts_HC(net_data, ACPPowerModel, Ipopt.Optimizer, ref,gen_ID, feeder_ID_1, feeder_HC, gen_step, healthy_feeders, feeder_names)
        
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
    
    return feeder_HC, gen_ID

end

function solve_branch_voltage_cuts_HC(data::Dict{String,<:Any}, model_type::Type, optimizer,ref,gen_ID, feeder_ID_1, feeder_HC, gen_step,healthy_feeders, feeder_names; kwargs...)
    return solve_branch_voltage_cuts_HC!(data, model_type, optimizer, ref, gen_ID, feeder_ID_1, feeder_HC, gen_step, healthy_feeders, feeder_names; kwargs...)
end

function solve_branch_voltage_cuts_HC!(data::Dict{String,<:Any}, model_type::Type, optimizer, ref,gen_ID, feeder_ID_1, feeder_HC, gen_step, healthy_feeders, feeder_names; solution_processors=[], max_iter::Int=100, time_limit::Float64=3600.0)

    Memento.info(_LOGGER, "maximum cut iterations set to value of $max_iter")

    for (i,branch) in data["branch"]
        if haskey(branch, "rate_a")
            branch["rate_a_inactive"] = branch["rate_a"]
        end
    end

    start_time = time()

    pm = instantiate_model(data, model_type, ref)
    result = optimize_model!(pm, optimizer=JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0), solution_processors=solution_processors)

    iteration = 1
    violated = true
    while violated && iteration < max_iter && (time() - start_time) < time_limit 
        violated = false
       
        # Congestion constraint
        for (i,branch) in data["branch"]
            
            if haskey(branch, "rate_a_inactive")
                rate_a = branch["rate_a_inactive"]
                branch_sol = result["solution"]["branch"][i]

                cap = branch["cong_cap"]

                mva_fr = abs(branch_sol["pf"])
                mva_to = abs(branch_sol["pt"])

                if !isnan(branch_sol["qf"]) && !isnan(branch_sol["qt"])
                    mva_fr = sqrt(branch_sol["pf"]^2 + branch_sol["qf"]^2)
                    mva_to = sqrt(branch_sol["pt"]^2 + branch_sol["qt"]^2)
                end

                #println(branch["index"], rate_a, mva_fr, mva_to)

                if mva_fr > rate_a * cap || mva_to > rate_a * cap
                    Memento.info(_LOGGER, "activate rate_a * cap on branch $(branch["index"])")

                    branch["rate_a"] = branch["rate_a_inactive"]
                    delete!(branch, "rate_a_inactive")

                    idx = branch["index"]
                    constraint_thermal_limit_from_me(pm, idx)
                    constraint_thermal_limit_to_me(pm, idx)

                    data["branch"][i]["violated"] = true 
                    data["branch"][i]["violation_type"] = "thermal"

                    violated = true
                else
                    data["branch"][i]["violated"] = false
                end
            end
            
        end

        # Voltage constraint 
        for (i,bus) in data["bus"]

            
            vm = round(result["solution"]["bus"][i]["vm"],digits = 5)
            v_min = bus["vmin"]
            v_max = bus["vmax"]
            idx = bus["index"]
            if vm < v_min
    
                Memento.info(_LOGGER, "activate v_min on bus $i")
                constraint_voltage_magnitude_lower_me(pm,idx)
                data["bus"][i]["violated"] = true 
                data["bus"][i]["violation_type"] = "down_v"
                
                violated = true

            elseif vm > v_max

                Memento.info(_LOGGER, "activate v_max on bus $i")
                constraint_voltage_magnitude_upper_me(pm,idx)
                data["bus"][i]["violated"] = true 
                data["bus"][i]["violation_type"] = "up_v"

                violated = true
                
            else
                data["bus"][i]["violated"] = false
            end 
            

        end

        violated_bus = Vector()
        violated_branch = Vector()
        [push!(violated_branch,id) for (id, branch) in data["branch"] if branch["violated"]]
        [push!(violated_bus,id) for (id, bus) in data["bus"] if bus["violated"]]

        violated_feeders = find_violated_feeders(violated_bus,violated_branch, feeder_ID_1)

        if violated

            println("Trying to solve network constraints by exploiting flexibility")
            iteration += 1
            result = optimize_model!(pm, solution_processors=solution_processors)
            println(result["termination_status"])

            if result["termination_status"] != LOCALLY_SOLVED  #Flexibility insufficient to solve network problems on all feeders

                violated_feeders = violated_feeders[in(healthy_feeders).(violated_feeders)]  # Find which feeder among the healthy feeders present network issues

                if length(violated_feeders)>1  # try to see if we can get an optimal solution by decreasing generation feeder by feeder
                    
                    print("\n Could not globally solve network issues at: ")
                    [println(f_name) for f_name in violated_feeders]
                  

                    for (id, gen) in gen_ID  #reducing capacity of violated feeders to previous safe operation 
                        if gen["feeder"] in violated_feeders 
                            gen_ID[id]["p_nominal"] -= gen_step
                        end
                    end
                    update_gens(data, gen_ID)

                    for f_name in violated_feeders
                        
                        println("\n Trying to check feasibility when increasing generation only at $f_name \n")
                    
                        for (id, gen) in gen_ID
                            if gen["feeder"] == f_name
                                gen_ID[id]["p_nominal"] += gen_step
                            end
                        end
                        update_gens(data, gen_ID)

                        println(" Capacity installed at unhealthy feeders: ", )
                        [println(" - $feeder ", sum([gen["p_nominal"] for (id, gen) in gen_ID if gen["feeder"] == feeder]), " MW") for feeder in sort(violated_feeders)]
                        
                        println("\n Running PF")
                        update_gens(data, gen_ID)
                        result = run_pf_simplified(data, ACPPowerModel, ref) #Run basically solve_branch_voltage_cuts_HC! without all these if statements

                        if result["termination_status"] != LOCALLY_SOLVED  # flex insufficient to save f_name
                            
                            for (id, gen) in gen_ID  #safe capacity for feeder f_name was the previous size
                                if gen["feeder"] == f_name
                                    gen_ID[id]["p_nominal"] -= gen_step
                                end
                            end
                            update_gens(data, gen_ID)

                            feeder_HC[f_name] = Dict("HC"=> sum([gen["p_nominal"] for (id, gen) in gen_ID if gen["feeder"] == f_name]))
                            println("\n Hosting capacity for $f_name reached at ",feeder_HC[f_name]["HC"], " MW installed")

                            deleteat!(healthy_feeders, findall(x->x==f_name, healthy_feeders))
                            
                        else
                            println("Flexibility can be exploited to solve network issues at $f_name")
                            violated = false
                            
                        end

                    end

                    if violated #if flexibility cannot solve the network issues for none of the violated feeders
                        #=
                        for f_name in violated_feeders

                            #Reduce Pg of each gen to the previous safe generation value
                            feeder_HC[f_name] = Dict("HC"=> sum([gen["pg"] - gen_step for (id, gen) in gen_ID if gen["feeder"] == f_name]))
                            println("\n Hosting capacity for $f_name reached at ",feeder_HC[f_name]["HC"], " MW installed")
            
                            for (id, gen) in gen_ID
                                if gen["feeder"] == f_name
                                    gen_ID[id]["pg"] -= gen_step
                                end
                            end
                            
                            deleteat!(healthy_feeders, findall(x->x==f_name, healthy_feeders))
                        
                        end
                    
                        update_gens(data, gen_ID)
                        =#

                        println("Running a check PF with: ")
                        [println(f_name," ", round(sum([gen["p_nominal"] for (id, gen) in gen_ID if gen["feeder"] == f_name]), digits = 4), " MW") for f_name in feeder_names]
                        
                        # in this case the solution should be feasible without going in the if statements
                        result = run_pf_simplified(data, ACPPowerModel, ref)
                        #result = solve_branch_voltage_cuts_HC(data, ACPPowerModel, Ipopt.Optimizer, ref, gen_ID, feeder_ID_1, feeder_HC, gen_step, healthy_feeders)
                        update_data!(data, result) 
                        println(result["termination_status"])

                    end
                else
                    
                    f_name = violated_feeders[1]
                    feeder_HC[f_name] = Dict("HC"=> sum([gen["p_nominal"] - gen_step for (id, gen) in gen_ID if gen["feeder"] == f_name]))
                    println("\n Hosting capacity for $f_name reached at ",feeder_HC[f_name]["HC"], " MW installed")

                    for (id, gen) in gen_ID
                        if gen["feeder"] == f_name
                            gen_ID[id]["p_nominal"] -= gen_step
                        end
                    end
                    update_gens(data, gen_ID)
                    deleteat!(healthy_feeders, findall(x->x==f_name, healthy_feeders))
                    
                    println("Running a check PF with: ")
                    [println(f_name," ", round(sum([gen["p_nominal"] for (id, gen) in gen_ID if gen["feeder"] == f_name]), digits = 4), " MW") for f_name in feeder_names]
  
                    #result = solve_branch_voltage_cuts_HC(data, ACPPowerModel, Ipopt.Optimizer, ref, gen_ID, feeder_ID_1, feeder_HC, gen_step, healthy_feeders)
                    result = run_pf_simplified(data, ACPPowerModel, ref)
                    update_data!(data, result) 
                    
                    if result["termination_status"] == LOCALLY_SOLVED
                        violated=false
                    end

                end

            else
                println("Network violations successuflly solved with flexibility at: \n")
                [println("- ",f_name) for f_name in violated_feeders] 
                violated = false
            end
        else
            Memento.info(_LOGGER, "flow cuts converged in $iteration iterations")
        end
    end
    

    result["solve_time"] = time() - start_time
    result["iterations"] = iteration

    return result
end

# At the end of the unconstrained PF, identify which feeders show network violations 
function find_violated_feeders(violated_bus,violated_branch, feeder_ID_1)
    violated_feeders = Vector()

    for (id,feeder) in feeder_ID_1
        for bus in violated_bus
            bus = parse(Int64,bus)
            if bus in feeder["Buses"][2:end] && feeder["Name"] ∉ violated_feeders
                push!(violated_feeders, feeder["Name"])
            end
        end
        for branch in violated_branch
            if branch in feeder["Branches"] && feeder["Name"] ∉ violated_feeders
                push!(violated_feeders, feeder["Name"])
            end
        end

    end
    return violated_feeders
end

# One iteration of the unconstrained PF, used in solve_branch_voltage_cuts_HC as a feasibility check PF.
function run_pf_simplified(data::Dict{String,<:Any}, model_type::Type, ref; solution_processors=[], max_iter::Int=100, time_limit::Float64=3600.0)

    Memento.info(_LOGGER, "maximum cut iterations set to value of $max_iter")

    for (i,branch) in data["branch"]
        if haskey(branch, "rate_a")
            branch["rate_a_inactive"] = branch["rate_a"]
        end
    end

    pm = instantiate_model(data, model_type, ref)
    result = optimize_model!(pm, optimizer=JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0), solution_processors=solution_processors)
    
    violated = false
    
    # Congestion constraint
    for (i,branch) in data["branch"]
        if haskey(branch, "rate_a_inactive")
            rate_a = branch["rate_a_inactive"]
            branch_sol = result["solution"]["branch"][i]

            cap = branch["cong_cap"]

            mva_fr = abs(branch_sol["pf"])
            mva_to = abs(branch_sol["pt"])

            if !isnan(branch_sol["qf"]) && !isnan(branch_sol["qt"])
                mva_fr = sqrt(branch_sol["pf"]^2 + branch_sol["qf"]^2)
                mva_to = sqrt(branch_sol["pt"]^2 + branch_sol["qt"]^2)
            end

            #println(branch["index"], rate_a, mva_fr, mva_to)

            if mva_fr > rate_a * cap || mva_to > rate_a * cap
                Memento.info(_LOGGER, "activate rate_a * cap on branch $(branch["index"])")

                branch["rate_a"] = branch["rate_a_inactive"]
                delete!(branch, "rate_a_inactive")

                idx = branch["index"]
                constraint_thermal_limit_from_me(pm, idx)
                constraint_thermal_limit_to_me(pm, idx)

                data["branch"][i]["violated"] = true 
                data["branch"][i]["violation_type"] = "thermal"

                violated = true
            else
                data["branch"][i]["violated"] = false
            end
        end
    end

    # Voltage constraint 
    for (i,bus) in data["bus"]

        vm = round(result["solution"]["bus"][i]["vm"],digits = 5)
        v_min = bus["vmin"]
        v_max = bus["vmax"]
        idx = bus["index"]
        if vm < v_min

            Memento.info(_LOGGER, "activate v_min on bus $i")
            constraint_voltage_magnitude_lower_me(pm,idx)
            data["bus"][i]["violated"] = true 
            data["bus"][i]["violation_type"] = "down_v"
            
            violated = true

        elseif vm > v_max

            Memento.info(_LOGGER, "activate v_max on bus $i")
            constraint_voltage_magnitude_upper_me(pm,idx)
            data["bus"][i]["violated"] = true 
            data["bus"][i]["violation_type"] = "up_v"

            violated = true
            
        else
            data["bus"][i]["violated"] = false
        end 

    end

    violated_bus = Vector()
    violated_branch = Vector()
    [push!(violated_branch,id) for (id, branch) in data["branch"] if branch["violated"]]
    [push!(violated_bus,id) for (id, bus) in data["bus"] if bus["violated"]]

    violated_feeders = find_violated_feeders(violated_bus,violated_branch, feeder_ID_1)

    if violated

        println("Trying to solve network constraints by exploiting flexibility")
        result = optimize_model!(pm, solution_processors=solution_processors)
        if result["termination_status"] == LOCALLY_SOLVED
            println("Network constraints solved by contracting ", round(result["objective"], digits = 5), " MW of flexibility in the following feeders: \n")
            [println(f_name, "\n") for f_name in sort(violated_feeders)]
        else
            println("--> Unfeasible solution")
        end
    end

    return result
end

# Used at the end to evaluate which kind of violation has occurrred in all feeders
function add_violation_type_flex(net_data, feeder_ID, feeder_HC)

    [feeder_HC[f_name]["Congestion"] = false for (f_name,x) in feeder_HC]
    [feeder_HC[f_name]["V_up"] = false for (f_name,x) in feeder_HC]
    [feeder_HC[f_name]["V_down"] = false for (f_name,x) in feeder_HC]
    
    for (id,branch) in net_data["branch"] 
        count = []
        if haskey(branch,"violation_type")
            for (x,feed) in feeder_ID
                if id in feed["Branches"] && feed["Name"] ∉ count
                    feeder_HC[feed["Name"]]["Congestion"] = true
                    push!(count, feed["Name"])
                end
            end
        end
    end
    #=
    for (id,bus) in net_data["bus"] 
        count = []
        if haskey(bus,"violation_type")
            type = bus["violation_type"]
            for (x,feed) in feeder_ID
                if id in parse.(Int64,feed["Buses"]) && feed["Name"] ∉ count
                    push!(count, feed["Name"])
                    if type == "up_v"
                        feeder_HC[feed["Name"]]["V_up"] = true
                    elseif type == "down_v"
                        feeder_HC[feed["Name"]]["V_down"] = true
                    end
                end
            end
        end
    end
    =#
    for (id,bus) in net_data["bus"] 
        count = []
        id = parse(Int64,id)
        if haskey(bus,"violation_type")
            type = bus["violation_type"]
            for (x,feed) in feeder_ID
                if id in feed["Buses"] && feed["Name"] ∉ count
                    push!(count, feed["Name"])
                    if type == "up_v"
                        feeder_HC[feed["Name"]]["V_up"] = true
                    elseif type == "down_v"
                        feeder_HC[feed["Name"]]["V_down"] = true
                    end
                end
            end
        end
    end

end

function printing_statements_HC_flex(feeder_HC, feeder_ID, flex_loads, p_load, q_load)
    for (x, feeder) in feeder_ID
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
    
    
    
    # Limiting factors for feeders HC and flex offered
    for (id,feeder) in sort(feeder_HC)
        println("\n $id HC: ", round(feeder["HC"], digits = 4), " MW")
        
    
        if Bool(feeder["Congestion"]) && Bool(feeder["V_up"])
            println(" Limiting factor for $id: voltage upper bound and thermal limits\n")
    
        elseif !Bool(feeder["V_down"]) && !Bool(feeder["V_up"])
            println(" Limiting factor for $id: thermal limits\n")
        elseif !Bool(feeder["Congestion"]) && Bool(feeder["V_up"]) || !Bool(feeder["Congestion"]) && Bool(feeder["V_down"])
            println(" Limiting factor for $id: voltage bounds \n")
        end
    
        println(" Upwards Active Flexibility contracted: ", round(feeder["up_p_flex"], digits = 4), " MW")
        println(" Upwards Reactive Flexibility contracted: ", round(feeder["up_q_flex"], digits = 4), " MVar")
        println(" Downwards Active Flexibility contracted: ", round(feeder["down_p_flex"], digits = 4), " MW")
        println(" Downwards Reactive Flexibility contracted: ", round(feeder["down_q_flex"], digits = 4) , " MVar")
    end   
    
        
    
    print("\n Slack bus power generated: ", result["solution"]["gen"]["1"]["pg"], "MW \n")
    print("\n Immaginary power generated: ", result["solution"]["gen"]["1"]["qg"], "MVar \n")
    
    print("\n MV busbar voltage : ", round(result["solution"]["bus"]["$mv_busbar"]["vm"],digits =5), "p.u. \n")
    
    print("\n Real power requested: ", p_load, "MW \n")
    print("\n Immaginary power requested: ", q_load, "MVar \n")
    
    print("\n Real power losses: ", sum(p_loss), "MW \n")
    print("\n Immaginary power losses: ", sum(q_loss), "MVar \n")
    
    
    #=

    #[println("\n Vmin: ", round(feeder_ID[gen["ref"]]["vmin"], digits = 5)," Vmax: ",round(feeder_ID[gen["ref"]]["vmax"],digits = 5)," for ",gen["feeder"]," (where generator $i is connected)") for (i,gen) in gen_ID]
    #[println("\n Feeder loading: ", round(gen["max_branch_load"],digits = 3),"% in branch ",gen["Critical_branch"]) for (i,gen) in gen_ID]



    for (id, gen) in sort(gen_ID)
        println(" Generator $id is connected at bus ", gen["bus"], " in ", gen["feeder"], " with ",gen["pg"],"MW installed")
    end

    println()

    feeder_names = Vector{String}()
    [push!(feeder_names, feeder["Name"]) for (i,feeder) in feeder_ID]
    sort!(feeder_names)

    for f_name in feeder_names

        for (ref,data) in feeder_ID
            if data["Name"] == f_name

                println(" ",f_name," --> Vmin: ", round(data["vmin"], digits = 5), " Vmax: ", round(data["vmax"], digits = 5))
                
            end
        end
    end

    println()


    vmin = Vector()
    vmax = Vector()
    [push!(vmin,feeder["vmin"]) for (id, feeder) in feeder_ID ]
    println("\n Lowest voltage magnitude in the grid: ", minimum(vmin))
    [push!(vmin,feeder["vmax"]) for (id, feeder) in feeder_ID ]
    println("\n Highest voltage magnitude in the grid: ", maximum(vmin))

    =#


end

function create_df_HC_flex(feeder_HC, save = false)
    
    df = DataFrame()

    feeder_HC = sort(feeder_HC)

    HC = Vector()
    [push!(HC, data["HC"]) for (x,data) in feeder_HC]

    df.HC_flex = HC

    up_p_flx = Vector()
    [push!(up_p_flx, data["up_p_flex"]) for (x,data) in feeder_HC]
    up_q_flx = Vector()
    [push!(up_q_flx, data["up_q_flex"]) for (x,data) in feeder_HC]
    dwn_p_flx = Vector()
    [push!(dwn_p_flx, data["down_p_flex"]) for (x,data) in feeder_HC]
    dwn_q_flx = Vector()
    [push!(dwn_q_flx, data["down_q_flex"]) for (x,data) in feeder_HC]

    df.up_p_flex = round.(up_p_flx, digits = 4)
    df.up_q_flex = round.(up_q_flx, digits = 4)
    df.dwn_p_flex = round.(dwn_p_flx, digits = 4)
    df.dwn_q_flex = round.(dwn_q_flx, digits = 4)

    vio = Vector()
    
    for (f_name, data) in feeder_HC
        if Bool(data["V_up"]) && Bool(data["Congestion"])
            push!(vio, "Thermal & Voltage")
        elseif  Bool(data["V_up"]) && !Bool(data["Congestion"]) || Bool(data["V_down"]) && !Bool(data["Congestion"]) 
            push!(vio, "Voltage")
        else 
            push!(vio, "Thermal")
        end
    end

    df.violation = vio

    if save


        XLSX.openxlsx("C:\\Users\\u0152683\\Desktop\\Networks\\Hosting Capacity\\Results\\Results.xlsx", mode="rw") do xf
            sheet = xf["DR"]

            column_names = names(df)

            for c in 1:length(column_names)
                sheet[XLSX.CellRef(1,c)] = column_names[c]
            end

            for r in 1:size(df,1), c in 1:size(df,2)
                sheet[XLSX.CellRef(r+1,c)] = df[r,c]
            end
        end
    end
    return df
end

function create_df_HC(feeder_HC, save = false)
    
    df = DataFrame()

    HC = Vector()
    [push!(HC, data["HC"]) for (x,data) in feeder_HC]

    df.HC_flex = HC

    vio = Vector()
    
    for (f_name, data) in feeder_HC
        if Bool(data["Volt_violation"]) && Bool(data["Congestion_violation"])
            push!(vio, "Thermal & Voltage")
        elseif  Bool(data["Volt_violation"]) && !Bool(data["Congestion_violation"]) 
            push!(vio, "Voltage")
        else 
            push!(vio, "Thermal")
        end
    end

    df.violation = vio
    
    return df
end