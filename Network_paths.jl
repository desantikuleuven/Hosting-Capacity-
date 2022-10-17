using PowerModels
using PowerPlots
import Pkg
using FileIO
#Pkg.add("LightGraphs")
#Pkg.add("GraphPlot")
#Pkg.add("Graphs")

#=
    Computes all the nodes belonging to a path (lines) in a network. 
    Retrives also the longest path in terms of km. 
=#

using GraphPlot
using LightGraphs

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

function downstreamcalcs(feeder_data::Dict)

    ref = PowerModels.build_ref(feeder_data)[:it][:pm][:nw][0]
    arcs_from = sort(ref[:arcs_from])
    n_busses = length(keys(ref[:bus]))

    map = map_nodes(net_data)

    Gdir = LightGraphs.SimpleDiGraph(n_busses)

    # Add arcs in directed Graph

    for (l, i, j) in arcs_from
        i = map[i]
        j = map[j]
        LightGraphs.add_edge!(Gdir, i, j)   
    end

    # For each node, add all the downstream nodes (i.e. all possible path nodes)  

    downstream_bynode = sort(Dict(b => [LightGraphs.dst(e) for e in collect(LightGraphs.edges(LightGraphs.dfs_tree(Gdir, b; dir=:in)))] for b in 1:n_busses))
    
    return downstream_bynode, Gdir, map
end

function paths_calc(g, extremes, len, dist, tot_paths, inv_map, file, node_dist)

    for i in extremes

        l = 0
        dist_vector = [0.0]
        flag = true
        path = [(inv_map[i],inv_map[i])]
    
        while flag 
    
            j = outneighbors(g,i)[1]
            push!(path,(inv_map[i], inv_map[j])) 
    
            if j == outneighbors(g,1)[1]   #Tutte le strade portano a Roma  (j==45 for offical_rural)
                l = l + len[(i,j)]
                push!(dist, l)
                push!(tot_paths, path)
                push!(dist_vector, round(l, digits = 4))
                dist_vector = round.(maximum(dist_vector).-reverse(dist_vector), digits = 4)
                push!(node_dist, dist_vector)
                flag = false
            else
                l = l + len[(i,j)]
                push!(dist_vector, round(l, digits = 4))
                i = j
            end
        end
    end
    
    popfirst!(tot_paths)
    [dist[ind] = round(value, digits=3) for (ind,value) in enumerate(dist)]
    
    for (ind,vec) in enumerate(tot_paths)
        [tot_paths[ind][j] = reverse(tot_paths[ind][j]) for j in 1:length(vec)]
    end
    [tot_paths[i] = reverse(tot_paths[i]) for i in keys(tot_paths)]
    #=
    for i in 1:length(tot_paths)
        println("Path $i: ", tot_paths[i], " has a length of ", dist[i], " km \n")
    end

    # Alternative printing 
    
    for i in 1:length(tot_paths)
        print("Path $i: ")
        for (j,k) in tot_paths[i][2:end]
            print(" --> ",j, " --> ",k)
        end
        print("\n\n")
    end
    
    replace(".m"-->"distances.txt", file)
    open(replace(".m"-->"distances.txt", file), "w")
    =#
    
    max, ind = findmax(dist)
    print("max length is: ", max, " km correspdonding to path $ind\n")

    return tot_paths, dist, node_dist
end
#=

file = "Test_rural.m"

net_data = parse_file(file)

# INPUT FOR PATH CALC

down_node, g, map = downstreamcalcs(net_data)

extremes = []
[push!(extremes, i) for (i,j) in down_node if j ==[]]

tot_paths = [Tuple[]]
dist = AbstractFloat[]

len = Dict((map[j["from_b"]],map[j["to_b"]]) => j["length"] for (ind,j) in net_data["distance"])

inv_map = inverse_map(net_data)

node_dist = []

# RUN PATH CALC

paths_calc(g, extremes, len, dist, tot_paths, inv_map, file, node_dist)
=#