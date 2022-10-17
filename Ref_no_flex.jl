# POWER FLOW WITH DGs BUT NO FLEXIBILITY AND NO VOLT/CURRENT CONSTRAINTS 
var(pm::AbstractPowerModel, nw::Int=nw_id_default) = _IM.var(pm, pm_it_sym, nw)
var(pm::AbstractPowerModel, nw::Int, key::Symbol) = _IM.var(pm, pm_it_sym, nw, key)
var(pm::AbstractPowerModel, nw::Int, key::Symbol, idx) = _IM.var(pm, pm_it_sym, nw, key, idx)
var(pm::AbstractPowerModel, key::Symbol; nw::Int=nw_id_default) = _IM.var(pm, pm_it_sym, key; nw = nw)
var(pm::AbstractPowerModel, key::Symbol, idx; nw::Int=nw_id_default) = _IM.var(pm, pm_it_sym, key, idx; nw = nw)

function build_ref_no_flex(pm::AbstractPowerModel)
    
    variable_bus_voltage(pm, bounded = false)
    variable_gen_power(pm, bounded = false)
    variable_branch_power(pm, bounded = false)

    constraint_model_voltage(pm)  #do nothing

    for (i,bus) in ref(pm, :ref_buses)
        @assert bus["bus_type"] == 3

        constraint_theta_ref(pm, i)
        constraint_voltage_magnitude_setpoint(pm, i)

    end

    for (i,bus) in ref(pm, :bus)
        constraint_power_balance(pm, i)

        # PQ Bus Constraints
        if length(ref(pm, :bus_gens, i)) > 0 && !(i in ids(pm,:ref_buses))
            # this assumes inactive generators are filtered out of bus_gens
            @assert bus["bus_type"] == 2  # Should be ==1 but PowerModels is stupid

            for j in ref(pm, :bus_gens, i)
                constraint_gen_setpoint_active(pm, j)
                constraint_gen_setpoint_reactive(pm, j)
            end
        end
    end

    for i in ids(pm, :branch)
        constraint_ohms_yt_to(pm,i)
        constraint_ohms_yt_from(pm,i)

        #constraint_thermal_limit_from_me(pm, i)
        #constraint_thermal_limit_to_me(pm, i)
    end

end

function constraint_thermal_limit_from_me(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)

    if haskey(branch, "rate_a")
        p_fr = var(pm, nw, :p, f_idx)
        q_fr = var(pm, nw, :q, f_idx)
        
        JuMP.@constraint(pm.model, p_fr^2 + q_fr^2 <= (branch["cong_cap"]*branch["rate_a"])^2)  #80% loading 
    end
end

function constraint_thermal_limit_to_me(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    t_idx = (i, t_bus, f_bus)

    if haskey(branch, "rate_a")
        p_to = var(pm, nw, :p, t_idx)
        q_to = var(pm, nw, :q, t_idx)

        JuMP.@constraint(pm.model, p_to^2 + q_to^2 <= (branch["cong_cap"] * branch["rate_a"])^2)  #80% loading 
    end
end
