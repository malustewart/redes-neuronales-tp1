include("HH.jl")
using .HH

using Plots
using ArgParse
using JSON
using DifferentialEquations
using Statistics

const default_length = 100

function fill_defaults(cfg::Dict{String,Any})
    return Dict(
        "V0" => get(cfg, "V0", default_V0),
        "model" => get(cfg, "model", :HH_standard),
        "k" => get(cfg, "k", default_k)
    )
end

function read_config(filename::String)
    open(filename, "r") do io
        return JSON.parse(io)
    end
end

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--config", "-c"
        help = "path to config file"
        arg_type = String
        default = "config.json"
        "--figpath", "-f"
        help = "path where to save the resulting image"
        arg_type = String
        default = "imgs/lastfig"
    end
    return parse_args(s)
end

function solve_system(HH_model, V0, I, g, v, C, k, tmax)
    (; m_inf, h_inf, n_inf) = MhnParameters(V0)
    Y0 = [V0, m_inf, h_inf, n_inf]
    params = HHParams(g, v, I, C, k)
    tspan = (0, tmax)

    prob = ODEProblem(HH_model, Y0, tspan, params)
    solve(prob, Rosenbrock23())
end

function plot_solution(I, frequencies, fig_path, initial_conditions)
    initial_conditions_text = join(["$(k) = $(v)" for (k, v) in pairs(initial_conditions)], ", ")
    plt = plot(
        xlabel="I(uA)",
        ylabel="f(Hz)",
        legend=false,
        title="Initial conditions: $initial_conditions_text",
    )

    plot!(plt, I, frequencies)

    savefig(plt, fig_path)
    gui()
end

function get_frequency(time_solution)
    trigger = mean(time_solution.u[div(end,2)+1:end][1])

    edges_indices = filter(i -> time_solution[i][1] < trigger && time_solution[i+1][1] >= trigger, 1:(length(time_solution) - 1))  # -1 to skip last iteration
    
    if length(edges_indices) < 6 #assume system did not oscillate permanently, only in transitory
        return 0
    end

    edges_times = map(i -> time_solution.t[i], edges_indices)./1000 #convert from ms to s
    
    1.0./mean(diff(edges_times))
end

args = parse_commandline()
config_path = args["config"]
fig_path = args["figpath"]

config = read_config(config_path)
tmax = get(config, "tmax", default_tmax)
I_min = get(config, "I_min", default_I_min)
I_max = get(config, "I_max", default_I_max)
len = get(config, "I_max", default_length)

for (i, plot_cfg) in enumerate(config["plots"])
    cfg = fill_defaults(plot_cfg)

    HH_model = HH_models[Symbol(cfg["model"])]

    V0 = cfg["V0"]
    I = range(I_min, I_max, length = len)
    k = cfg["k"]
    g = Conductances()
    v = InversionV()

    time_solutions = map(i -> solve_system(HH_model, V0, i, g, v, C, k, tmax), I)
    frequencies = map(solution -> get_frequency(solution), time_solutions)

    plot_solution(I, frequencies, fig_path * "_fi_" * string(i) * ".pdf", (; V0, k))
end

println("Press the enter key to quit:")
readline()