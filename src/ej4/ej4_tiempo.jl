include("HH.jl")
using .HH

using Plots
using ArgParse
using JSON
using DifferentialEquations


function fill_defaults(cfg::Dict{String,Any})
    return Dict(
        "V0" => get(cfg, "V0", default_V0),
        "I" => get(cfg, "I", default_I),
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

function plot_solution(sol, fig_path, initial_conditions)
    initial_conditions_text = join(["$(k) = $(v)" for (k, v) in pairs(initial_conditions)], ", ")
    plt = plot(
        layout=(4, 1),
        size=(400, 800),
        xlabel="Time (ms)",
        legend=false,
        title="Initial conditions: $initial_conditions_text",
    )

    ylabels = ["V(mV)", "m", "h", "n"]

    for i in 1:4
        plot!(plt, sol.t, sol[i, :],
            subplot=i,
            ylabel=ylabels[i]
        )
    end

    savefig(plt, fig_path)
    gui()
end

args = parse_commandline()
config_path = args["config"]
fig_path = args["figpath"]

config = read_config(config_path)
tmax = get(config, "tmax", default_tmax)

for (i, plot_cfg) in enumerate(config["plots"])
    cfg = fill_defaults(plot_cfg)

    HH_model = HH_models[Symbol(cfg["model"])]

    V0 = cfg["V0"]
    I = cfg["I"]
    k = cfg["k"]
    g = Conductances()
    v = InversionV()

    sol = solve_system(HH_model, V0, I, g, v, C, k, tmax)
    plot_solution(sol, fig_path * string(i) * ".pdf", (; V0, I))
end

gui()

println("Press the enter key to quit:")
readline()