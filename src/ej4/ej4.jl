using Plots
using ArgParse
using JSON
using DifferentialEquations

default_V0 = -80
default_I = 1
default_tmax = 100
g_Na = 120
g_K = 36
g_L = 0.3
V_Na = 50
V_K = -77
V_L = -54.4
C = 1e-6

struct Conductances
    g_Na
    g_K
    g_L
end

function Conductances()
    Conductances(g_Na, g_K, g_L)
end

struct MhnParameters
    m_inf
    h_inf
    n_inf
    m_tau
    h_tau
    n_tau
end

function MhnParameters(V)
    a_m = 0.1 * (V + 40) / (1 - exp((-V - 40) / 10))
    b_m = 4 * exp((-V - 65) / 18)
    a_h = 0.07 * exp((-V - 65) / 20)
    b_h = 1 / (1 + exp((-V - 35) / 10))
    a_n = 0.01 * (V + 55) / (1 - exp((-V - 55) / 10))
    b_n = 0.125 * exp((-V - 65) / 80)

    m_tau = 1 / (a_m + b_m)
    h_tau = 1 / (a_h + b_h)
    n_tau = 1 / (a_n + b_n)

    m_inf = m_tau * a_m
    h_inf = h_tau * a_h
    n_inf = n_tau * a_n

    MhnParameters(m_inf, h_inf, n_inf, m_tau, h_tau, n_tau)
end

struct InversionV
    V_Na
    V_K
    V_L
end

function InversionV()
    InversionV(V_Na, V_K, V_L)
end

struct ODEParams
    conductances::Conductances
    inversion_v::InversionV
    I
    C
end

function f(dyk, yk, params::ODEParams, tk)
    V, m, h, n = yk
    (; V_Na, V_K, V_L) = params.inversion_v
    (; g_Na, g_K, g_L) = params.conductances
    (; m_inf, h_inf, n_inf, m_tau, h_tau, n_tau) = MhnParameters(V)

    dyk[1] = (params.I - g_Na * m^3 * h * (V - V_Na) - g_K * n^4 * (V - V_K) - g_L * (V - V_L)) / params.C
    dyk[2] = (m_inf - m) / m_tau
    dyk[3] = (h_inf - h) / h_tau
    dyk[4] = (n_inf - n) / n_tau
end

function fill_defaults(cfg::Dict{String,Any})
    return Dict(
        "V0" => get(cfg, "V0", default_V0),
        "I" => get(cfg, "I", default_I),
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

function solve_system(V0, I, g, v, C, tmax)
    (; m_inf, h_inf, n_inf) = MhnParameters(V0)
    Y0 = [V0, m_inf, h_inf, n_inf]
    params = ODEParams(g, v, I, C)
    tspan = (0, tmax)

    prob = ODEProblem(f, Y0, tspan, params)
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

    V0 = cfg["V0"]
    I = cfg["I"]
    g = Conductances()
    v = InversionV()

    sol = solve_system(V0, I, g, v, C, tmax)
    PlotSolution(sol, fig_path * string(i) * ".pdf", (; V0, I))
    plot(sol.t)
end

gui()

println("Press the enter key to quit:")
readline()