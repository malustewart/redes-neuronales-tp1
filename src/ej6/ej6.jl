using Plots
using JSON
using ArgParse
using DifferentialEquations

const default_V_init = 0
const default_A_init = 0.5
const default_A_0 = 0.1
const default_V_thres = 1
const default_tau = 1
const default_tauA = 10
const default_V_rest = 0
const default_R = 1
const default_I = 1

const default_tmax = 300

struct Params
    A_0
    V_thres
    tau
    tauA
    V_rest
    R
    I
    reset_flag
end

function fill_defaults(cfg::Dict{String,Any})
    return Dict(
        "V_init"        => get(cfg, "V_init", default_V_init),
        "A_init"        => get(cfg, "A_init", default_A_init),
        "A_0"           => get(cfg, "A_0", default_A_0),
        "V_thres"       => get(cfg, "V_thres", default_V_thres),
        "tau"           => get(cfg, "tau", default_tau),
        "tauA"          => get(cfg, "tauA", default_tauA),
        "V_rest"        => get(cfg, "V_rest", default_V_rest),
        "R"             => get(cfg, "R", default_R),
        "I"             => get(cfg, "I", default_I),
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

condition(u, t, integrator) = u[1] - 1

function affect!(integrator) 
    println(integrator.t[1])
    integrator.u[1] = 0
    integrator.u[2] += integrator.p.A_0 / integrator.p.tauA
end

function adap_leaky_iaf(dyk, yk, params::Params, tk)
    V = yk[1]
    A = yk[2]
    dyk[1] = (-(V - params.V_rest) + params.R * (params.I + A)) / params.tau
    dyk[2] = -A / params.tauA
end

args = parse_commandline()
config_path = args["config"]
fig_path = args["figpath"]

config = read_config(config_path)
tmax = get(config, "tmax", default_tmax)

tspan = (0, tmax)

for (i, plot_cfg) in enumerate(config["plots"])
    plt = plot(layout=(1,2))

    cfg = fill_defaults(plot_cfg)
    V_init = cfg["V_init"]
    A_init = cfg["A_init"]
    A_0 = cfg["A_0"]
    V_thres = cfg["V_thres"]
    tau = cfg["tau"]
    tauA = cfg["tauA"]
    V_rest = cfg["V_rest"]
    R = cfg["R"]
    I = cfg["I"]
    reset_flag = false
    
    params = Params(A_0, V_thres, tau, tauA, V_rest, R, I, reset_flag)

    Y0 = [V_init, A_init]

    cb = ContinuousCallback(condition, affect!)

    prob = ODEProblem(adap_leaky_iaf, Y0, tspan, params)
    sol = solve(prob, callback = cb)

    plot!(sol[1, :], sol[2,:], subplot=1, xlabel="V", ylabel="A")
    plot!(sol.t, sol[1, :], ylabel="V", subplot=2, xlabel="time")
    plot!(sol.t, sol[2, :], ylabel="A", subplot=2, xlabel="time")
    savefig(plt, fig_path * string(i) * ".pdf")
end
gui()

println("Press the enter key to quit:")
readline()
