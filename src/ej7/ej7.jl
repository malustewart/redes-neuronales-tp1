using Plots
using JSON
using ArgParse

default_a = 0.7
default_b = 0.8
default_gamma = 0.5
default_I = 0.5
default_step = 0.001
default_tmax = 100
default_tau = 1
default_tauw = 12.5
default_V0 = 2
default_w0 = 1

struct Params
    a::Float64
    b::Float64
    gamma::Float64
    I::Float64
    tau::Float64
    tauw::Float64
    V0::Float64
    w0::Float64
end

function fill_defaults(cfg::Dict{String,Any})
    return Dict(
        "a"         => get(cfg, "a", default_a),
        "b"         => get(cfg, "b", default_b),
        "gamma"     => get(cfg, "gamma", default_gamma),
        "I"         => get(cfg, "I", default_I),
        "V0"        => get(cfg, "V0", default_V0),
        "w0"        => get(cfg, "w0", default_w0),
        "tau"       => get(cfg, "tau", default_tau),
        "tauw"       => get(cfg, "tauw", default_tauw),
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
            default = "imgs/lastfig.pdf"
    end
    return parse_args(s)
end

f(V, a) = V*(a-V)*(V-1)

function fhn(tk, yk, params)
    [(f(yk[1], params.a)) + params.I - yk[2], params.gamma * yk[2] + params.b*yk[1]]./[params.tau, params.tauw]
end

function rk4(yk, tk, h, f, params)
    q1 = f(tk,       yk,          params)
    q2 = f(tk + h/2, yk + h/2*q1, params)
    q3 = f(tk + h/2, yk + h/2*q2, params)
    q4 = f(tk + h,   yk + h*q3,   params)

    yk + h * (q1 + 2*q2 +  2*q3 + q4) / 6
end

args = parse_commandline()
config_path = args["config"]
fig_path = args["figpath"]

config = read_config(config_path)
tmax = get(config, "tmax", default_tmax)
s = get(config, "step", default_step)
len = convert(Int, round(tmax/s))

t = range(0, tmax, length=len)

plt = plot(layout=(1,2))

for (i, plot_cfg) in enumerate(config["plots"])
    cfg = fill_defaults(plot_cfg)
    V0 = cfg["V0"]
    w0 = cfg["w0"]
    a = cfg["a"]
    I = cfg["I"]
    gamma = cfg["gamma"]
    tau = cfg["tau"]
    tauw = cfg["tauw"]
    b = cfg["b"]
    
    params = Params(a, b, gamma, I, tau, tauw, V0, w0)

    Y = zeros(Float64, len + 1, 2)
    Y[1,:] = [V0,w0]
    
    for (k, tk) in enumerate(t)
      Y[k+1, :] =  rk4(Y[k, :], tk, s, fhn, params)
    end
    
    plot!(Y[1:len,1], Y[1:len,2], layout=(1,2), subplot=1, xlabel="V", ylabel="w")
    plot!(t, [Y[1:len,1] Y[1:len, 2]], label=["V" "w"],subplot=2, xlabel="time")
end
gui()

println("Press the enter key to quit:")
readline()
savefig(plt, fig_path)