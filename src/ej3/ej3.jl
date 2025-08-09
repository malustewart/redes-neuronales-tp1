using Plots
using JSON
using ArgParse

default_temp = 37.5 # C
default_conc_in = 1 # mol/l
default_conc_out = 2    # mol/l
default_permeability = 2e-8 # cm/s
default_valence = 1
default_vmin = -200
default_vmax = 200
default_step = 1
default_name = "Unnamed"

R = 8.314
F = 96485

function fill_defaults(cfg::Dict{String,Any})
    return Dict(
        "name"         => get(cfg, "name", default_name),
        "temp"         => get(cfg, "temp", default_temp) + 273, # convert from C to K
        "conc_in"      => get(cfg, "conc_in", default_conc_in) / 1000, # convert from mol/L to mol/cm3
        "conc_out"     => get(cfg, "conc_out", default_conc_out) / 1000, # convert from mol/L to mol/cm3
        "permeability" => get(cfg, "permeability", default_permeability),   # cm/s
        "valence"      => get(cfg, "valence", default_valence),
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

function GHK_i(Vk, temp, conc_in, conc_out, permeability, valence)
    zF_RT = valence * F / R / temp
    [permeability .* valence .* F .* zF_RT .* Vk .* (conc_in .- conc_out * exp.(-zF_RT.*Vk))./ (1 .- exp.(-zF_RT.*Vk))]
end

args = parse_commandline()
config_path = args["config"]
fig_path = args["figpath"]

config = read_config(config_path)
vmin = get(config, "Vmin", default_vmin) / 1000
vmax = get(config, "Vmax", default_vmax) / 1000
s = get(config, "step", default_step) / 1000

len = convert(Int,round((vmax-vmin)/s))
v = range(vmin, step=s, length=len)

plt = plot(xlabel="V (mV)", ylabel="I (A/cm2)")

for (i, plot_cfg) in enumerate(config["plots"])
    cfg = fill_defaults(plot_cfg)
    name          = cfg["name"]
    temp          = cfg["temp"]
    conc_in       = cfg["conc_in"]
    conc_out      = cfg["conc_out"]
    permeability  = cfg["permeability"]
    valence       = cfg["valence"]

    i = GHK_i(v, temp, conc_in, conc_out, permeability, valence)

    plot!(1000 .*v, i, label=name)
end

gui()

println("Press the enter key to quit:")
readline()
savefig(plt, fig_path)