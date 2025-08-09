using Plots
using ArgParse

default_alpha = 0.66666
default_beta = 1.33333
default_gamma = 1
default_delta = 1
default_r0 = 3.0
default_f0 = 1.5
default_steps = 1000
default_last_t = 25
default_r0 = 1.0
default_f0 = 1.0

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--alpha", "-a"
            help = "alpha"
            arg_type = Float64
            default = default_alpha
        "--beta", "-b"
            help = "beta"
            arg_type = Float64
            default = default_beta
        "--gamma", "-c"
            help = "gamma"
            arg_type = Float64
            default = default_gamma
        "--delta", "-d"
            help = "delta"
            arg_type = Float64
            default = default_delta
        "--steps", "-s"
            help = "Number of steps to run"
            arg_type = Int
            default = default_steps
        "-t"
        help = "End time"
            arg_type = Float64
            default = default_last_t
        "-r"
            help = "Rabbit initial population"
            arg_type = Float64
            default = default_r0
        "-f"
            help = "Fox initial population"
            arg_type = Float64
            default = default_f0
    end

    return parse_args(s)
end

parsed_args = parse_commandline()

alpha = parsed_args["alpha"]
beta = parsed_args["beta"]
gamma = parsed_args["gamma"]
delta = parsed_args["delta"]
r0 =parsed_args["r"]
f0 = parsed_args["f"]
steps = parsed_args["steps"]
last_t = parsed_args["t"]
h = last_t / steps

function lotka_volterra(tk, yk)
  [(alpha - beta*yk[2])*yk[1], (gamma * yk[1] - delta)*yk[2]]
end

function rk4(yk, tk, h, f)
  q1 = f(tk,       yk)
  q2 = f(tk + h/2, yk + h/2*q1)
  q3 = f(tk + h/2, yk + h/2*q2)
  q4 = f(tk + h,   yk + h*q3)

  yk + h * (q1 + 2*q2 +  2*q3 + q4) / 6
end

t = range(0, convert(Int, last_t), length=steps)
Y = zeros(Float64, steps + 1, 2)

Y[1,:] = [r0,f0]

for (k, tk) in enumerate(t)
  Y[k+1, :] =  rk4(Y[k, :], tk, h, lotka_volterra)
end


p = plot!(Y[1:steps,1], Y[1:steps,2], layout=(1,2), subplot=1, xlabel="rabbits", ylabel="foxes")
p = plot!(t, [Y[1:steps,1] Y[1:steps, 2]], label=["rabbit" "foxes"],subplot=2, xlabel="time")
gui()
println("Press the enter key to quit:")
readline()