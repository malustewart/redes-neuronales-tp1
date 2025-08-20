module HH

export default_k, default_I, default_V0, default_tmax, default_I_min, default_I_max, C, Conductances, InversionV, MhnParameters, HHParams, HH_models

default_k = 1.0 # h+n = k in aprox 2
default_V0 = -80
default_I = 1
default_I_min = -100
default_I_max = 300
default_tmax = 300

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

struct HHParams
    conductances::Conductances
    inversion_v::InversionV
    I
    C
    k
end

function HH_standard(dyk, yk, params::HHParams, tk)
    V, m, h, n = yk
    (; V_Na, V_K, V_L) = params.inversion_v
    (; g_Na, g_K, g_L) = params.conductances
    (; m_inf, h_inf, n_inf, m_tau, h_tau, n_tau) = MhnParameters(V)

    dyk[1] = (params.I - g_Na * m^3 * h * (V - V_Na) - g_K * n^4 * (V - V_K) - g_L * (V - V_L)) / params.C
    dyk[2] = (m_inf - m) / m_tau
    dyk[3] = (h_inf - h) / h_tau
    dyk[4] = (n_inf - n) / n_tau
end

function HH_aprox_1(dyk, yk, params::HHParams, tk)
    V, _, h, n = yk
    (; V_Na, V_K, V_L) = params.inversion_v
    (; g_Na, g_K, g_L) = params.conductances
    (; m_inf, h_inf, n_inf, m_tau, h_tau, n_tau) = MhnParameters(V)

    m = m_inf

    dyk[1] = (params.I - g_Na * m^3 * h * (V - V_Na) - g_K * n^4 * (V - V_K) - g_L * (V - V_L)) / params.C
    dyk[2] = (m_inf - m) / m_tau    # =0
    dyk[3] = (h_inf - h) / h_tau
    dyk[4] = (n_inf - n) / n_tau
end

function HH_aprox_2(dyk, yk, params::HHParams, tk)
    V, _, _, n = yk
    (; V_Na, V_K, V_L) = params.inversion_v
    (; g_Na, g_K, g_L) = params.conductances
    (; m_inf, h_inf, n_inf, m_tau, h_tau, n_tau) = MhnParameters(V)

    m = m_inf
    h = params.k - m

    dyk[1] = (params.I - g_Na * m^3 * h * (V - V_Na) - g_K * n^4 * (V - V_K) - g_L * (V - V_L)) / params.C
    dyk[2] = (m_inf - m) / m_tau    # =0
    dyk[3] = (h_inf - h) / h_tau
    dyk[4] = (n_inf - n) / n_tau
end

HH_models = Dict(
    :HH_standard => HH_standard,
    :HH_aprox_1 => HH_aprox_1,
    :HH_aprox_2 => HH_aprox_2,
)

end