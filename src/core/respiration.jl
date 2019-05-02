abstract type AbstractRespiration end
@mix @columns struct Resp{pK,K,F,μMoM2S}
    q10f::pK        | 0.67   | K^-1           | _                | (0.0, 1.0)     | "Logarithm of the Q10"
    dayresp::F      | 1.0    | _              | Beta(5, 1)       | (0.0, 1.0)     | "Respiration in the light as fraction of that in the dark"
    rd0::μMoM2S     | 0.001  | μmol*m^-2*s^-1 | Gamma(10, 0.9/10)| (0.0, 0.1)     | "Dark respiration at the reference temperature"
    tbelow::K       | 173.15 | K              | _                | (250.0, 300.0) | "Temperature below which no respiration occurs"
    tref::K         | 298.15 | K              | _                | (250.0, 350.0) | "Reference temperature at which rd0 was measured"
end

@Resp struct Respiration{} <: AbstractRespiration end
@Resp struct AcclimatizedRespiration{pK,K} <: AbstractRespiration
    k10f::pK        | 10.0   | K^-1           | _             | (0.0, 10.0) | _
    tmove::K        | 10.0   | K              | Gamma(2, 1/2) | (0.0, 10.0) | _
end

struct FractionalRespiration{R} <: AbstractRespiration 
    formulation::R
end

""" Calculates respiration from temperature using a Q10 (exponential) formulation """
respiration!(f::Nothing, v, rd) = zero(rd) 

# What is this 10 number???
function respiration(f::Respiration, tleaf)
    tleaf < f.tbelow && return zero(f.rd0)

    f.rd0 * exp(f.q10f * (tleaf - f.tref)) * f.dayresp
end

function respiration(f::AcclimatizedRespiration, tleaf)
    tleaf < f.tbelow && return zero(v.rd)

    rd0acc = f.rd0 * exp(f.k10f * (f.tmove - f.tref))
    rd0acc * exp(f.q10f * (tleaf - f.tref)) * f.dayresp 
end 

# dayresp(f) = f.dayresp

# Make sure light suppression of dark respiration only occurs when it is light.
# See Atkin et al. 1998 (or 2001?). From yplantmc
# TODO: The cutoff should be a parameter
# lightfrac = v.par < 100oneunit(v.par) ? oneunit(dayresp(f.formulation)) : dayresp(f.formulation)
# respiration(f.formulation, v) * lightfrac 
