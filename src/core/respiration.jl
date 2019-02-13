abstract type AbstractRespiration end
@mix @columns struct Resp{pK,K,F,μMoM2S}
    q10f::pK        | 0.67   | K^-1           | _                | (0.0, 1.0)     | "logarithm of the Q10 (Equation for respiration)"
    dayresp::F      | 1.0    | _              | Beta(5, 1)       | (0.0, 1.0)     | "Respiration in the light as fraction of that in the dark."
    rd0::μMoM2S     | 0.9    | μmol*m^-2*s^-1 | Gamma(10, 0.9/10)| (0.0, 1.0)     | _
    tbelow::K       | 173.15 | K              | _                | (250.0, 300.0) | "No respiration occurs below this temperature."
    tref::K         | 298.15 | K              | _                | (250.0, 350.0) | "Reference temperature (T at which RD0 was measured)"
end

@Resp struct Respiration{} <: AbstractRespiration end
@Resp struct AcclimatizedRespiration{pK,K} <: AbstractRespiration
    k10f::pK        | 0.0    | K^-1           | _             | (0.0, 10.0) | _
    tmove::K        | 1.0    | K              | Gamma(2, 1/2) | (0.0, 10.0) | _
end

""" Calculates respiration from temperature using a Q10 (exponential) formulation """
respiration!(f::Nothing, v, rd) = zero(rd) 

function respiration(f::Respiration, v)
    v.tleaf < f.tbelow && return zero(v.rd)
    # Make sure light suppression of dark respiration only occurs when it is light.
    # See Atkin et al. 1998 (or 2001?). From yplantmc
    resplightfrac = v.par < 100oneunit(v.par) ? oneunit(f.dayresp) : f.dayresp

    f.rd0 * exp(f.q10f * (v.tleaf - f.tref)/10) * resplightfrac 
end

function respiration(f::AcclimatizedRespiration, v)
    v.tleaf < f.tbelow && return zero(v.rd)
    resplightfrac = v.par < 100oneunit(v.par) ? oneunit(f.dayresp) : f.dayresp

    rd0acc = f.rd0 * exp(f.k10f * (f.tmove - f.tref))
    rd0acc * exp(f.q10f * (v.tleaf - f.tref)) * resplightfrac 
end 
