
"""
Respiration models, run in a [`respiration`](@ref) method.
"""
abstract type AbstractRespiration end

""" 
    respiration(f::Union{Nothing,Respiration}, tleaf)

Respiration for formulation `f` at temperature `tleaf` in `u"K"`.
Returns respiration in `u"μmol*m^-2*s^-1"`.
"""
function respiration end

@mix @columns struct MixinResp{pK,K,F,μMoM2S}
#   Field       | Default | Unit           | Bonds          | Description
    q10f::pK    | 0.67    | K^-1           | (0.0, 1.0)     | "Logarithm of the Q10"
    dayresp::F  | 1.0     | _              | (0.0, 1.0)     | "Respiration in the light as fraction of that in the dark"
    rd0::μMoM2S | 0.001   | μmol*m^-2*s^-1 | (0.0, 0.1)     | "Dark respiration at the reference temperature"
    tbelow::K   | 173.15  | K              | (250.0, 300.0) | "Temperature below which no respiration occurs"
    tref::K     | 298.15  | K              | (250.0, 350.0) | "Reference temperature at which rd0 was measured"
end

respiration(f::Nothing, tleaf) = 0.0μmol*m^-2*s^-1

"""
    Respiration(q10f, dayresp, rd0, tbelow, tref)

Standard respiration model.
Calculates respiration from temperature using a Q10 (exponential) formulation 

TODO: specify origin of formulation
"""
@MixinResp struct Respiration{} <: AbstractRespiration end

respiration(f::Respiration, tleaf) = begin
    tleaf < f.tbelow && return zero(f.rd0)
    f.rd0 * exp(f.q10f * (tleaf - f.tref)) * f.dayresp
end

"""
    AcclimatizedRespiration(k10f, tmove, q10f, dayresp, rd0, tbelow, tref)

Respiration with acclimatization parameters `k10f` and `tmove`. 

TODO test this. The formulation appears to need a variable `tmove`, 
not a parameter.
"""
@MixinResp struct AcclimatizedRespiration{pK,K} <: AbstractRespiration
#   Field    | Default | Unit | Bonds       | Description
    k10f::pK | 10.0    | K^-1 | (0.0, 10.0) | _
    tmove::K | 10.0    | K    | (0.0, 10.0) | _
end

respiration(f::AcclimatizedRespiration, tleaf) = begin
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
#
