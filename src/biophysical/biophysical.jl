
@inline leaftemp(p, v) = v.tair + (v.rnet - v.et * v.lhv) / (CPAIR * AIRMA * v.gh)

"""
    latent_heat_water_vapour(tair)

Caculates the late hear of water vapour from the air temperature.
"""
@inline latent_heat_water_vapour(tair) = (H2OLV0 - 2.365e3J*kg^-1*K^-1 * tair) * H2OMW


"""
Calculate saturated water vapour pressure (Pa) at temperature Ta (Celsius)
from Jones 1992 p 110 (note error in a - wrong units)
TODO: what are thes magic numbers? Name them.
"""
@inline saturated_vapour_pressure(tair) = 613.75Pa * exp(17.502 * tair / (K(240.97°C) + tair))

vapour_pressure_deficit(tair, rh) = begin
    es = saturated_vapour_pressure(tair)
    ea = rh * es
    ea - es
end


"""
    arrhenius(kt, ea, t, tref)
The Arrhenius function.
kT is the value at tref deg # 
Ea the activation energy (j mol - 1) and 
T the temp (deg #).
Standard form and temperature difference form.
"""
arrhenius(A, Ea, T) = arrhenius(A, Ea, T |> K)
arrhenius(A, Ea, T::typeof(1.0K)) = A * exp(Ea / (R * T))
arrhenius(kref, Ea, T, Tref) = arrhenius(kref, Ea, T |> K, Tref |> K)
arrhenius(kref, Ea, T::typeof(1.0K), Tref::typeof(1.0K)) = kref * exp(Ea * (T - Tref) / (R * T * Tref))
