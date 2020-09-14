
"""
    leaftemp(p, v)

Calculate the leaf temperature given variables `v`.

Returns a value in u"K".
"""
@inline leaftemp(p, v) = v.tair + (v.rnet - v.et * v.lhv) / (4CPAIR * AIRMA * v.gh)

"""
    latent_heat_water_vapour(tair)

Caculates the late hear of water vapour from the air temperature.

Returns a value in u"K".
"""
@inline latent_heat_water_vapour(tair) = begin
    T = ustrip(°C, tair)
    (H2OLV0 - 2.365e3J*kg^-1 * T) * H2OMW
end


"""
    vapour_pressure_deficit(tair, rh)

Calculate vapour pressure deficit given air temperature
`tair` and relative humidity `rh`

Returns a value in u"kPa"
"""
vapour_pressure_deficit(tair, rh) = begin
    es = saturated_vapour_pressure(tair)
    ea = rh * es
    ea - es
end


"""
    arrhenius(kref, Ea, T, Tref)

The Arrhenius function.

-`kref`: is the value at Tref deg
-`Ea`: the activation energy (j mol - 1) and 
-`T`: the temp in °C or K
-`Tref`: the reference temp in °C or K
Standard form and temperature difference form.
"""
arrhenius(kref, Ea, T, Tref) = arrhenius(kref, Ea, T |> K, Tref |> K)
arrhenius(kref, Ea, T::typeof(1.0K), Tref::typeof(1.0K)) = kref * exp(Ea * (T - Tref) / (R * T * Tref))
arrhenius(A, Ea, T) = arrhenius(A, Ea, T |> K)
arrhenius(A, Ea, T::typeof(1.0K)) = A * exp(Ea / (R * T))
