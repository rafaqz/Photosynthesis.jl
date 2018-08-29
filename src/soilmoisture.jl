
"""
Calculate the effect of soil water content on stomatal conductance

These methods dispatch on AbstractSoilMethod or AbstractSoilData types.
"""
function soil_soil_water_conductance end

soil_water_conductance!(v, f::ConstantSoilMethod, p) = v.fsoil = 1.0

"Convert volumetric data to deficit"
soil_water_conductance!(v, f::DeficitSoilMethod{<:AbstractContentSoilData}, p) = begin
    s = f.soildata
    soilmoist = (s.swmax - v.soilmoist) / (s.swmax - s.swmin)
    v.fsoil = soil_water_conductance!(s, soilmoist, v, p)
end

soil_water_conductance!(v, f::DeficitSoilMethod{<:DeficitSoilData}, p) = 
    v.fsoil = soil_water_conductance(f.soildata, v.soilmoist, v, p)

soil_water_conductance!(v, f::PotentialSoilMethod{<:PotentialSoilData}, p) = 
    v.fsoil = soil_water_conductance(f.soildata, v.soilmoist, v, p)

soil_water_conductance!(v, f::VolumetricSoilMethod{<:SimulatedSoilData}, p) =
    v.soilmoist = f.soilroot / f.soildepth # TODO fix this, where is fsoil set??

soil_water_conductance!(v, f::VolumetricSoilMethod{<:AbstractContentSoilData}, p) = begin
    # TODO: check in constructor f.wc1 > f.wc2 && error("error: wc1 needs to be smaller than wc2.")
    fsoil = -f.wc1 / (f.wc2 - f.wc1) + v.soilmoist / (f.wc2 - f.wc1)
    v.fsoil = max(min(fsoil, 1.0), 0.0)
end

soil_water_conductance!(v, f::TuzetSoilMethod, p) = begin
    # What is this from? v.psif = (1 + e^(v.sf * v.psiv)) / (1 + e^(v.sf * (v.psiv - v.psil)))
    vjmod = vjmax_water(f.soildata, v, p)
    v.vcmax *= vjmod
    v.jmax *= vjmod
    v.fsoil = 1.0
end

soil_water_conductance!(v, f::EmaxSoilMethod, p) = begin
    vjmod = vjmax_water(f.soildata, v, p)
    v.vcmax *= vjmod
    v.jmax *= vjmod
    v.fsoil = 1.0
end

"""
    soil_water_conductance(soil::PotentialSoilData, soilmoist, v, p)
Negative exponential function of soil water potential. 
"""
soil_water_conductance(soil::PotentialSoilData, soilmoist, v, p) = begin
    # Exponential relationship with potential: parameter = SWPEXP
    # TODO: default value for effect? 
    if soil.swpexp > zero(soil.swpexp)
        effect = e^((soil.swpexp/oneunit(soil.swpexp) * soilmoist/oneunit(soilmoist)))
    end
    return max(effect, zero(effect))
end

"""
    soil_water_conductance(soil::AbstractDeficitSoilData, soilmoist, v, p)
GranierLoustau 1994 Fs = 1 - a exp(b SMD) where SMD is soil
moisture deficit, dimless """
soil_water_conductance(soil::AbstractDeficitSoilData, soilmoist, v, p) = begin
    effect = 1.0

    # Exponential relationship with deficit: params SMD1, SMD2
    if soil.smd1 > zero(soil.smd1)
        effect = 1.0 - soil.smd1 / oneunit(soil.smd1) * e^(soil.smd2 * soilmoist / oneunit(soilmoist)^2)
    # Linear decline with increasing deficit: pUT SMD1 < 0, param SMD2
    elseif soil.smd2 > zero(soil.smd2)
        if oneunit(soilmoist) - soilmoist < soil.smd2
            effect = (oneunit(soilmoist) - soilmoist) / soil.smd2
        end
    end
    return max(effect, zero(effect))
end


"""
Dependance of assimilation on soil water potential (SWP).
Modifier function (return value between 0.0 and 1.0) 
"""
function vjmax_water end

"""
    vjmax_water(f::NoPotentialDependence, v, p)
Returns 1.0
"""
vjmax_water(f::NoPotentialDependence, v, p) = 1.0

"""
    vjmax_water(f::LinearPotentialDependence, v, p)
Simple linear dependance.
vpara is swp where vcmax is zero, vparb is swp where vcmax is one.
"""
vjmax_water(f::LinearPotentialDependence, v, p) =
    if v.weightedswp < f.vpara 
        0.0
    elseif v.weightedswp > f.vparb 
        1.0
    else 
        (v.weightedswp - f.vpara) / (f.vparb - f.vpara)
    end

"""
    vjmaxw(f::ZhouPotentialDependence, v, p)
Zhou, et al. Agricultural and Forest Meteorology. 2013.0
"""
vjmax_water(f::ZhouPotentialDependence, v, p) = (1 + e^(f.vpara * f.vparb)) /
    (1 + e^(f.vpara * (f.vparb - v.weightedswp)))

