# Energy balance interface

abstract type AbstractEnergyBalance end

"""
    model_init!(v, f)
Runs any model initialisation that needs to happen at the start of energy balance
"""
function enbal_init! end

"""
    model_update!(f, v, tleaf)
Runs any model specific variable updates that need to happen at the end of
the leaf temperature convergene loop
"""
function enbal_update! end

"""
    enbal!(p, v)
This subroutine calculates leaf photosynthesis and transpiration.

These may be calculated by:
(1) assuming leaf temperature = air temperature, cs = ca and ds = da
(2) using iterative scheme of Leuning et al (1995) (PCE 18:1183-1200) to calculate leaf temp, CsCa.

Setting itermax = 0 gives (1); itermax > 0 (suggest 100) gives (2).
"""
function enbal! end
