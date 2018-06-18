
abstract type AbstractSoilData end
abstract type AbstractDeficitSoilData <: AbstractSoilData end
abstract type AbstractContentSoilData <: AbstractDeficitSoilData end

" @Deficit mixin macro adds water deficit fields to a struct"
@mix @with_kw struct Deficit{T} 
    smd1::T = 1.0
    smd2::T = 1.0
end

"@Content mixin macro adds water content fields to a struct"
@Deficit @mix struct Content{W}
    swmax::W = 1.0
    swmin::W = 0.0
end

""" 
Soil water is measured as deficit 
"""

@Deficit mutable struct DeficitSoilData{} <: AbstractDeficitSoilData end
""" Soil water is measured as volumetric content 
"""
@Content mutable struct ContentSoilData{} <: AbstractContentSoilData end
""" Simulated soil volumetric content 
"""
@Content mutable struct SimulatedSoilData{} <: AbstractContentSoilData end
""" Soil water is measured as water potential 
"""
@with_kw mutable struct PotentialSoilData{T} <: AbstractSoilData
    swpexp::T = 1.0
end
mutable struct NoSoilData <: AbstractSoilData end


abstract type AbstractPotentialDependence end

"@Potdep mixin macro adds potential dependence fields to a struct"
@mix @with_kw struct Potdep
    vpara::typeof(1.0u"kPa") = 1.0u"kPa"
    vparb::typeof(1.0u"kPa") = 2.0u"kPa"
end

@Potdep mutable struct LinearPotentialDependence{} <: AbstractPotentialDependence end
@Potdep mutable struct ZhouPotentialDependence{} <: AbstractPotentialDependence end
type NoPotentialDependence{} <: AbstractPotentialDependence end


abstract type AbstractSoilMethod{T} end

@with_kw mutable struct VolumetricSoilMethod{WC,R,D,T<:AbstractContentSoilData
                                            } <: AbstractSoilMethod{T}
    soildata::T  = ContentSoilData()
    wc1::WC      = 0.5
    wc2::WC      = 0.5
    soilroot::R  = 0.5
    soildepth::D = 0.5
end

@with_kw mutable struct ConstantSoilMethod{T<:NoSoilData} <: AbstractSoilMethod{T}
    soildata::T = NoSoilData()
end
@with_kw mutable struct DeficitSoilMethod{T<:AbstractDeficitSoilData} <: AbstractSoilMethod{T}
    soildata::T = ContentSoilData()
end
@with_kw mutable struct PotentialSoilMethod{T<:PotentialSoilData} <: AbstractSoilMethod{T}
    soildata::T = PotentialSoilData()
end

@mix @with_kw struct MaespaSoil{T<:AbstractPotentialDependence}
    soildata::T = LinearPotentialDependence()
end

@MaespaSoil mutable struct EmaxSoilMethod{} <: AbstractSoilMethod{T} end
@MaespaSoil mutable struct TuzetSoilMethod{} <: AbstractSoilMethod{T} end


"Physiological constants for CO2 snd rubisco comprensation points"
abstract type AbstractCompensation end

@with_kw struct BadgerCollatzCompensation <: AbstractCompensation 
    Kc25::typeof(1.0u"μmol*mol^-1") = 404.0u"μmol*mol^-1" # MM coefft of Rubisco for CO2 (umol mol - 1)
    Ko25::typeof(1.0u"μmol*mol^-1") = 248.0u"mmol*mol^-1" # MM coefft of Rubisco for O2 (umol mol - 1)
    ΔHa_Kc::typeof(1.0u"kJ*mol^-1") = 59.4u"kJ*mol^-1"  # Temp. response of Kc (J mol-1)
    ΔHa_Ko::typeof(1.0u"kJ*mol^-1") = 36.0u"kJ*mol^-1"  # Temp. response of Ko (J mol-1)
    tref::typeof(1.0u"°C")          = 25.0u"°C"
end

"""
Bernacchi et al 2001 PCE 24: 253 - 260
Extra deactivation terms may be required above 40°C.
"""
@with_kw struct BernacchiCompensation <: AbstractCompensation 
    Kc25::typeof(1.0u"μmol*mol^-1") = 404.9u"μmol*mol^-1" # MM coefft of Rubisco for CO2
    Ko25::typeof(1.0u"μmol*mol^-1") = 278.4u"mmol*mol^-1" # MM coefft of Rubisco for O2
    Γ☆25::typeof(1.0u"μmol*mol^-1") = 42.75u"μmol*mol^-1"
    ΔHa_Kc::typeof(1.0u"kJ*mol^-1") = 79.43u"kJ*mol^-1" # Temp. response of Kc
    ΔHa_Ko::typeof(1.0u"kJ*mol^-1") = 36.38u"kJ*mol^-1" # Temp. response of Ko
    ΔHa_Γ☆::typeof(1.0u"kJ*mol^-1") = 37.83u"kJ*mol^-1"
    tref::typeof(1.0u"°C")          = 25.0u"°C"
end


""" Electron flux parameters """
abstract type AbstractJmax end

@with_kw mutable struct Jmax <: AbstractJmax 
    jmax25::typeof(1.0u"μmol*m^-2*s^-1") = 184.0u"μmol*m^-2*s^-1"
    delsj::typeof(1.0u"J*mol^-1*K^-1")   = 640.02u"J*mol^-1*K^-1" # DELTAS in Medlyn et al. (2002)
    eavj::typeof(1.0u"J*mol^-1")         = 37259.0u"J*mol^-1" # Ha in Medlyn et al. (2002)
    edvj::typeof(1.0u"J*mol^-1")         = 200000.0u"J*mol^-1" # Hd in Medlyn et al. (2002)
end


abstract type AbstractVcmax end

"@Vcmax mixin macro adds vcmax fields to a struct"
@mix @with_kw struct Vcmax
    vcmax25::typeof(1.0u"μmol*m^-2*s^-1") = 110.0u"μmol*m^-2*s^-1" 
    eavc::typeof(1.0u"J*mol^-1")          = 47590.0u"J*mol^-1" # Ha in Medlyn et al. (2002)
end

@Vcmax mutable struct NoOptimumVcmax <: AbstractVcmax end
@Vcmax mutable struct OptimumVcmax <: AbstractVcmax
    edvc::typeof(1.0u"J*mol^-1")       = 1.0u"J*mol^-1" # Hd in Medlyn et al. (2002)
    delsc::typeof(1.0u"J*mol^-1*K^-1") = 629.26u"J*mol^-1*K^-1" # DELTAS in Medlyn et al. (2002)
end


abstract type AbstractVcJmax end

@with_kw mutable struct VcJmax{J,V} <: AbstractVcJmax where {J<:AbstractVcmax,V<:AbstractJmax}
    jmaxformulation::J      = Jmax()
    vcmaxformulation::V     = OptimumVcmax()
end
@with_kw mutable struct DukeVcJmax{J,V} <: AbstractVcJmax where {J<:AbstractVcmax,V<:AbstractJmax}
    jmaxformulation::J      = Jmax()
    vcmaxformulation::V     = OptimumVcmax()
    tvjup::typeof(1.0u"°C") = 10.0u"°C"
    tvjdn::typeof(1.0u"°C") = 0.0u"°C"
end


abstract type AbstractRubiscoRegen end

@with_kw mutable struct RubiscoRegen <: AbstractRubiscoRegen
    theta::Float64  = 0.4   # Shape parameter of the non-rectangular hyperbola.
    ajq::Float64    = 0.324 # Quantum yield of electron transport.
end


abstract type AbstractRespiration end

@mix @with_kw struct Resp
    q10f::typeof(1.0u"K^-1")          = 0.67u"K^-1" # logarithm of the Q10 (Equation for respiration)
    rtemp::typeof(1.0u"°C")           = 25.0u"°C"   # Reference temperature (T at which RD0 was measured)
    tbelow::typeof(-100.0u"°C")       = -100.0u"°C" # No respiration occurs below this temperature (degC).
	dayresp::typeof(1.0)              = 1.0         # Respiration in the light as fraction of that in the dark.
    rd0::typeof(1.0u"μmol*m^-2*s^-1") = 0.9u"μmol*m^-2*s^-1"
end

@Resp mutable struct Respiration <: AbstractRespiration end
@Resp mutable struct AcclimatizedRespiration <: AbstractRespiration 
    k10f::typeof(1.0u"K^-1") = 0.0u"K^-1"
    tmove::typeof(1.0u"K")   = 1.0u"K"
end


abstract type AbstractStomatalConductance end

" @Gs mixin macro adds stomatal conductance fields to a struct"
@mix @with_kw struct Gs 
    gamma::typeof(1.0u"μmol*mol^-1") = 0.0u"μmol*mol^-1" # Gamma for all Ball-Berry type models # 	
    g1::typeof(1.0)                  = 7.0 # Slope parameter 
end

""" 
Ball-Berry stomaratl conductance formulation parameters
`gamma` in μmol*mol^-1 and `g1`scalar, also used for all Ball-Berry type models.
(modelgs = 2 in maestra)
"""
@Gs mutable struct BallBerryStomatalConductance <: AbstractStomatalConductance end

""" 
Leuning stomatal conductance formulation 
Has the extra `d0l` paramater in Pa. 
"""
@Gs mutable struct LeuningStomatalConductance <: AbstractStomatalConductance 
    d0l::typeof(1.0u"Pa")    = 1500.0u"Pa" 
end

""" 
Medlyn stomatal conductance formulation parameters 
Has the extra `vpdmin` paramater in Pa. 
(modelgs = 4 in maestra)
"""
@Gs mutable struct MedlynStomatalConductance <: AbstractStomatalConductance
    vpdmin::typeof(1.0u"Pa") = 1500.0u"Pa"
end

""" Three parameter Ball-Berry stomatal conductance formulation parameters 
Has the extra `vpdmin` paramater in Pa and gk scalar parameter. 
(modelgs = 5 in maestra)
"""
@Gs mutable struct ThreeParStomatalConductance <: AbstractStomatalConductance
    gk::typeof(1.0)          = 0.3 # Medlyn et al. 2011 model
    vpdmin::typeof(1.0u"Pa") = 1500.0u"Pa"
end
 
""" Tuzet stomatal conductance formulation parameters 
(modelgs = 5 in maestra)
"""
@Gs mutable struct TuzetStomatalConductance <: AbstractStomatalConductance end


abstract type AbstractGSShape end

struct HardMinimumGS <: AbstractGSShape end 

@with_kw mutable struct HyperbolicMinimumGS <: AbstractGSShape  
    hmshape::Float64 = 0.5
end


abstract type AbstractPhotoModel end
abstract type AbstractBallBerryModel{GS,SM} <: AbstractPhotoModel end
abstract type AbstractMaespaModel{GS,SM} <: AbstractBallBerryModel{GS,SM} end

" @BB mixin macro adds base Balll-Berry fields to a struct"
@mix @with_kw struct BB
    g0::typeof(1.0u"mol*m^-2*s^-1") = 0.03u"mol*m^-2*s^-1" # Stomatal leakiness (gs when photosynthesis is zero).  
end

" @Maespa mixin macro adds base maespa fields to a struct"
@BB @mix struct Maespa{SH}
    plantk::typeof(1.0u"mmol*m^-2*s^-1*MPa^-1")     = 3.0u"mmol*m^-2*s^-1*MPa^-1" 
    totsoilres::typeof(1.0u"m^2*s^1*MPa^1*mmol^-1") = 0.5u"m^2*s^1*MPa^1*mmol^-1"
    gsshape::SH                                     = HardMinimumGS()
end

"""
General Ball-Berry stomatal conductance model.

THis model has various submodels and soil water methods defined in gsmodel
and soilmethod.
"""
@BB mutable struct BallBerryModel{GS, SM} <: AbstractBallBerryModel{GS,SM}
    gsmodel::GS    = BallBerryStomatalConductance()
    soilmethod::SM = PotentialSoilMethod()
end

"""
Tuzet stomatal conductance model.

This model is limited to using `TuzetStomatalConductance` and `TuzetSoilMethods`.
"""
@Maespa mutable struct TuzetModel{GS<:TuzetStomatalConductance, 
                                           SM<:TuzetSoilMethod{<:AbstractPotentialDependence}
                                          } <: AbstractMaespaModel{GS,SM}
    gsmodel::GS    = TuzetStomatalConductance()
    soilmethod::SM = TuzetSoilMethod()
end


"""
Emax stomatal conductance model from maespa.
Similar options for specialised stomatal conducance, 
but using emax soil water methods.
"""
@Maespa mutable struct EmaxModel{GS<:AbstractStomatalConductance, 
                                          SM<:EmaxSoilMethod{<:AbstractPotentialDependence},
                                         } <: AbstractMaespaModel{GS,SM}
    gsmodel::GS    = BallBerryStomatalConductance()
    soilmethod::SM = EmaxSoilMethod()
end

"""
Fixed parameters for photosynthesis. 

THe majority of these are "composed" from submodels that both hold
the parameters and use their own specific methods during the photosynthesis 
routines.

Calling `PhotoParams()` will give the default values for all of these submodels.
Any parameters and submodels can be overridden with keyword arguments:

`PhotoParams(model=TuzetModel, ca= 450u"μmol*mol^-1")` 
"""
@with_kw mutable struct FvCBPhoto{F<:AbstractPhotoModel,
                                  V<:AbstractVcJmax,
                                  KM<:AbstractCompensation,
                                  Ru<:AbstractRubiscoRegen,
                                  Re<:AbstractRespiration}
    model::F          = BallBerryModel()
    vcjmax::V         = VcJmax()
    compensation::KM  = BernacchiCompensation()
    rubisco_regen::Ru = RubiscoRegen()
    respiration::Re   = Respiration()
end


abstract type AbstractBoundaryConductance end

@with_kw mutable struct BoundaryConductance{Me} <: AbstractBoundaryConductance
    leafwidth::Me = 0.05u"m"
end


abstract type AbstractRadiationConductance end

@with_kw mutable struct YingPingRadiationConductance{T} <: AbstractRadiationConductance
    rdfipt::T = 1.0
    tuipt::T  = 1.0
    tdipt::T  = 1.0
end


abstract type AbstractDecoupling end

mutable struct McNaughtonJarvisDecoupling <: AbstractDecoupling end
struct NoDecoupling <: AbstractDecoupling end


@with_kw mutable struct EnergyBalance{Ph,
                                    Ra<:AbstractRadiationConductance,
                                    Bo<:AbstractBoundaryConductance,
                                    De<:AbstractDecoupling}
    radiation_conductance::Ra = YingPingRadiationConductance()
    boundary_conductance::Bo  = BoundaryConductance()
    decoupling::De            = McNaughtonJarvisDecoupling()
    photo::Ph                 = FvCBPhoto()
    itermax::Int              = 100
end


"@Environment mixin macro"

@mix @with_kw struct Environment
    tair::typeof(1.0u"°C")                  = 25.0u"°C"
    windspeed::typeof(1.0u"m*s^-1")         = 1.0u"m*s^-1"
    par::typeof(1.0*250.0u"μmol*m^-2*s^-1") = 4.575*250.0u"μmol*m^-2*s^-1"
    rnet::typeof(1.0u"W*m^-2")              = 250.0u"W*m^-2"
    soilmoist::typeof(1.0)                  = 0.2
    pressure::typeof(1.0u"kPa")             = 101.25u"kPa"
    tleaf::typeof(1.0u"°C")                 = 25.0u"°C"
    swp::typeof(-100.0u"kPa")               = -100.0u"kPa"
    swpshade::typeof(-100.0u"kPa")          = -100.0u"kPa"
    vpd::typeof(1.0u"kPa")                  = 0.5u"kPa"
    ca::typeof(1.0u"μmol*mol^-1")           = 400.0u"μmol*mol^-1"
    rh::typeof(1.0)                         = 0.5
end


"Vars mixin macro"
@Environment @mix struct Vars
    # shared
    cs::typeof(1.0u"μmol*mol^-1")          = 400.0u"μmol*mol^-1"
    vpdleaf::typeof(1.0u"kPa")             = 0.8u"kPa"
    rhleaf::typeof(1.0)                    = 0.99 # Ball-Berry Stomatal Conductance
    # energy balance
    fheat::typeof(1.0u"W*m^-2")            = 1.0u"W*m^-2"
    gbhu::typeof(1.0u"m*mol*Pa*J^-1*s^-1") = 1.0u"m*mol*Pa*J^-1*s^-1"
    gbhf::typeof(1.0u"m*mol*Pa*J^-1*s^-1") = 1.0u"m*mol*Pa*J^-1*s^-1"
    gh::typeof(1.0u"m*mol*Pa*J^-1*s^-1")   = 1.0u"m*mol*Pa*J^-1*s^-1"
    gradn::typeof(1.0u"mol*m^-2*s^-1")     = 1.0u"mol*m^-2*s^-1"
    lhv::typeof(1.0u"J*mol^-1")            = 1.0u"J*mol^-1"
    et::typeof(1.0u"mol*m^-2*s^-1")        = 1.0u"mol*m^-2*s^-1"
    slope::typeof(1.0u"Pa*K^-1")           = 1.0u"Pa*K^-1"
    decoup::typeof(1.0)                    = 0.0
    # photosynthesis
    gsdiva::typeof(1.0u"mol*μmol^-1")      = 1.0u"mol*μmol^-1"
    km::typeof(1.0u"μmol*mol^-1")          = 1.0u"μmol*mol^-1"
    ci::typeof(1.0u"μmol*mol^-1")          = 1.0u"μmol*mol^-1"
    gammastar::typeof(1.0u"μmol*mol^-1")   = 1.0u"μmol*mol^-1"
    gs::typeof(1.0u"mol*m^-2*s^-1")        = 1.0u"mol*m^-2*s^-1"
    jmax::typeof(1.0u"μmol*m^-2*s^-1")     = 1.0u"μmol*m^-2*s^-1"
    vcmax::typeof(1.0u"μmol*m^-2*s^-1")    = 1.0u"μmol*m^-2*s^-1"
    rd::typeof(1.0u"μmol*m^-2*s^-1")       = 1.0u"μmol*m^-2*s^-1"
    ac::typeof(1.0u"μmol*m^-2*s^-1")       = 1.0u"μmol*m^-2*s^-1"
    aleaf::typeof(1.0u"μmol*m^-2*s^-1")    = 1.0u"μmol*m^-2*s^-1"
    vj::typeof(1.0u"μmol*m^-2*s^-1")       = 1.0u"μmol*m^-2*s^-1"
    aj::typeof(1.0u"μmol*m^-2*s^-1")       = 1.0u"μmol*m^-2*s^-1"
    # soil
    fsoil::typeof(1.0)                     = 1.0
end

"@MaespaVars mixin macro"
@Vars @mix struct MaespaVars
    ktot::typeof(1.0u"mmol*m^-2*s^-1*MPa^-1") = 2.0u"mmol*m^-2*s^-1*MPa^-1" 
    weightedswp::typeof(1.0u"kPa")            = 0.0u"kPa"
    psil::typeof(-111.0u"kPa")                = -111.0u"kPa"
end

@MaespaVars @mix struct TuzetVars
    psilin::typeof(-999.0u"kPa") = -999.0u"kPa"
  	psiv::typeof(-1.9u"kPa")     = -1.9u"kPa"
    sf::typeof(1.0u"kPa^-1")     = 3.2u"kPa^-1"
end

@MaespaVars @mix struct EmaxVars
    minleafwp::typeof(1.0u"kPa")           = 0.1u"kPa"
    emaxleaf::typeof(1.0u"mmol*m^-2*s^-1") = 400.0u"mmol*m^-2*s^-1"
end

@Vars @mix struct JarvisVars
    vmleaf::typeof(1.0u"mmol*mol^-1") = 1.0u"mmol*mol^-1"
end

@Vars mutable struct PhotoVars{} end
@Vars mutable struct BallBerryVars{} end
@EmaxVars mutable struct EmaxVars{} end
@TuzetVars mutable struct TuzetVars{} end
@JarvisVars mutable struct JarvisVars{} end
