
abstract type AbstractSoilData end
abstract type AbstractDeficitSoilData <: AbstractSoilData end
abstract type AbstractContentSoilData <: AbstractDeficitSoilData end

@chain columns @description @prior @units @default_kw

###########################################################################################
# Parameters

" @Deficit mixin macro adds water deficit fields to a struct"
@mix @columns struct Deficit{T}
    smd1::T   | 1.0 | _      | _ | _
    smd2::T   | 1.0 | _      | _ | _
end

"@Content mixin macro adds water content fields to a struct"
@Deficit @mix struct Content{W}
    swmax::W  | 1.0 | _      | _ | _
    swmin::W  | 0.0 | _      | _ | _
end

"Soil water is measured as deficit"
@Deficit struct DeficitSoilData{} <: AbstractDeficitSoilData end
"Soil water is measured as volumetric content"
@Content struct ContentSoilData{} <: AbstractContentSoilData end
"Simulated soil volumetric content"
@Content struct SimulatedSoilData{} <: AbstractContentSoilData end
"Soil water is measured as water potential"
@columns struct PotentialSoilData{T} <: AbstractSoilData
    swpexp::T | 1.0 | _      | _ | _
end
struct NoSoilData <: AbstractSoilData end


abstract type AbstractPotentialDependence end

"@Potdep mixin macro adds potential dependence fields to a struct"
@mix @columns struct Potdep{Pa}
    vpara::Pa | 1.0 | kPa | Gamma(10, 1/10) | _
    vparb::Pa | 2.0 | kPa | Gamma(10, 2/10) | _
end
      
@units struct Potdep{Pa}
    vpara::Pa | kPa
    vparb::Pa | kPa
end

@Potdep struct LinearPotentialDependence{} <: AbstractPotentialDependence end
@Potdep struct ZhouPotentialDependence{} <: AbstractPotentialDependence end
struct NoPotentialDependence{} <: AbstractPotentialDependence end   


abstract type AbstractSoilMethod{T} end

@columns struct VolumetricSoilMethod{WC,R,D,T<:AbstractContentSoilData
                                            } <: AbstractSoilMethod{T}
    soildata::T  | ContentSoilData() | _ | _ | _
    wc1::WC      | 0.5               | _ | Beta(2.0, 2.0) | _
    wc2::WC      | 0.5               | _ | Beta(2.0, 2.0) | _
    soilroot::R  | 0.5               | _ | _ | _ # Not sure these are parameters
    soildepth::D | 0.5               | _ | _ | _ # One must be a variable...
end

@default_kw struct ConstantSoilMethod{T<:NoSoilData} <: AbstractSoilMethod{T}
    soildata::T | NoSoilData()
end
@default_kw struct DeficitSoilMethod{T<:AbstractDeficitSoilData} <: AbstractSoilMethod{T}
    soildata::T | ContentSoilData()
end
@default_kw struct PotentialSoilMethod{T<:PotentialSoilData} <: AbstractSoilMethod{T}
    soildata::T | PotentialSoilData()
end

@mix @default_kw struct MaespaSoil{T<:AbstractPotentialDependence}
    soildata::T | LinearPotentialDependence()
end

@MaespaSoil struct EmaxSoilMethod{} <: AbstractSoilMethod{T} end
@MaespaSoil struct TuzetSoilMethod{} <: AbstractSoilMethod{T} end


"Physiological constants for CO2 snd rubisco comprensation points"
abstract type AbstractCompensation end

@columns struct BadgerCollatzCompensation{μMoMo,kJMo,C} <: AbstractCompensation
    Kc25::μMoMo     | 404.0    | μmol*mol^-1    | Gamma(100, 404/100)    | "MM coefft of Rubisco for CO2"
    Ko25::μMoMo     | 248000.0 | μmol*mol^-1    | Gamma(100, 248000/100) | "MM coefft of Rubisco for O2"
    ΔHa_Kc::kJMo    | 59.4     | kJ*mol^-1      | Gamma(100, 59.4/100)   | "Temp. response of Kc"
    ΔHa_Ko::kJMo    | 36.0     | kJ*mol^-1      | Gamma(100, 36/100)     | "Temp. response of Ko"
    tref::C         | 25.0     | °C             | _                      | "Temperature reference, usually 25.0° C"
end

"""
Bernacchi et al 2001 PCE 24: 253 - 260
Extra deactivation terms may be required above 40°C.
"""
@columns struct BernacchiCompensation{μMoMo,kJMo,C} <: AbstractCompensation
    Kc25::μMoMo     | 404.9    | μmol*mol^-1    | Gamma(10, 404.9/10)  | "MM coefft of Rubisco for CO2"
    Ko25::μMoMo     | 278400.0 | μmol*mol^-1    | Gamma(10, 278400/10) | "MM coefft of Rubisco for O2"
    Γ☆25::μMoMo     | 42.75    | μmol*mol^-1    | Gamma(10, 42.75/10)  | _
    ΔHa_Kc::kJMo    | 79.43    | kJ*mol^-1      | Gamma(10, 79.43/10)  | "Temp. response of Kc"
    ΔHa_Ko::kJMo    | 36.38    | kJ*mol^-1      | Gamma(10, 36.38/10)  | "Temp. response of Ko"
    ΔHa_Γ☆::kJMo    | 37.83    | kJ*mol^-1      | Gamma(10, 37.83/10)  | _
    tref::C         | 25.0     | °C             | _                    | _
end

"Electron flux parameters"
abstract type AbstractJmax end
 
@columns struct Jmax{μMoM2S,JMoK,JMo} <: AbstractJmax
    jmax25::μMoM2S  | 184.0    | μmol*m^-2*s^-1 | Gamma(100, 184/100)    | _
    delsj::JMoK     | 640.02   | J*mol^-1*K^-1  | Gamma(100, 640.02/100) | "DELTAS in Medlyn et al. (2002)"
    eavj::JMo       | 37259.0  | J*mol^-1       | Gamma(100, 37259/100)  | "Ha     in Medlyn et al. (2002)"
    edvj::JMo       | 200000.0 | J*mol^-1       | Gamma(100, 200000/100) | "Hd     in Medlyn et al. (2002)"
end

abstract type AbstractVcmax end

"@Vcmax mixin macro adds vcmax fields to a struct"
@mix @columns struct Vcmax{μMoM2S,JMo}
    vcmax25::μMoM2S | 110.0    | μmol*m^-2*s^-1 | Gamma(100, 114/100)   | _
    eavc::JMo       | 47590.0  | J*mol^-1       | Gamma(100, 47590/100) | "Ha    in Medlyn et al. (2002)"
end

@Vcmax struct NoOptimumVcmax{} <: AbstractVcmax end
@Vcmax struct OptimumVcmax{JMo,JMoK} <: AbstractVcmax
    edvc::JMo       | 1.0      | J*mol^-1       | Gamma(10, 1/10)        | "Hd in Medlyn et al. (2002)"
    delsc::JMoK     | 629.26   | J*mol^-1*K^-1  | Gamma(100, 629.26/100) | "DELTAS in Medlyn et al. (2002)"
end


abstract type AbstractVcJmax end

@default_kw struct VcJmax{J,V} <: AbstractVcJmax where {J<:AbstractVcmax,V<:AbstractJmax}
    jmaxformulation::J  | Jmax()
    vcmaxformulation::V | NoOptimumVcmax()
end
@columns struct DukeVcJmax{J,V,C} <: AbstractVcJmax where {J<:AbstractVcmax,V<:AbstractJmax}
    jmaxformulation::J  | Jmax()           | _  | _ | _
    vcmaxformulation::V | NoOptimumVcmax() | _  | _ | _
    tvjup::C            | 10.0             | °C | _ | _
    tvjdn::C            | 0.0              | °C | _ | _
end


abstract type AbstractRubiscoRegen end
@columns struct RubiscoRegen{} <: AbstractRubiscoRegen
    theta::Float64  | 0.4    | _              | Beta(8, 12)  | "Shape parameter of the non-rectangular hyperbola."
    ajq::Float64    | 0.324  | _              | Beta(4, 8.3) | "Quantum yield of electron transport."
end

abstract type AbstractRespiration end
@mix @columns struct Resp{pK,C,F,μMoM2S}
    q10f::pK        | 0.67   | K^-1           | _                 | "logarithm of the Q10 (Equation for respiration)"
    dayresp::F      | 1.0    | _              | Beta(5, 1)        | "Respiration in the light as fraction of that in the dark."
    rd0::μMoM2S     | 0.9    | μmol*m^-2*s^-1 | Gamma(10, 0.9/10) | _
    tbelow::C       | -100.0 | °C             | _                 | "No respiration occurs below this temperature."
    tref::C         | 25.0   | °C             | _                 | "Reference temperature (T at which RD0 was measured)"
end

@Resp struct Respiration{} <: AbstractRespiration end
@Resp struct AcclimatizedRespiration{pK,K} <: AbstractRespiration
    k10f::pK        | 0.0    | K^-1           | _             | _
    tmove::K        | 1.0    | K              | Gamma(2, 1/2) | _
end

abstract type AbstractStomatalConductance end
" @Gs mixin macro adds stomatal conductance fields to a struct"
@mix @columns struct Gs{μMoMo,F}
    gamma::μMoMo    | 0.0    | μmol*mol^-1    | Gamma(1, 2)     | "Gamma for all Ball-Berry type models"
    g1::F           | 7.0    | _              | Gamma(10, 7/10) | "Slope parameter"
end

"""
Ball-Berry stomaratl conductance formulation parameters
`gamma` in μmol*mol^-1 and `g1`scalar, also used for all Ball-Berry type models.
(modelgs = 2 in maestra)
"""
@Gs struct BallBerryStomatalConductance{} <: AbstractStomatalConductance end

"""
Leuning stomatal conductance formulation
Has the extra `d0l` paramater in Pa.
"""
@Gs struct LeuningStomatalConductance{Pa} <: AbstractStomatalConductance
    d0l::Pa    | 1500.0 | Pa | Gamma(10, 1500/10) | _
end

"""
Medlyn stomatal conductance formulation parameters
Has the extra `vpdmin` paramater in Pa.
(modelgs = 4 in maestra)
"""
@Gs struct MedlynStomatalConductance{Pa} <: AbstractStomatalConductance
    vpdmin::Pa | 1500.0 | Pa | Gamma(10, 1500/10) | _
end

""" Three parameter Ball-Berry stomatal conductance formulation parameters
Has the extra `vpdmin` paramater in Pa and gk scalar parameter.
(modelgs = 5 in maestra)
"""
@Gs struct ThreeParStomatalConductance{F,Pa} <: AbstractStomatalConductance
    gk::F      | 0.3    | _  | Gamma(2, 0.3/2)    | _
    vpdmin::Pa | 1500.0 | Pa | Gamma(10, 1500/10) | _
end

""" Tuzet stomatal conductance formulation parameters
(modelgs = 5 in maestra)
"""
@Gs struct TuzetStomatalConductance{} <: AbstractStomatalConductance end


abstract type AbstractGSShape end

struct HardMinimumGS <: AbstractGSShape end

@columns struct HyperbolicMinimumGS <: AbstractGSShape
    hmshape::Float64 | 0.999 | _ | Beta(20, 1.01) | _
end


abstract type AbstractPhotoModel end
abstract type AbstractBallBerryModel{GS,SM} <: AbstractPhotoModel end
abstract type AbstractMaespaModel{GS,SM} <: AbstractBallBerryModel{GS,SM} end

" @BB mixin macro adds base Balll-Berry fields to a struct"
@mix @columns struct BB{MoM2S}
    g0::MoM2S | 0.03 | mol*m^-2*s^-1 | Gamma(10, 0.03/10) | "Stomatal leakiness (gs when photosynthesis is zero)."
end

" @Maespa mixin macro adds base maespa fields to a struct"
@BB @mix struct Maespa{mMoM2S,M2SMPaMo,SH}
    plantk::mMoM2S       | 3.0             | mmol*m^-2*s^-1*MPa^-1 | _ | _
    totsoilres::M2SMPaMo | 0.5             | m^2*s^1*MPa^1*mmol^-1 | _ | _
    gsshape::SH          | HardMinimumGS() | _                     | _ | _
end

"""
General Ball-Berry stomatal conductance model.

THis model has various submodels and soil water methods defined in gsmodel
and soilmethod.
"""
@BB struct BallBerryModel{GS,SM} <: AbstractBallBerryModel{GS,SM}
    gsmodel::GS    | BallBerryStomatalConductance() | _ | _ | _
    soilmethod::SM | PotentialSoilMethod()          | _ | _ | _
end

"""
Tuzet stomatal conductance model.

This model is limited to using `TuzetStomatalConductance` and `TuzetSoilMethods`.
"""
@Maespa struct TuzetModel{GS<:TuzetStomatalConductance,
                          SM<:TuzetSoilMethod{<:AbstractPotentialDependence}
                         } <: AbstractMaespaModel{GS,SM}
    gsmodel::GS    | TuzetStomatalConductance() | _ | _ | _
    soilmethod::SM | TuzetSoilMethod()          | _ | _ | _
end


"""
Emax stomatal conductance model from maespa.
Similar options for specialised stomatal conducance,
but using emax soil water methods.
"""
@Maespa struct EmaxModel{GS<:AbstractStomatalConductance,
                         SM<:EmaxSoilMethod{<:AbstractPotentialDependence},
                        } <: AbstractMaespaModel{GS,SM}
    gsmodel::GS    | BallBerryStomatalConductance() | _ | _ | _
    soilmethod::SM | EmaxSoilMethod()               | _ | _ | _
end

"""
Fixed parameters for photosynthesis.

THe majority of these are "composed" from submodels that both hold
the parameters and use their own specific methods during the photosynthesis
routines.

Calling `PhotoParams()` will give the default values for all of these submodels.
Any parameters and submodels can be overridden with keyword arguments:

`PhotoParams(model=TuzetModel, ca= 450 | μmol*mol^-1)`
"""
@default_kw struct FvCBPhoto{F<:AbstractPhotoModel,
                             V<:AbstractVcJmax,
                             KM<:AbstractCompensation,
                             Ru<:AbstractRubiscoRegen,
                             Re<:AbstractRespiration}
    model::F          | BallBerryModel()       
    vcjmax::V         | VcJmax()               
    compensation::KM  | BernacchiCompensation()
    rubisco_regen::Ru | RubiscoRegen()         
    respiration::Re   | Respiration()          
end


abstract type AbstractBoundaryConductance end

@columns struct BoundaryConductance{Me} <: AbstractBoundaryConductance
    leafwidth::Me | 0.05 | m | Gamma(2, 0.05/2) | "Mean width of leaves"
end


abstract type AbstractRadiationConductance end

@columns struct YingPingRadiationConductance{T} <: AbstractRadiationConductance
    rdfipt::T | 1.0 | _ | _ | _
    tuipt::T  | 1.0 | _ | _ | _
    tdipt::T  | 1.0 | _ | _ | _
end


abstract type AbstractDecoupling end

struct McNaughtonJarvisDecoupling <: AbstractDecoupling end
struct NoDecoupling <: AbstractDecoupling end


@flattenable @default_kw struct EnergyBalance{Ph,
                                    Ra<:AbstractRadiationConductance,
                                    Bo<:AbstractBoundaryConductance,
                                    De<:AbstractDecoupling}
    radiation_conductance::Ra | YingPingRadiationConductance() | true
    boundary_conductance::Bo  | BoundaryConductance()          | true
    decoupling::De            | McNaughtonJarvisDecoupling()   | true
    photo::Ph                 | FvCBPhoto()                    | true
    itermax::Int              | 100                            | false
end



###########################################################################################
# Variables

"@Environment mixin macro"
@mix @description @units @default_kw struct Environment{C,MS,μMoM2S,WM2,F,kPa,μMoMo}
    tair::C          | 25.0        | °C                    | _
    windspeed::MS    | 1.0         | m*s^-1                | _
    par::μMoM2S      | 4.575*250.0 | μmol*m^-2*s^-1        | _
    rnet::WM2        | 250.0       | W*m^-2                | _
    soilmoist::F     | 0.2         | _                     | _
    pressure::kPa    | 101.25      | kPa                   | _
    tleaf::C         | 25.0        | °C                    | _
    swp::kPa         | -100.0      | kPa                   | _
    swpshade::kPa    | -100.0      | kPa                   | _
    vpd::kPa         | 0.5         | kPa                   | _
    ca::μMoMo        | 400.0       | μmol*mol^-1           | _
    rh::F            | 0.5         | _                     | _
end

"Vars mixin macro"
@Environment @mix struct Vars{μMoMo,kPa,F,WM2,MMoPaJS,μMoM2S,JMo,PaK,MoM2S,MoμMo,μMoMo}
    # shared
    cs::μMoMo        | 400.0       | μmol*mol^-1           | _
    vpdleaf::kPa     | 0.8         | kPa                   | _
    rhleaf::F        | 0.99        | _                     |  "Only in Ball-Berry Stomatal Conductance"
    # energy balance
    fheat::WM2       | 1.0         | W*m^-2                | _
    gbhu::MMoPaJS    | 1.0         | m*mol*Pa*J^-1*s^-1    | _
    gbhf::MMoPaJS    | 1.0         | m*mol*Pa*J^-1*s^-1    | _
    gh::MMoPaJS      | 1.0         | m*mol*Pa*J^-1*s^-1    | _
    gradn::MoM2S     | 1.0         | mol*m^-2*s^-1         | _
    lhv::JMo         | 1.0         | J*mol^-1              | _
    et::MoM2S        | 1.0         | mol*m^-2*s^-1         | _
    slope::PaK       | 1.0         | Pa*K^-1               | _
    decoup::F        | 0.0         | _                     | _
    # photosynthesis
    gsdiva::MoμMo    | 1.0         | mol*μmol^-1           | _
    km::μMoMo        | 1.0         | μmol*mol^-1           | _
    ci::μMoMo        | 1.0         | μmol*mol^-1           | _
    gammastar::μMoMo | 1.0         | μmol*mol^-1           | _
    gs::MoM2S        | 1.0         | mol*m^-2*s^-1         | _
    jmax::μMoM2S     | 1.0         | μmol*m^-2*s^-1        | _
    vcmax::μMoM2S    | 1.0         | μmol*m^-2*s^-1        | _
    rd::μMoM2S       | 1.0         | μmol*m^-2*s^-1        | _
    ac::μMoM2S       | 1.0         | μmol*m^-2*s^-1        | _
    aleaf::μMoM2S    | 1.0         | μmol*m^-2*s^-1        | _
    vj::μMoM2S       | 1.0         | μmol*m^-2*s^-1        | _
    aj::μMoM2S       | 1.0         | μmol*m^-2*s^-1        | _
    # soil
    fsoil::F         | 1.0         | _                     | _
end

"@MaespaVars mixin macro"
@Vars @mix struct MaespaVars{mMoM2SMPa,kPa}
    ktot::mMoM2SMPa  | 2.0         | mmol*m^-2*s^-1*MPa^-1 | _
    weightedswp::kPa | 0.0         | kPa                   | _
    psil::kPa        | -111.0      | kPa                   | _
end

@MaespaVars mutable struct TuzetVars{kPa,pkPa}
    psilin::kPa      | -999.0      | kPa                   | _
    psiv::kPa        | -1.9        | kPa                   | _
    sf::pkPa         | 3.2         | kPa^-1                | _
end

@MaespaVars mutable struct EmaxVars{kPa,mMoM2S}
    minleafwp::kPa   | 0.1         | kPa                   | _
    emaxleaf::mMoM2S | 400.0       | mmol*m^-2*s^-1        | _
end

@Vars mutable struct JarvisVars{mMoMo}
    vmleaf::mMoMo    | 1.0         | mmol*mol^-1           | _
end

@Vars mutable struct PhotoVars{} end
@Vars mutable struct BallBerryVars{} end
