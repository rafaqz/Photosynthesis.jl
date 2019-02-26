###########################################################################################
# Variables

"@Environment mixin macro"
@mix @description @units @udefault_kw struct Environment{K,MS,μMoM2S,WM2,F,kPa,μMoMo}
    tair::K          | (273.15 + 25.0) | K                 | _
    windspeed::MS    | 1.0         | m*s^-1                | _
    par::μMoM2S      | 4.575*250.0 | μmol*m^-2*s^-1        | _
    rnet::WM2        | 250.0       | W*m^-2                | _
    soilmoist::F     | 0.2         | _                     | _
    pressure::kPa    | 101.25      | kPa                   | _
    tleaf::K         | (273.15 + 25.0) | K                 | _
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
    gbh::MMoPaJS     | 1.0         | m*mol*Pa*J^-1*s^-1    | _
    gv::MMoPaJS      | 1.0         | m*mol*Pa*J^-1*s^-1    | _
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
    gsv::MoM2S       | 1.0         | mol*m^-2*s^-1         | _
    gbv::MoM2S       | 1.0         | mol*m^-2*s^-1         | _
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

@Vars mutable struct PhotoVars{} end
