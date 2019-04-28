###########################################################################################
# Environmental Variables

"@Environment mixin macro"
@mix @description @units @udefault_kw struct MixinEnviroVars{K,MS,μMoM2S,WM2,F,kPa,μMoMo}
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
