###########################################################################################
# Environmental Variables

@chain vars @udefault_kw @units @description 

"""
@MixinEnviroVars mixin macro.
"""
@mix @vars struct MixinEnviroVars{TA,WS,PA,RN,SM,PR,SWP,VPD,CA,RH}
    tair::TA       | (273.15 + 25.0) | K              | _
    windspeed::WS  | 1.0             | m*s^-1         | _
    par::PA        | 4.575*250.0     | μmol*m^-2*s^-1 | _
    rnet::RN       | 250.0           | W*m^-2         | _
    soilmoist::SM  | 0.2             | _              | _
    pressure::PR   | 101250.0        | Pa             | _
    swp::SWP       | -0.1            | MPa            | _
    vpd::VPD       | 500.0           | Pa             | _
    ca::CA         | 400.0           | μmol*mol^-1    | _
    rh::RH         | 0.5             | _              | _
end
