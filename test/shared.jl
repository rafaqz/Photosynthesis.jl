using Photosynthesis, Unitful, Test, Libdl

using Unitful: °C, K

const BALLBERRY_GS = 2
const LEUNING_GS = 3
const MEDLYN_GS = 4

const BADGERCOLLATZ = 1
const BERNACCI = 0

const SOILDATA_POTENTIAL = 1     # SOIL MOISTURE DATA IS POTENTIAL (MPa)
const SOILDATA_DEFICIT   = 2     # SOIL MOISTURE DATA IS DEFICIT (DIMNLESS)
const SOILDATA_CONTENT   = 3     # SOIL MOISTURE IS CONTENT (m3 m-3), as when simulated.
const SOILDATA_SIMULATED = 4     # 
const SOILDATA_NONE      = 0     # NO SOIL MOISTURE DATA

photosynlib = dlopen(joinpath(ENV["MAESPA"], "physiol"))
