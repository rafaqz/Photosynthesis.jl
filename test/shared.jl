using Photosynthesis, Unitful, Test, Libdl, Flatten, Combinatorics

using Unitful: °C, K

const BALLBERRY_GS = 2
const LEUNING_GS = 3
const MEDLYN_GS = 4

const BERNACCI = 0
const BADGERCOLLATZ = 1

const SOILMETHOD_NONE = 0
const SOILMETHOD_EMAX = 1
const SOILMETHOD_VOLUMETRIC = 2
const SOILMETHOD_POTENTIAL = 2
const SOILMETHOD_DEFICIT = 4

const SOILDATA_POTENTIAL = 1     # SOIL MOISTURE DATA IS POTENTIAL (MPa)
const SOILDATA_DEFICIT   = 2     # SOIL MOISTURE DATA IS DEFICIT (DIMNLESS)
const SOILDATA_CONTENT   = 3     # SOIL MOISTURE IS CONTENT (m3 m-3), as when simulated.
const SOILDATA_SIMULATED = 4     # 
const SOILDATA_NONE      = 0     # NO SOIL MOISTURE DATA

# Download and compile maespa and maepa if they're not present
maespa_dir = joinpath(dirname(pathof(Photosynthesis)), "../test/maespa")
maestra_dir = joinpath(dirname(pathof(Photosynthesis)), "../test/maestra")
if !isdir(maespa_dir) 
    run(`git clone --depth 1 https://github.com/rafaqz/maespa.git $maespa_dir`)
    run(`git clone --depth 1 https://github.com/rafaqz/Maestra.git $maestra_dir`)
end

if Sys.islinux()
    if !isfile(joinpath(maespa_dir, "physiol.so"))
        cd(maespa_dir)
        println("Compiling maespa...")
        run(`gfortran -ffree-line-length-200 -shared -O2 maestcom.f90 -o maestcom.so -fPIC`)
        run(`gfortran -ffree-line-length-200 -shared -O2 metcom.f90 -o metcom.so -fPIC`)
        run(`gfortran -shared -O2 physiol.f90 -o physiol.so -fPIC`)
        # cd(maestra_dir)
        # println("Compiling maestra...")
        # run(`gfortran -ffree-line-length-200 -shared -O2 maestcom.f90 -o maestcom.so -fPIC`)
        # run(`gfortran -ffree-line-length-200 -shared -O2 metcom.f90 -o metcom.so -fPIC`)
        # run(`gfortran -shared -O2 physiol.f90 -o physiol.so -fPIC`)
    end
elseif Sys.isapple()
    if !isfile(joinpath(maestra_dir, "physiol.dylib"))
        cd(maespa_dir)
        println("Compiling maespa...")
        run(`gfortran -ffree-line-length-200 -shared -O2 maestcom.f90 -o maestcom.dylib -fpic`) 
        run(`gfortran -ffree-line-length-200 -shared -O2 metcom.f90 -o metcom.dylib -fPIC`)
        run(`gfortran -shared -O2 physiol.f90 -o physiol.dylib -fPIC`)
    end
else
    error("Only linux and apple can be used to test")
end

maespa_photosynlib = dlopen(joinpath(maespa_dir, "physiol"))
# maestra_photosynlib = dlopen(joinpath(maestra_dir, "physiol"))
