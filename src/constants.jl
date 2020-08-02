const TOL = 0.02K         # Tolerance for leaf temp iteration
const GBHGBC = 1.32       # Ratio of Gbh:Gbc
const GSVGSC = 1.57       # Ratio of Gsw:Gsc
const GBVGBH = 1.075      # Ratio of Gbw:Gbh
const OI = 205000μmol*mol^-1 # Oxygen partial pressure (umol mol-1)
const CPAIR = 1010.0J*kg^-1*K^-1  # heat capacity of air (J kg-1 K-1)
const AIRMA = 29.0g*mol^-1        # mol mass air (kg / mol)
const CPAIR = 1010.0J*kg^-1*K^-1  # heat capacity of air (J kg-1 K-1)
const EMLEAF = 0.95       # Emissivity of thermal radiation by leaf TODO: should not be constant, and is often higher
const H2OLV0 = 2.501e6J*kg^-1    # latent heat H2O (J/kg)
const H2OMW = 18.0e-3kg*mol^-1# mol mass H2O (kg/mol)
const DHEAT = 21.5e-6m^2*s^-1 # molecular diffusivity for heat


# Values of physical constants
# const DEFWIND = 2.5m*s^-1       # Default wind speed (m s-1)
# const UMOLPERJ = 4.57μmol*J^-1     # Conversion from J to umol quanta
# const FPAR = 0.5          # Fraction of global radiation that is PAR
# const ABSZERO = -273.15   # Absolute zero in degrees Celsius
# const FREEZE = 273.15     # Zero degrees Celsius in Kelvin
# const TAU = 0.76          # Transmissivity of atmosphere
# const PID2 = pi / 2.0     # Pi divided by two
# const PID180 = pi / 180.0 # Pi divided by 180 degrees
# const AIRMA = 29.0g*mol^-1   # mol mass air (kg / mol)
# const PATM = 1.0125e5Pa     # atmospheric pressure - standard condns (Pa)
# const CPH2O = 4.186e06J*kg^-1*K^-1    # heat capacity of water (J kg-1 K-1)
# const CPQUARTZ = 1.942e06J*kg^-1*K^-1 # heat capacity of quartz (J kg-1 K-1)
# const TCQUARTZ = 7.7W*m^-1*K^-1      # thermal conductivity of quartz (W m-1 K-1)
# const TCH2O = 0.594W*m^-1*K^-1       # thermal conductivity of water (W m-1 K-1)
# const TCORG = 0.25W*m^-1*K^-1        # thermal conductivity of organic matter (W m-1 K-1)
# const EMLEAF = 0.95       # Emissivity of thermal radiation by leaf
# const EMSOIL = 0.95       # Emissivity of thermal radiation by soil
# const H2OLV0 = 2.501e6J*kg^-1    # latent heat H2O (J/kg)
# const H2OMW = 18.0e-3kg*mol^-1# mol mass H2O (kg/mol)
# const H2OVW = 18.05e-6m^3*mol^-1    # partial molal volume of water at 20C (m3 mol-1)
# const RCONST = 8.314J*mol^-1*K^-1      # universal gas constant (J/mol/K)
# const SIGMA = 5.67e-8W*m^-2*K^-4     # Steffan Boltzman constant (W/m2/K4)
# const ALPHAQ = 0.425mol*mol^-1      # Quantum yield of RuBP regen (mol mol-1)
# const SOLARC = 1370J*m^-2*s^-1       # Solar constant (J m-2 s-1)
# const GCPERMOL = 12.0     # grams # per mol #
# const CPERDW = 0.5        # fraction per DW
# const VONKARMAN = 0.41    # von Karman"s constant
# const EMLEAF = 0.95       # Emissivity of thermal radiation by leaf
# const H2OLV0 = 2.501e6J*kg^-1    # latent heat H2O (J/kg)
# const H2OMW = 18.0e-3kg*mol^-1# mol mass H2O (kg/mol)
