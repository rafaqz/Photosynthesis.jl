using Documenter, DocStringExtensions, Photosynthesis

using Photosynthesis: rubisco_limited_rate,
                      latent_heat_water_vapour,
                      max_rubisco_activity,
                      psil_gs!,
                      calc_decoupling,
                      converge_tleaf!,
                      calc_fco2,
                      calc_flight,
                      quad,
                      forced_boundary_conductance,
                      transport_limited_rate,
                      extremes!,
                      rubisco_compensation_point,
                      soil_water_conductance,
                      assimilation_gs_unknown!,
                      co2_compensation_point,
                      soil_soil_water_conductance,
                      model_init!,
                      max_electron_transport_rate,
                      leaf_water_potential_finder,
                      arrhenius,
                      shape_gs,
                      yingping_radiation_conductance,
                      free_boundary_conductance,
                      grashof_number,
                      assimilation_gs_known!,
                      model_update!,
                      fpsil,
                      soil_water_conductance!,
                      calc_fvpd,
                      calc_ftemp,
                      rubisco_regeneration,
                      vjmax_water

@template (FUNCTIONS, METHODS, MACROS) =
    """
    $(SIGNATURES)
    $(DOCSTRING)
    $(METHODLIST)
    """

@template (TYPES,) =
    """
    $(TYPEDEF)
    $(DOCSTRING)
    $(FIELDS)
    """


makedocs(
    modules = [Photosynthesis],
    doctest = false,
    clean = false,
    sitename = "Photosynthesis.jl",
    format = :html,
    pages = Any[
        "Introduction" => "index.md",
    ]
)

# deploydocs(
#     repo = "github.com/rafaqz/Photosynthesis.jl.git",
#     osname = "linux",
#     julia = "0.6",
#     target = "build",
#     deps = nothing,
#     make = nothing
# )
