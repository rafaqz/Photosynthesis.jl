
abstract type AbstractPhotosynthesis end

abstract type AbstractPhotoModel end

@mix @columns struct PhotoModel{V<:AbstractVcJmax,
                                KM<:AbstractCompensation,
                                Ru<:AbstractRubiscoRegen,
                                Re<:Union{Nothing,AbstractRespiration}
                               }
    vcjmax::V         | VcJmax()                | _ | _ | _ | _
    compensation::KM  | BernacchiCompensation() | _ | _ | _ | _
    rubisco_regen::Ru | RubiscoRegen()          | _ | _ | _ | _
    respiration::Re   | Respiration()           | _ | _ | _ | _       
end
