struct Lower end
struct Upper end

"""
Solves quadratic equations
"""
function quad(typ, a, b, c)
    x = (b * b - 4 * a * c) 
    x < zero(x) && error("imaginary roots in quadratic")
    if a == zero(a) 
        if b == zero(b)
            q = zero(-c / b)
            c != zero(c) && error("error: cant solve quadratic")
        else
            q = -c / b
        end
    else 
       q = side(typ, x, a, b)
    end
    q
end

@inline side(::Upper, x, a, b) = (-b + sqrt(x)) / 2a
@inline side(::Lower, x, a, b) = (-b - sqrt(x)) / 2a
