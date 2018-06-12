"""
Solves quadratic equations - finds smaller root.
"""
function quad(typ, a, b, c)
    x = (b * b - 4.0 * a * c) 
    x < zero(x) && error("imaginary roots in quadratic")
    if a == zero(a) 
        if b == zero(b)
            q = 0.0
            c != zero(c) && error("error: cant solve quadratic")
        else
            q = -c / b
        end
    else 
       q = side(typ, x, a, b)
    end
    q
end

@inline side(::Type{Val{:upper}}, x, a, b) = (-b + sqrt(x)) / (2.0 * a)
@inline side(::Type{Val{:lower}}, x, a, b) = (-b - sqrt(x)) / (2.0 * a)
