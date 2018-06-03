
# Create the subset {x0,f(x0),...,f^N(x0)} of an orbit...
function createOrbit(f::Function, x0::Real, N::Integer=10)
  orbit = [x0]
  
  for i in 1:N
    push!(orbit, f(orbit[i-1]))
  end
  
  orbit
end


# A Newton-Raphson method implementation
function findRoot(f::Function, f_der::Function, x0::Real, N::Integer=100, ɛ::Real=0.001)

  x = x0
  y = f(x0)
  iters = 0

  while abs(y) > ɛ && iters < N

    y_der = f_der(x)

    if abs(y_der) < ɛ
      return Inf
    end

    x -= y/y_der
    y = f(x)

    iters += 1

  end # while

  x

end # function


# Numerical (central) approximation to derivative
function numDer(f::Function, ɛ::Real=0.001)
  x -> (f(x+ɛ)-f(x-ɛ))/(2ɛ)
end


# A Newton-Raphson method implementation, with numerical approximation to derivative
function findRootNumDer(f::Function, x0::Real, N::Integer=100, ɛ::Real=0.001)
  findRoot(f, numDer(f), x0, N, ɛ)
end # function


# Create a piecewise linear function with points of period n odd, but without
# points of period m<n odd. Lets take n=2k+1...
function createPWLOdd(k::Integer)
  function pwlOdd(x::Real)
    if x < k
      return -x+2k+2
    elseif x < k+1
      return -2x+3k+2
    elseif x < 2k
      return -x+2k+1
    else
      return k*x-2k^2+1
    end # if elseif else
  end # function
end # function


# Create the double of a function
function createDouble(f::Function)
  function F(x::Real)
    if x < 1/3
      return f(3x)/3 + 2/3
    elseif x < 2/3
      return (2.+f(1.))*(2/3-x)
    else
      return x - 2/3
    end # if elseif else
  end # function
end # function


