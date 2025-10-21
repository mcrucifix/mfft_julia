using FFTW
using LinearAlgebra

function Q(μ::Float64, ν::Float64, 𝜏::Float64)
    if (μ == ν) 
        return 1
    else 
        return  ( sin(μ-ν)/((μ-ν)*𝜏) * (π^2)/(π^2-((μ-ν)^2 * 𝜏^2)) )
  end
end 




# Function to apply Hanning window, perform FFT, and find the index of the highest component
function find_max_fft_index(input_vector::Vector{Complex{Float64}})
    # Apply Hanning window
    n = length(input_vector)
    # Perform FFT
    fft_result = fft(input_vector)

    # Find the index of the component with the highest absolute value
    max_index = argmax(abs.(fft_result)) - 1

    return (max_index)
end

# Define the function to compute the weighted sum with cis(omega * i)
function proj_cis_omega(input_vector::Vector{Complex{Float64}}, omega::Float64)::Complex{Float64}
    N = length(input_vector)
    sum(input_vector[i] * cis(- omega * (i-1)) for i in 1:(N))
end


# Perform Golden Section Search
function golden_section_search(objective::Function, range = [Float64,Float64],  tol::Float64 = 1e-6)
    φ = (1 + sqrt(5)) / 2  # Golden ratio
    inv_φ = 1 / φ

    a, b = range
    c = b - (b - a) * inv_φ
    d = a + (b - a) * inv_φ

    while abs(b - a) > tol
        if objective(c) > objective(d)
            b = d
        else
            a = c
        end

        # Update c and d
        c = b - (b - a) * inv_φ
        d = a + (b - a) * inv_φ
    end

    # Return the omega that maximizes the absolute value of the weighted sum
    omega_opt = (a + b) / 2
    value_opt = objective(omega_opt)
    return omega_opt, value_opt
end

using Statistics

# V = 0.02 * randn(1024) + im* 0.0004 * randn(1024) .+ 
#              3.0.*(cis(0.142 * t) for t in (0:1023)) .+
#              (2.93 + 0.14im) .*(cis(0.192 * t) for t in (0:1023)) .+
#              (0.93 + 0.54im) .*(cis(0.242 * t) for t in (0:1023)) .+
#              (0.55 + 0.93im) .*(cis(0.441 * t) for t in (0:1023)) 

function mfft(V::Vector{ComplexF64}, M = 5):: Tuple{Vector{ComplexF64}, Vector{Float64}}
    ## function mfft(V::Vector{ComplexF64}, M = 5)
  
  n = length(V)
  𝜏 =  n / 2.0 
  hanning_window = 2. * ( 0.5 .- 0.5 .* cos.(2.0 .* π .* (0:n-1) ./ n) )
  
  fm = V

  ν = Vector{Float64}(undef, M)
  F = Vector{ComplexF64}(undef, M)
  A = Vector{ComplexF64}(undef, M)
  α = zeros(ComplexF64, M, M)
  B = zeros(ComplexF64, M, M)
  
  
  α[1,1] = 1.0 + 0.0im 
  
  for m in 1:M
    input_vector  = fm .* hanning_window
    
    # look for frequency of maximum power
    maxI = find_max_fft_index(input_vector)
    objective = omega -> abs(proj_cis_omega(input_vector, omega))
    range = 2 * π * ( maxI .+  [-0.5, 0.5]) / n 
    
    omega, amp = golden_section_search(objective, range) 
    F[m] = proj_cis_omega( input_vector, omega ) / n
    
    #  fm1 = V .- Fm .* (cis( omega * i) for i in 0:(n-1))
    # eq. 21
    
    ν[m] = omega
   
    if (m == 1)
        α[m,m] = 1.0 
    else
        for j in (1:m)
            B[j,m] = - sum( α[j,s] * Q(ν[m], ν[s], 𝜏) for s in 1:j)
        end
        α[m,m] = 1.0/ ( sqrt( 1 - sum( abs2(B[j,m]) for j in 1:(m-1) ) ) )
        for j in (m-1):-1:1
            α[m,j] = α[m,m] * sum( B[s,m] * α[s,j] for s in j:(m-1))
        end
    end
    
    for j in 1:m
        fm =  fm .- (α[m,m] * F[m]) .* 
             α[m,j] .* (cis.( ν[j] * t for t in 0:(n-1))) * cis( ν[m] - ν[j] ) 
     end 
  
  end 
  
  # final amplitudes 
  
   for s in 1:M
       A[s] = sum( α[j,j]*α[j,s] * F[j] * cis( ν[j] - ν[s] ) for j in (s:M))
   end
  
   return A, ν

end
