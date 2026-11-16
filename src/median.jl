# Spatial (geometric) median: the point minimising the sum of Euclidean
# distances to the rows of P. Weiszfeld's iteration with the Vardi & Zhang
# (2000) modification, which stays correct when the median falls on a data
# point. Started from the coordinate-wise median, iterated to convergence.

function spatial_median(P; tol::Float64=1e-12, maxit::Int=10_000)
    n, k = size(P)
    k == 0 && return Float64[]
    n == 1 && return vec(P[1, :])
    y = [median(view(P, :, j)) for j in 1:k]
    T = zeros(k)
    R = zeros(k)
    scale = max(maximum(abs, P), eps())
    for _ in 1:maxit
        fill!(T, 0.0)
        wsum = 0.0
        ties = 0
        for i in 1:n
            d = 0.0
            for j in 1:k
                d += (P[i, j] - y[j])^2
            end
            d = sqrt(d)
            if d <= 1e-14 * scale
                ties += 1
                continue
            end
            w = 1 / d
            wsum += w
            for j in 1:k
                T[j] += w * P[i, j]
            end
        end
        wsum == 0 && return y                       # all points coincide with y
        T ./= wsum                                 # Weiszfeld step
        if ties == 0
            ynew = T
        else
            # Vardi-Zhang: blend towards the data point the iterate sits on
            for j in 1:k
                R[j] = (T[j] - y[j]) * wsum
            end
            r = sqrt(sum(abs2, R))
            r <= ties && return y                  # y is the median
            γ = ties / r
            ynew = (1 - γ) .* T .+ γ .* y
        end
        step = sqrt(sum(abs2, ynew .- y))
        y = copy(ynew)
        step <= tol * scale && break
    end
    return y
end
