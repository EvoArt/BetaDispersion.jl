
using Turing

function param(α,β)
    μ = α/(α+β)
    σ² = (α*β)/((α+β)^2 *(α+β+1))
    return [μ σ²]
end
    @model function BayesDisp(R,μ_scale = 1)
        k = length(R)      
        μ ~ filldist(truncated(Normal(0,μ_scale),0,Inf),k)
        σ ~ filldist(Exponential(1),k)
        for i in 1:k
            R[i] .~ Normal(μ[i],σ[i])
        end
    end

    @model function BetaDisp(R)
        k = length(R) 
        α ~ filldist(Uniform(1,100),k)#Gamma(2,2)
        β ~ filldist(Uniform(1,100),k)#Gamma(2,2)
        for i in 1:k
            R[i] .~ Beta(α[i],β[i])
        end     
    end
function bayesdisp(disp,μ_scale,n_samples = 1000,n_chains = 4, sampler = NUTS())
    mod = BayesDisp(disp.residuals, μ_scale)
    chns = sample(mod,sampler, MCMCThreads(),n_samples,n_chains)
    D = Dict()
    for (i,level) in enumerate(string.(disp.levels))
        D["μ[$(i)]"] = level * "_μ"
        D["σ[$(i)]"] = level * "_σ"
    end
    
    return replacenames(chns, D)
end

function bayesdisp(disp,n_samples = 1000,n_chains = 4, sampler = NUTS())
    mod = BetaDisp(disp.residuals)
    chns = sample(mod,sampler, MCMCThreads(),n_samples,n_chains)
    D = Dict()
    for (i,level) in enumerate(string.(disp.levels))
        D["α[$(i)]"] = level * "_μ"
        D["β[$(i)]"] = level * "_σ"
    end
    
    return replacenames(chns, D)
end

export bayesdisp, param

    
