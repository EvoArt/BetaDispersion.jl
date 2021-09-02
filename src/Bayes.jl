
using Turing
    @model function BayesDisp(R,μ_scale = 1)
        k = length(R)      
        μ ~ filldist(truncated(Normal(0,μ_scale),0,Inf),k)
        σ ~ filldist(Exponential(1),k)
        for i in 1:k
            R[i] .~ Normal(μ[i],σ[i])
        end
    end


function bayesdisp(disp,μ_scale = 1.0,n_samples = 1000,n_chains = 4, sampler = NUTS())
    mod = BayesDisp(disp.residuals, μ_scale)
    chns = sample(mod,sampler, MCMCThreads(),n_samples,n_chains)
    D = Dict()
    for (i,level) in enumerate(string.(disp.levels))
        D["μ[$(i)]"] = level * "_μ"
        D["σ[$(i)]"] = level * "_σ"
    end
    
    return replacenames(chns, D)
end

export bayesdisp

    
