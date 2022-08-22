function mv_mul!(r::AbstractArray, A::AbstractArray, v::AbstractArray)
    fill!(r, 0.0)
    @inbounds for j in 1:size(A,2)
        @inbounds for i in 1:size(A,1)
            r[i]+=A[i,j]*v[j]
        end
    end
end

function infest(;df::DataFrame, P::Array{Float64}, P_log::Array{Float64}, T_matrix::Array{Float64}, area::Int, sites::Array{String}, tmax::Int, sim::Int, α::Float64)
    
    P_i=zeros(size(P,1))
    di=zeros(tmax, size(P,1))
    dit=zeros(Int, size(di, 2))
    i=zeros(size(di, 2))
    init_id=rand(1:1089)

    di[1, init_id]=1.0

    infest!(di=di, dit=dit, i=i, P_i=P_i, P_log=P_log, tmax=tmax)

    @views prop_inf=sum(di[tmax,:])/size(di,2)
    @views T_in=T_matrix[:,init_id]
    @views T_out=mean(T_matrix[:,init_id])
    @views P_in=P[:,init_id]
    @views P_out=mean(P[:,init_id])

    df[!,:area]=repeat([area], size(di, 2))
    df[!,:alpha]=repeat([α], size(di, 2))
    df[!,:site]=sites
    df[!,:inf_time]=dit
    df[!,:avg_inftime]=repeat([mean(dit)], size(di, 2))
    df[!,:prop_inf]=repeat([prop_inf], size(di, 2))
    df[!,:infest_focus]=repeat([init_id], size(di, 2))
    df[!,:T_in]=T_in
    df[!,:T_out]=repeat([T_out], size(di, 2))
    df[!,:P_in]=P_in
    df[!,:P_out]=repeat([P_out], size(di, 2))
    df[!,:nsim]=repeat([sim], size(di, 2))

    return df

end

function infest_sim(P::AbstractArray, Plog::AbstractArray, T::AbstractArray, s::Array{String}, tmax::Int, α::Float64, scenario::String, nsim::Int)

    for k in 1:length(P)

        results=[DataFrame([Int64,Float64,String,Int64,Float64,Float64,Int64,Float64,Float64,Float64,Float64,Int64],
        [:area, :alpha, :site, :inf_time, :avg_inftime, :prop_inf, :infest_focus, :T_in, :T_out, 
        :P_in, :P_out, :nsim], 1089) for _ in 1:nsim]

        Threads.@threads for n in 1:nsim

            infest(df=results[n], P=P[k], P_log=Plog[k], T_matrix=T[k], area=k, sites=s, tmax=tmax, α=α, sim=n)

        end

        CSV.write(datadir("$(scenario)_simulations_D$(k)_alpha_$α.csv"), vcat(results...))

    end

end