function prob_matrix(D::AbstractArray, scaling, α)

    P=zeros(size(D,1), size(D,2))
    A=zeros(size(D,1), size(D,2))

    @. A[D>0.0]=1.0
    ∑aij=dropdims(sum(A, dims=2), dims=2)
    @. A=A/∑aij
    @. A[isnan(A)]=0.0

    @. P=exp(-α*(D/scaling)^2)
    @. P=P*A
    
    return P

end

function T_prob(P)

    T=inv(I-P)

    return T

end

function infection!(;di::Array{Float64}, dit::Array{Int}, P_i::Array{Float64}, time::Int)
    @inbounds for k in 1:length(P_i)
        P_i[k]=(1.0 - (exp(P_i[k])))*(1.0-di[time,k])
        r=rand()
        if P_i[k]>=r
            di[time+1,k]=1.0
            dit[k]=time
        elseif P_i[k]<=r
            di[time+1,k]=di[time,k]
        end
    end
end

function infest!(;di::Array{Float64}, dit::Array{Int}, i::Array{Float64}, P_i::Array{Float64}, P_log::Array{Float64}, tmax::Int)

    @inbounds for t in 1:(tmax-1)

        i_vector!(i=i, di=di, time=t)
        mul!(P_i, P_log, i)
        infection!(di=di, dit=dit, P_i=P_i, time=t)

    end
end

function P_log(P)
    return @. log(1.0-P)
end

function i_vector!(;i::Array{Float64}, di::Array{Float64}, time::Int64)
    @inbounds for j in 1:size(di,2)
        i[j]=di[time,j]
    end
end