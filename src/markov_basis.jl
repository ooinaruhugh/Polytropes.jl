using Oscar
import Oscar: IntegerUnion

function markov_basis(S:: Union{MatElem{U}, <: AbstractMatrix{U}};
                      use_kernel::Bool=true) where U <: IntegerUnion

    # S = ZZRingElem_mat(S)
    if !use_kernel
        S = transpose(S)
    end

    pm_M = Polymake.fulton.markov_basis(
        S,
        use_kernel=use_kernel
    )

    return matrix(ZZ, pm_M)*(-1)
end
