"Make a mmaped array of type A"
function typemmap(
    ::Type{A},
    dims::NTuple{N,Int};
    basedir::String = tempdir(),
    suffix::String = "",
    fpath::String = joinpath(basedir, basename(tempname()) * suffix),
    autoclean::Bool = true,
    readonly::Bool = false,
) where {A<:AbstractArray,N}
    arr = open(fpath, ifelse(readonly, "r+", "w+")) do io
        Mmap.mmap(io, A, dims; grow = true)
    end
    if autoclean
        finalizer(arr) do a
            try
                rm(fpath; force = true)
            catch
            end
        end
    end
    arr, fpath
end

function typemmap(a::AbstractArray{T,N}, args...; kwargs...) where {T,N}
    typemmap(Array{T,N}, size(a), args...; kwargs...)
end

function to_mmap(a::AbstractArray, arrtype::DataType = typeof(a); kwargs...)
    mma, path = typemmap(arrtype, size(a); kwargs...)
    copyto!(mma, a)
    mma, path
end

function only_matches(reg::Regex, strs::AbstractArray{<:AbstractString})
    n_s = length(strs)
    matches = Vector{RegexMatch}(undef, n_s)
    out_no = 0
    for str in strs
        maybe_match = match(reg, str)
        if !isnothing(maybe_match)
            out_no += 1
            matches[out_no] = maybe_match
        end
    end
    resize!(matches, out_no)
    matches
end
