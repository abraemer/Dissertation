### A Pluto.jl notebook ###
# v0.19.43

using Markdown
using InteractiveUtils

# ╔═╡ fe94fe12-3e9a-11ef-3a89-7f23e9876f8e
begin
	import Pkg
	Pkg.activate(".")
	using SpinModels, LinearAlgebra, CairoMakie,PlutoUI,SparseArrays,Statistics,SpinSymmetry, JLD2
	include("helpers.jl")
end

# ╔═╡ a2341a96-f390-4045-a327-53c422074516
TableOfContents()

# ╔═╡ 6cded40a-265c-4079-b8e6-a06fd847339f
N=12

# ╔═╡ a486a8cf-190a-4ec7-8a53-cb736b6791f5
# generic_H = PowerLaw(3)(
# 	NoisyChain(N,0.1))*XXZ(-0.4) +
# 	0.8*Y(ones(N)+0.2rand(N)) -
# 	0.3*X(ones(N)+0.2rand(N))
generic_H = NN(Chain(N))*XXZ(-0.4) +
	0.8*Y(ones(N)+0.2rand(N)) -
	0.3*X(ones(N)+0.2rand(N))

# ╔═╡ 563bf443-fd28-4d61-ace6-062488fba075
eigenvals, eigenvects = eigen(Hermitian(Matrix(generic_H)))

# ╔═╡ d4f6675a-298b-4df4-9f7b-20cc18a86677
md"""
# Energies
"""

# ╔═╡ a4a21eda-eab2-4697-9dc2-18318ac23892
hist(eigenvals)

# ╔═╡ 623ebbaa-690c-4251-8397-51c7343d2268
md"""
# EEV
"""

# ╔═╡ 3321be84-53da-4d37-82c9-2ea1398b301b
# O = sparse(NN(Chain(N))*ZZ())
O = let couplings = zeros(N,N)
	couplings[3,4] = 1
	sparse(couplings*ZZ())
end

# ╔═╡ 7f0514fe-3567-4f3e-8342-646822d74fa4
EEV = [real(dot(ψ,O,ψ)) for ψ in eachcol(eigenvects)]

# ╔═╡ 8f3a4e8b-9d92-4101-9611-c76dfe72aa4a
scatter(EEV)

# ╔═╡ 689ab82e-d3e3-407c-bb3e-a5512fed95e8
windowed_evals, windowed_EEV = let window=10, r = 0:length(eigenvals)-window
	[mean(eigenvals[i+1:i+window]) for i in r],
	[mean(EEV[i+1:i+window]) for i in r]
end

# ╔═╡ bc2e651a-9561-4770-af21-2df7cc624b29
let (fig, ax,s) = scatter(eigenvals,EEV)
	lines!(ax, windowed_evals,windowed_EEV; color=:orange, label="moving average")
	ax.xlabel = L"Eigenstate energy $E_k$"
	ax.ylabel = L"Eigenstate expectation value $\langle k|O|k \rangle$"
	fig
end

# ╔═╡ 35d9cf54-fabf-4810-85a7-1978eddfba38
md"""
# EON
"""

# ╔═╡ fdfdac0a-e91e-405b-aef5-55089c6b080f
ψ0 = [1;zeros(2^N-1)]
# ψ0 = let v = zeros(2^N)
# 	v[1+sum(n->4^n, 0:N÷2-1)] = 1
# 	v
# end
# ψ0 = let v = zeros(2^N)
# 	v[1] = 1
# 	v
# end
#ψ0 = normalize(ones(2^N))
#ψ0 = normalize(kron([[1,-im] for _ in 1:N]...))

# ╔═╡ ffb54136-d1d0-46fa-a126-9b239e08a49d
EON = [abs2(dot(ψ0,ψ)) for ψ in eachcol(eigenvects)]

# ╔═╡ bb545fe2-0274-4000-8cba-9972e5fb3d19
scatter(EON)

# ╔═╡ 6044befd-5691-4487-b026-36ab74f77038
md"""
# Time evolution
"""

# ╔═╡ 988735ee-9c03-4015-8c7f-544c41de869e
tspan = range(0,50;length=201)

# ╔═╡ 74d1bb38-af93-424e-93ec-1569f143b86c
Ot_data = let Uψ0 = eigenvects'ψ0, D = Diagonal(eigenvals)
	[let ψt = eigenvects*(cis(-D*t)*Uψ0)
		real(dot(ψt, O, ψt))
	end for t in tspan]
end

# ╔═╡ f475a627-eb50-4218-9fae-2681417a6aba
diag_ensemble = dot(EEV, EON)

# ╔═╡ 0a06af21-67a0-4e2c-b272-406c6770598b
dot(EON,eigenvals)

# ╔═╡ 62ffbd0c-283c-45ca-8345-d898e047c8ca
searchsortedfirst(eigenvals, dot(EON,eigenvals))

# ╔═╡ 95110bf2-50d1-4e70-a32f-d5a22400b678
eigenvals[709:711]

# ╔═╡ a72e13f3-55a4-4590-b270-3b739ea4bf9d
argmax(EON)

# ╔═╡ 10138dfe-3c81-4273-9fc7-765c23a8afa2
microcanonical_ensemble = let maxind = argmax(EON),
		Δs = [5,10,20,50],
		indranges = [maxind-Δ:maxind+Δ for Δ in Δs]
	(; expvals=[mean(EEV[r]) for r in indranges], indranges)
end

# ╔═╡ cf8c4190-a439-416f-a1b0-94ec95831a2f
let (fig, ax, s) = scatter(eigenvals, EON)
	ax.xlabel = L"Eigenstate energy $E_k$"
	ax.ylabel = L"Eigenstate occupation $|\langle k| \psi \rangle|^2$"
	for r in microcanonical_ensemble.indranges
		min,max = extrema(r)
		vlines!(ax, eigenvals[[min,max]]; linestyle=:dash)
	end
	vlines!(ax, [-10.325])
	fig
end

# ╔═╡ 7f75e93e-7ec1-4cce-aef6-5f04ea4eeefd
let (fig, ax, l) = lines(tspan, Ot_data)
	hlines!(ax, [diag_ensemble]; color=:grey, linestyle=:dash, label="diagonal ensemble")
	for e in microcanonical_ensemble.expvals
		hlines!(ax, [e]; linestyle=:dash, label="diagonal ensemble")
	end
	fig
end

# ╔═╡ c8453242-276c-4816-8723-2b69424e6084
microcanonical_ensemble2 = let meanE = dot(EON,eigenvals),
		ΔEs = [0.005,0.01,0.02,0.05] .* Ref(-reduce(-,extrema(eigenvals))),
		indranges = [searchsortedfirst(eigenvals,meanE-ΔE):searchsortedlast(eigenvals,meanE+ΔE) for ΔE in ΔEs]
	(; expvals=[mean(EEV[r]) for r in indranges], indranges, meanE, ΔEs)
end

# ╔═╡ ec222f8e-1d25-447b-9f26-4a714a466246
argmax(EON)

# ╔═╡ Cell order:
# ╠═fe94fe12-3e9a-11ef-3a89-7f23e9876f8e
# ╠═a2341a96-f390-4045-a327-53c422074516
# ╠═6cded40a-265c-4079-b8e6-a06fd847339f
# ╠═a486a8cf-190a-4ec7-8a53-cb736b6791f5
# ╠═563bf443-fd28-4d61-ace6-062488fba075
# ╟─d4f6675a-298b-4df4-9f7b-20cc18a86677
# ╠═a4a21eda-eab2-4697-9dc2-18318ac23892
# ╟─623ebbaa-690c-4251-8397-51c7343d2268
# ╠═3321be84-53da-4d37-82c9-2ea1398b301b
# ╠═7f0514fe-3567-4f3e-8342-646822d74fa4
# ╠═8f3a4e8b-9d92-4101-9611-c76dfe72aa4a
# ╠═689ab82e-d3e3-407c-bb3e-a5512fed95e8
# ╠═bc2e651a-9561-4770-af21-2df7cc624b29
# ╠═35d9cf54-fabf-4810-85a7-1978eddfba38
# ╠═fdfdac0a-e91e-405b-aef5-55089c6b080f
# ╠═ffb54136-d1d0-46fa-a126-9b239e08a49d
# ╠═bb545fe2-0274-4000-8cba-9972e5fb3d19
# ╠═cf8c4190-a439-416f-a1b0-94ec95831a2f
# ╠═6044befd-5691-4487-b026-36ab74f77038
# ╠═988735ee-9c03-4015-8c7f-544c41de869e
# ╠═74d1bb38-af93-424e-93ec-1569f143b86c
# ╠═f475a627-eb50-4218-9fae-2681417a6aba
# ╠═7f75e93e-7ec1-4cce-aef6-5f04ea4eeefd
# ╠═0a06af21-67a0-4e2c-b272-406c6770598b
# ╠═62ffbd0c-283c-45ca-8345-d898e047c8ca
# ╠═95110bf2-50d1-4e70-a32f-d5a22400b678
# ╠═a72e13f3-55a4-4590-b270-3b739ea4bf9d
# ╠═10138dfe-3c81-4273-9fc7-765c23a8afa2
# ╠═c8453242-276c-4816-8723-2b69424e6084
# ╠═ec222f8e-1d25-447b-9f26-4a714a466246
