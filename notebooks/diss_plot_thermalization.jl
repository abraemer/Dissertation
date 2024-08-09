### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# ╔═╡ fe94fe12-3e9a-11ef-3a89-7f23e9876f8e
begin
	import Pkg
	Pkg.activate(".")
	using SpinModels, LinearAlgebra, CairoMakie,PlutoUI,SparseArrays,Statistics,SpinSymmetry, JLD2, Random, Optim, LogExpFunctions
	include("helpers.jl")
end

# ╔═╡ be9a4606-ac35-45c2-b86f-04a886e7c6ff
#using Random

# ╔═╡ a2341a96-f390-4045-a327-53c422074516
TableOfContents()

# ╔═╡ 6cded40a-265c-4079-b8e6-a06fd847339f
#N=17
N=13

# ╔═╡ 1f84a47e-2d37-4287-99d4-7b1de21ef25b
binomial(13,7)

# ╔═╡ 3fefa6cb-1416-44ff-95f0-61d4fda577ce
begin
	geom1 = Blockaded(SpinModels.NoisyChain(N,0.1); blockade=0.95)
	geom2 = Blockaded(SpinModels.Box(N,[N]); blockade=0.3)
	pos1 = positions(geom1; rng=Xoshiro(2))
	pos2 = positions(geom2; rng=Xoshiro(1))
	J1 = interaction_matrix(PowerLaw(3), geom1, pos1)
	J2 = interaction_matrix(PowerLaw(3), geom2, pos2)
end;

# ╔═╡ 9a654b24-2a37-47b7-af75-e44a1ab36c16
let (f,ax,_) = scatter(vec(pos1), 1.5ones(N); alpha=0.5, markersize=20)
	scatter!(ax, vec(pos2), ones(N); alpha=0.5, markersize=20)
	ylims!(ax, 0.5,2.0)
	xlims!(ax, -0.5, N+0.5)
	f
end

# ╔═╡ a486a8cf-190a-4ec7-8a53-cb736b6791f5
begin
	sector = transformationmatrix(symmetrized_basis(N, N÷2))
	ordered_H1 = J1*XXZ(0.4)
	disordered_H2 = J2*XXZ(0.4)
	eigenvals_ordered, eigenvects_ordered = eigen(Hermitian(Matrix(sector*sparse(ordered_H1)*sector')))
	eigenvals_ordered_normalized = (eigenvals_ordered .- eigenvals_ordered[1]) ./ (eigenvals_ordered[end]-eigenvals_ordered[1])
	eigenvals_disordered, eigenvects_disordered = eigen(Hermitian(Matrix(sector*sparse(disordered_H2)*sector')))
	eigenvals_disordered_normalized = (eigenvals_disordered .- eigenvals_disordered[1]) ./ (eigenvals_disordered[end]-eigenvals_disordered[1])
end;

# ╔═╡ 3ad05f78-e58f-460a-97a2-3aaaeb827468
sector

# ╔═╡ d4f6675a-298b-4df4-9f7b-20cc18a86677
md"""
# Energies
"""

# ╔═╡ a4a21eda-eab2-4697-9dc2-18318ac23892
let (f,ax,_) = hist(eigenvals_ordered_normalized)
	hist!(eigenvals_disordered_normalized)
	f
end

# ╔═╡ 361bf4c3-93b5-4d38-82a5-63348bb37cdc
let (f,ax,_) = scatter(eigenvals_ordered_normalized)
	scatter!(eigenvals_disordered_normalized)
	ax.ylabel="Normalized eigenstate energy"
	ax.xlabel="Eigenstate number"
	f
end

# ╔═╡ 623ebbaa-690c-4251-8397-51c7343d2268
md"""
# EEV
"""

# ╔═╡ 3321be84-53da-4d37-82c9-2ea1398b301b
# O = sparse(NN(Chain(N))*ZZ())
O = let couplings = zeros(N,N)
	couplings[3,4] = 1
	#sector*sparse(couplings*ZZ())*sector'
	#sector*sparse(ones(N,N)/binomial(N,2)*XX())*sector'
	#sector*sparse(ones(N)/N*X())*sector'
	sector*sparse([zeros(N÷2);0;0;1;zeros(N-3-N÷2)]*Z())*sector'
end

# ╔═╡ 7f0514fe-3567-4f3e-8342-646822d74fa4
begin
	EEV_ordered = [real(dot(ψ,O,ψ)) for ψ in eachcol(eigenvects_ordered)]
	EEV_disordered = [real(dot(ψ,O,ψ)) for ψ in eachcol(eigenvects_disordered)]
end

# ╔═╡ 8f3a4e8b-9d92-4101-9611-c76dfe72aa4a
let (f,ax,_) = scatter(EEV_ordered)
	scatter!(EEV_disordered)
	ax.xlabel="Eigenstate number"
	f
end

# ╔═╡ 9966a21f-090b-4d06-855b-f9bffa50bf6a
let (f,ax,_) = scatter(eigenvals_ordered_normalized, EEV_ordered)
	scatter!(eigenvals_disordered_normalized, EEV_disordered)
	ax.xlabel="Normalized eigenstate energy"
	f
end

# ╔═╡ 689ab82e-d3e3-407c-bb3e-a5512fed95e8
# windowed_evals, windowed_EEV = let window=10, r = 0:length(eigenvals)-window
# 	[mean(eigenvals[i+1:i+window]) for i in r],
# 	[mean(EEV[i+1:i+window]) for i in r]
# end

# ╔═╡ bc2e651a-9561-4770-af21-2df7cc624b29
# let (fig, ax,s) = scatter(eigenvals,EEV)
# 	lines!(ax, windowed_evals,windowed_EEV; color=:orange, label="moving average")
# 	ax.xlabel = L"Eigenstate energy $E_k$"
# 	ax.ylabel = L"Eigenstate expectation value $\langle k|O|k \rangle$"
# 	fig
# end

# ╔═╡ 35d9cf54-fabf-4810-85a7-1978eddfba38
md"""
# EON
"""

# ╔═╡ fdfdac0a-e91e-405b-aef5-55089c6b080f
# ψ0 = let v = zeros(2^N)
# 	v[1+sum(n->4^n, 0:N÷2-1)] = 1
# 	v
# end
# ψ0 = let v = zeros(2^N)
# 	v[1] = 1
# 	v
# end
# ψ0 = sector*normalize(ones(2^N))
# ψ0 = let up = [1,0], down = [0,1], state = [1]
# 	for i in 1:N
# 		state = kron(state,iseven(i+1) ? up : down)
# 	end
# 	normalize!(sector*state)
# end
ψ0 = let up = [1,0], down = [0,1], state = [1]
	for i in 1:N
		state = kron(state, 2i <= N+1 ? up : down)
	end
	normalize!(sector*state)
end
# ψ0 = let up = [1,1], down = [1,-1], state = [1]
# 	for i in 1:N
# 		state = kron(state,iseven(i+1) ? up : down)
# 	end
# 	normalize!(sector*state)
# end
# ψ0 = sector*[1;zeros(2^N-1)]
#ψ0 = normalize(sector*ones(2^N))
# ψ0 = sector*normalize(kron([[1,-im] for _ in 1:N]...))

# ╔═╡ ffb54136-d1d0-46fa-a126-9b239e08a49d
begin
	EON_ordered = [abs2(dot(ψ0,ψ)) for ψ in eachcol(eigenvects_ordered)]
	EON_disordered = [abs2(dot(ψ0,ψ)) for ψ in eachcol(eigenvects_disordered)]
end

# ╔═╡ bdd01911-03cb-488e-b180-52e282820901
begin
	E0_ordered = dot(EON_ordered, eigenvals_ordered)
	E0_disordered = dot(EON_disordered, eigenvals_disordered)
	E0_ordered_normalized = dot(EON_ordered, eigenvals_ordered_normalized)
	E0_disordered_normalized = dot(EON_disordered, eigenvals_disordered_normalized)

	varE_ordered = dot(EON_ordered, eigenvals_ordered .^2)-E0_ordered^2
	varE_disordered = dot(EON_disordered, eigenvals_disordered .^2)-E0_disordered
	varE_ordered_normalized = dot(EON_ordered, eigenvals_ordered_normalized .^2)-E0_ordered_normalized^2
	varE_disordered_normalized = dot(EON_disordered, eigenvals_disordered_normalized .^2) - E0_disordered_normalized^2
end;

# ╔═╡ de9d9482-9395-4426-83d6-240f35bfcb53
varE_ordered,varE_disordered,varE_ordered_normalized,varE_disordered_normalized

# ╔═╡ 6b1eb731-cf69-4e8d-96df-77a81b404896
let (f,ax,_) = scatter(EON_ordered)
	scatter!(EON_disordered)
	ax.xlabel="Eigenstate number"
	ax.ylabel=L"Eigenstate occupation number $|\langle k|\hat{O}|k\rangle|^2$"
	f
end

# ╔═╡ bb545fe2-0274-4000-8cba-9972e5fb3d19
let (f,ax,_) = scatter(eigenvals_ordered_normalized, EON_ordered)
	scatter!(eigenvals_disordered_normalized, EON_disordered)
	ax.xlabel="Normalized eigenstate energy"
	ax.ylabel=L"Eigenstate occupation number $|\langle k|\hat{O}|k\rangle|^2$"
	f
end

# ╔═╡ a3fb9d7f-d360-4bb8-8cc7-ee63b40911be
md"""
# Ensembles
"""

# ╔═╡ 0b76e947-96ca-4fbe-b062-ea84953552c4
sqrt(varE_ordered)

# ╔═╡ 0d7b60d4-5573-4c84-8f04-ab54ad620b7d
sqrt(varE_disordered)

# ╔═╡ 9f95d12f-0933-49d9-a031-1e720033697e
md"""
## Diagonal Ensemble
"""

# ╔═╡ f475a627-eb50-4218-9fae-2681417a6aba
begin
	diag_ensemble_ordered = dot(EEV_ordered, EON_ordered)
	diag_ensemble_disordered = dot(EEV_disordered, EON_disordered)
end

# ╔═╡ e432060e-1a3c-4804-bb00-18fcc4fb50a5
md"""
## Microcanonical ensemble
"""

# ╔═╡ 10138dfe-3c81-4273-9fc7-765c23a8afa2
begin
	microcanonical_ensemble_ordered = let maxind = argmax(EON_ordered),
			Δs = [5,10,20,50],
			indranges = [maxind-Δ:maxind+Δ for Δ in Δs]
		(; expvals=[mean(EEV_ordered[r]) for r in indranges], indranges)
	end
	microcanonical_ensemble_disordered = let maxind = argmax(EON_disordered),
			Δs = [5,10,20,50],
			indranges = [maxind-Δ:maxind+Δ for Δ in Δs]
		(; expvals=[mean(EEV_disordered[r]) for r in indranges], indranges)
	end
end

# ╔═╡ c8453242-276c-4816-8723-2b69424e6084
microcanonical_ensemble2_ordered = let EON = EON_ordered, eigenvals=eigenvals_ordered_normalized, EEV = EEV_ordered
let meanE = dot(EON,eigenvals),
		ΔEs = [0.005,0.01,0.02,0.05],
		indranges = [searchsortedfirst(eigenvals,meanE-ΔE):searchsortedlast(eigenvals,meanE+ΔE) for ΔE in ΔEs]
	(; expvals=[mean(EEV[r]) for r in indranges], indranges, meanE, ΔEs)
end
end

# ╔═╡ eea9917a-5d84-4c18-bb27-c33761f10cde
microcanonical_ensemble2_disordered = let EON = EON_disordered, eigenvals=eigenvals_disordered, EEV = EEV_disordered
let meanE = dot(EON,eigenvals),
		ΔEs = [0.005,0.01,0.02,0.05,0.1] .* Ref(-reduce(-,extrema(eigenvals))),
		indranges = [searchsortedfirst(eigenvals,meanE-ΔE):searchsortedlast(eigenvals,meanE+ΔE) for ΔE in ΔEs]
	(; expvals=[mean(EEV[r]) for r in indranges], indranges, meanE, ΔEs)
end
end

# ╔═╡ 0562a901-2f5c-4c12-b074-840e6298e6be
energy_bin_values = range(0,nextfloat(1.0);length=101)

# ╔═╡ 7424f8e0-b6db-4e99-b655-42fd17ef43c8
energy_bins_ordered = searchsortedfirst.(Ref(eigenvals_ordered_normalized), energy_bin_values)

# ╔═╡ 3465eeba-1b9a-4728-8a2f-98ec0acee4b4
binned_EON_ordered = [sum(EON_ordered[s:e-1]; init=0.0) for (s,e) in zip(energy_bins_ordered,energy_bins_ordered[2:end])]

# ╔═╡ e29f2f01-1bb7-4545-b2c1-708ec278187f
energy_bins_disordered = searchsortedfirst.(Ref(eigenvals_disordered_normalized), energy_bin_values)

# ╔═╡ ba2c6b54-b594-43eb-a29f-24d8cb8f4349
binned_EON_disordered = [sum(EON_disordered[s:e-1]; init=0.0) for (s,e) in zip(energy_bins_ordered,energy_bins_ordered[2:end])]

# ╔═╡ 492e15cc-e83e-46d1-a912-0a209d2837c8
let (f,a,_) = lines(energy_bin_values[1:end-1], binned_EON_ordered)
	lines!(a, energy_bin_values[1:end-1], binned_EON_disordered)
	f
end

# ╔═╡ cf8c4190-a439-416f-a1b0-94ec95831a2f
let (fig, ax, s) = scatter(eigenvals_ordered, EON_ordered)
	ax.xlabel = L"Eigenstate energy $E_k$"
	ax.ylabel = L"Eigenstate occupation $|\langle k| \psi \rangle|^2$"
	for r in microcanonical_ensemble2_ordered.indranges
		min,max = extrema(r)
		vlines!(ax, eigenvals_ordered[[min,max]]; linestyle=:dash)
	end
	vspan!(ax, E0_ordered - sqrt(varE_ordered), E0_ordered + sqrt(varE_ordered); color=(:blue,0.2))
	vlines!(ax, [E0_ordered])
	fig
end

# ╔═╡ a9828d2f-0624-4c20-a7c2-ba4ed507bf87
let E0 = dot(EON_ordered,eigenvals_ordered)
	(dot(EON_ordered, eigenvals_ordered .^2) - E0^2)/E0^2
end

# ╔═╡ 5ca34394-3a13-4cd5-965a-cd87dca95711
let E0 = dot(EON_disordered,eigenvals_disordered)
	(dot(EON_disordered, eigenvals_disordered .^2) - E0^2)/E0^2
end

# ╔═╡ e4f8b637-240f-4c2e-a891-9ed3ac54712f
let (fig, ax, s) = scatter(eigenvals_disordered_normalized, EON_disordered)
	ax.xlabel = L"Eigenstate energy $E_k$"
	ax.ylabel = L"Eigenstate occupation $|\langle k| \psi \rangle|^2$"
	# for r in microcanonical_ensemble2_disordered.indranges
	# 	min,max = extrema(r)
	# 	vlines!(ax, eigenvals_disordered_normalized[[min,max]]; linestyle=:dash)
	# end
	vspan!(ax, E0_disordered_normalized-sqrt(varE_disordered_normalized)/2, E0_disordered_normalized+sqrt(varE_disordered_normalized)/2; color=(:blue,0.2))
	vlines!(ax, [E0_disordered_normalized])
	scatter!(ax, eigenvals_disordered_normalized, 100softmax(-1.9*eigenvals_disordered_normalized))
	fig
end

# ╔═╡ 1951feae-88c3-41d4-978a-6b9903ea37f0
E0_disordered_normalized

# ╔═╡ 879d0c93-617e-4be7-9bae-e73cb5380b6c
md"""
## Canonical Ensemble
"""

# ╔═╡ 34471754-7722-4dc6-873b-a85893b05326
function energy_expval(eigenvals, β)
	return dot(eigenvals, softmax(-β .* eigenvals))
end

# ╔═╡ de8cd0d1-5889-4b70-a01f-e913a6c54a89
eigenvals_disordered_normalized[1]

# ╔═╡ b6bc6d8e-45d2-4ab6-87ee-60eaef5d42af
dot(eigenvals_disordered_normalized, EON_disordered)

# ╔═╡ c3761a65-e8dd-43c5-a02f-84ddb437a21d
energy_expval(eigenvals_disordered_normalized, -1.90535)

# ╔═╡ 7c94177c-796c-41e0-bf8e-8b8a6af23b51
scatter(eigenvals_disordered_normalized, softmax(-1.90535*eigenvals_disordered_normalized))

# ╔═╡ 691c3299-c738-4f40-89cc-fa23a095e8f1
dot(softmax(-1.90535*eigenvals_disordered_normalized), EEV_disordered)

# ╔═╡ 6044befd-5691-4487-b026-36ab74f77038
md"""
# Time evolution
"""

# ╔═╡ 23ae0570-4ac0-4a47-b75b-ca58f0a2ed93
function make_Ot(evals, U; ψ0=ψ0, O=O)
	Uψ0 = U'ψ0
	D = Diagonal(evals)
	function Ot(t)
		ψt = U*(cis(-D*t)*Uψ0)
		return real(dot(ψt, O, ψt))
	end
end

# ╔═╡ 988735ee-9c03-4015-8c7f-544c41de869e
tspan = range(0,1000;length=201)

# ╔═╡ 6073fa3d-004f-4b2d-a918-7ddfb47157fd
begin 
	Ot_ordered = make_Ot(eigenvals_ordered, eigenvects_ordered)
	Ot_disordered = make_Ot(eigenvals_disordered, eigenvects_disordered)
end

# ╔═╡ b2e0d795-60c9-4cdf-871b-557ff9c5412b
Ot_data_disordered = let eigenvals = eigenvals_disordered, eigenvects=eigenvects_disordered
	let Uψ0 = eigenvects'ψ0, D = Diagonal(eigenvals)
	[let ψt = eigenvects*(cis(-D*t)*Uψ0)
		real(dot(ψt, O, ψt))
	end for t in tspan]
end
end

# ╔═╡ f1eb54e4-7a6a-4a31-88c8-4070602c7884
median(sum(J2; dims=1)),median(sum(J1; dims=1))

# ╔═╡ 7f75e93e-7ec1-4cce-aef6-5f04ea4eeefd
let tspan = range(0,40;length=201),
	fig = Figure(), ax = Axis(fig[1,1]) 
	l1 = lines!(ax, tspan, Ot_ordered)
	l2 = lines!(ax, tspan, Ot_disordered)
	hlines!(ax, [diag_ensemble_ordered]; color=l1.color, linestyle=:dash, label="diagonal ensemble")
	hlines!(ax, [diag_ensemble_disordered]; color=l2.color, linestyle=:dash, label="diagonal ensemble")
	ax.xlabel="Time t [1/J]"
	ax.ylabel=L"\langle\! \hat{O} \rangle"
	hlines!(ax, [microcanonical_ensemble2_ordered.expvals[end]]; color=l1.color, linestyle=:dot, label="microcan. ensemble")
	hlines!(ax, [microcanonical_ensemble2_disordered.expvals[end]]; color=l2.color, linestyle=:dot, label="microcan. ensemble")
	# for e in microcanonical_ensemble_ordered.expvals
	# 	hlines!(ax, [e]; linestyle=:dash, label="diagonal ensemble")
	# end
	fig
end

# ╔═╡ cef8c186-f93e-4dca-8418-999327a25afd
let (fig, ax, l) = lines(tspan, Ot_data_disordered)
	hlines!(ax, [diag_ensemble_disordered]; color=:grey, linestyle=:dash, label="diagonal ensemble")
	for e in microcanonical_ensemble_disordered.expvals
		hlines!(ax, [e]; linestyle=:dash, label="diagonal ensemble")
	end
	fig
end

# ╔═╡ 565161a7-aa79-4f3a-947f-f23826da3850
md"""
# Plots
"""

# ╔═╡ 7a26d100-9b5e-4cd0-82b0-1978082805f7
with_theme(theme(; height=2, width=2)) do
	fig = Figure()

	ax0 = Axis(fig[0,1:2]; title="Positions")
	scatter!(ax0, vec(pos1), 1.5ones(N); alpha=0.9, markersize=15, label="weak disorder")
	scatter!(ax0, vec(pos2), ones(N); alpha=0.9, markersize=15, label="strong disorder")
	ylims!(ax0, 0.0,1.75)
	xlims!(ax0, 0, N+0.5)
	hidedecorations!(ax0)
	axislegend(ax0; position=:cb, orientation=:horizontal, framevisible=false)
	Label(fig[0,1,TopLeft()],"(a)"; tellwidth=false, tellheight=false)
	
	ax1 = Axis(fig[1,1]; 
		ylabel = L"$\langle k|O|k\rangle$", 
		xlabel = L"Energy (normalized) $E_k$")
	scatter!(ax1, eigenvals_ordered_normalized,
		EEV_ordered)
	scatter!(ax1, eigenvals_disordered_normalized,
		EEV_disordered)
	Label(fig[1,1,TopLeft()],"(b)"; tellwidth=false, tellheight=false)
	
	ax2 = Axis(fig[1,2]; 
		ylabel = L"|c_k|^2",
		xlabel = L"Energy (normalized) $E_k$", )
	vspan!(ax2, E0_ordered_normalized-sqrt(varE_ordered_normalized), E0_ordered_normalized + sqrt(varE_ordered_normalized); color=(Makie.wong_colors()[1],0.2))
	vspan!(ax2, E0_disordered_normalized-sqrt(varE_disordered_normalized), E0_disordered_normalized + sqrt(varE_disordered_normalized); color=(Makie.wong_colors()[2],0.2))
	# scatter!(ax2, eigenvals_ordered_normalized, EON_ordered; label="EON")
	# scatter!(ax2, eigenvals_disordered_normalized, EON_disordered; label="EON")
	scatter!(ax2, energy_bin_values[1:end-1], binned_EON_ordered)
	scatter!(ax2, energy_bin_values[1:end-1], binned_EON_disordered)
	# ax2.xlabel = L"Eigenstate energy $E_k$"
	# ax2.ylabel = L"Eigenstate occupation $|\langle k| \psi \rangle|^2$"
	# vlines!(ax2, [E0_ordered_normalized])
	# vlines!(ax2, [E0_disordered_normalized])
	Label(fig[1,2,TopLeft()],"(c)"; tellwidth=false, tellheight=false)

	tspan = range(0,40;length=201)
	ax3 = Axis(fig[2,1:2])
	l1 = lines!(ax3, tspan, Ot_ordered; label="ordered")
	l2 = lines!(ax3, tspan, Ot_disordered; label="disordered")
	hlines!(ax3, [diag_ensemble_ordered]; color=l1.color, linestyle=:dash)
	hlines!(ax3, [diag_ensemble_disordered]; color=l2.color, linestyle=:dash)
	ax3.xlabel="Time t [1/J]"
	ax3.ylabel=L"\langle\! \hat{O}(t) \rangle"
	hlines!(ax3, [microcanonical_ensemble2_ordered.expvals[end]]; color=l1.color, linestyle=:dot)
	hlines!(ax3, [microcanonical_ensemble2_disordered.expvals[end]]; color=l2.color, linestyle=:dot)
	Label(fig[2,1,TopLeft()],"(d)"; tellwidth=false, tellheight=false)
	
	rowsize!(fig.layout, 0, Relative(0.3))
	fig
end

# ╔═╡ b9acd0d1-1cb5-4bf4-ba9d-c8035c5af28e
E0_ordered_normalized

# ╔═╡ 64dcc216-57dc-4e40-b71b-371884d01e30
md"""
# Unrelated stuff
"""

# ╔═╡ ef01e908-8b5c-419f-99cc-a28791d1aa3a
let r = range(0,2.01π/2;length=200), b = 0.5π/2
	f = Figure()
	a = Axis(f[1,1])
	lines!(a, r, x->-b^2/2/tan(x); label="center approx")
	vlines!(a,[b,π-b]; color=:grey, linestyle=:dot)
	lines!(a, filter(<(b),r), a->a-b-a^2/2/tan(b); label="left approx")
	lines!(a, filter(>(π-b),r), a->(-pi+a)+b+(pi-a)^2/2/tan(b); label="right approx")
	lines!(a, r, a->(1-cos(b))*(a-π/2); label="linear approx (center)")
	lines!(a, r, a->2b/π*(a-π/2); label="linear approx (whole)")
	lines!(a, r, x->acos(cos(x))-acos(cos(x)*cos(b));linestyle=:dash, color=:black)
	axislegend(a; position=:lt)
	#lines!(a, filter(<(b),r), a->pi/2-ab-a^2/2/tan(b))
	f
end

# ╔═╡ 38bab1fd-4246-49c5-8f91-95045f4162c8
theme_latexfonts()

# ╔═╡ 354dac55-fef0-4f49-931d-5ad97ba7194f
Makie.to_font("TeX Gyre Pagella Italic")

# ╔═╡ 3a248255-dcad-410c-9f1b-f5a1fef3af1f
Makie.to_font("Math Pazo")

# ╔═╡ fa2dc924-a72f-46ea-9979-7ffe55e23ed3
let r = range(0,1.95π/2;length=200), b = 0.1π/2
	f = Figure()
	a = Axis(f[1,1])
	lines!(a, r, x->x+cot(x)*b^2/2; label="center approx")
	vlines!(a,[b,π-b]; color=:grey, linestyle=:dot)
	lines!(a, filter(<(b),r), a->b+a^2/2/tan(b); label="left approx")
	lines!(a, filter(>(π-b),r), a->(pi-b)-(pi-a)^2/2/tan(b); label="right approx")
	lines!(a, r, a->π/2+(1-b^2/2)*(a-π/2); label="linear approx (center)")
	lines!(a, r, a->π/2+(π-2b)/π*(a-π/2); label="linear approx (whole)")
	lines!(a, r, x->acos(cos(x)*cos(b));linestyle=:dash, color=:black)
	axislegend(a; position=:lt)
	#lines!(a, filter(<(b),r), a->pi/2-ab-a^2/2/tan(b))
	f
end

# ╔═╡ 707d739d-e3e7-4db4-881a-0253bb1210bb
cot(pi/2)

# ╔═╡ 00e7e06d-634f-4bc3-adef-e4aa8f06d35c
let r = range(0,π;length=200), b = 0.2π/2
	f,a,_ = lines(r, x->x-acos(cos(x)*cos(b)))
	#lines!(a, r, x->-b^2/2/tan(x))
	f
end

# ╔═╡ 25c40d98-d67b-4bc6-994b-c343f3fb06fe
md"""
# Single particle localization
"""

# ╔═╡ a9cbbb3c-ec38-4742-8595-e4953c0a2c97
Hs2 = let N=100, W=3, Δ=1,
	fields = 2W*(rand(Xoshiro(2),N) .- 0.5)/2, # uniform[-W,W] * S_z
	total = 0sum(fields), # actually does not matter due to U(1)
	potentials = -2fields #.+ total .+ Δ/2*(N-5) # same here
	potentials[1] += Δ/4
	potentials[end] += Δ/4
	LinearAlgebra.SymTridiagonal(potentials, ones(N-1)/2)
end

# ╔═╡ bec82953-b272-4c90-96be-a4147e211c28
eigen(Hs2)

# ╔═╡ 0c362908-f122-4c9b-8ba1-909262cd80b4
with_theme(theme()) do
	vecs = abs2.(eigen(Hs2).vectors)
	fig = Figure()
	ax = Axis(fig[1,1]; yscale=log10, xlabel=L"site index $k$", ylabel=L"Overlap $|\langle \downarrow\ldots\downarrow|S_-^{(k)}|\psi\rangle|^2$")
	for v in eachcol(vecs)[1:6]
		lines!(ax, v)
	end
	ax.xticks=0:25:100
	ylims!(ax, 1e-20, 1e1)
	fig
end |> save_and_display("single-spin-localization-xxz", "part1")

# ╔═╡ Cell order:
# ╠═fe94fe12-3e9a-11ef-3a89-7f23e9876f8e
# ╠═be9a4606-ac35-45c2-b86f-04a886e7c6ff
# ╠═a2341a96-f390-4045-a327-53c422074516
# ╠═6cded40a-265c-4079-b8e6-a06fd847339f
# ╠═1f84a47e-2d37-4287-99d4-7b1de21ef25b
# ╠═3ad05f78-e58f-460a-97a2-3aaaeb827468
# ╠═3fefa6cb-1416-44ff-95f0-61d4fda577ce
# ╠═9a654b24-2a37-47b7-af75-e44a1ab36c16
# ╠═a486a8cf-190a-4ec7-8a53-cb736b6791f5
# ╟─d4f6675a-298b-4df4-9f7b-20cc18a86677
# ╠═a4a21eda-eab2-4697-9dc2-18318ac23892
# ╠═361bf4c3-93b5-4d38-82a5-63348bb37cdc
# ╟─623ebbaa-690c-4251-8397-51c7343d2268
# ╠═3321be84-53da-4d37-82c9-2ea1398b301b
# ╠═7f0514fe-3567-4f3e-8342-646822d74fa4
# ╠═8f3a4e8b-9d92-4101-9611-c76dfe72aa4a
# ╠═9966a21f-090b-4d06-855b-f9bffa50bf6a
# ╠═689ab82e-d3e3-407c-bb3e-a5512fed95e8
# ╠═bc2e651a-9561-4770-af21-2df7cc624b29
# ╠═35d9cf54-fabf-4810-85a7-1978eddfba38
# ╠═fdfdac0a-e91e-405b-aef5-55089c6b080f
# ╠═ffb54136-d1d0-46fa-a126-9b239e08a49d
# ╠═bdd01911-03cb-488e-b180-52e282820901
# ╠═de9d9482-9395-4426-83d6-240f35bfcb53
# ╠═6b1eb731-cf69-4e8d-96df-77a81b404896
# ╠═bb545fe2-0274-4000-8cba-9972e5fb3d19
# ╟─a3fb9d7f-d360-4bb8-8cc7-ee63b40911be
# ╠═0b76e947-96ca-4fbe-b062-ea84953552c4
# ╠═0d7b60d4-5573-4c84-8f04-ab54ad620b7d
# ╟─9f95d12f-0933-49d9-a031-1e720033697e
# ╠═f475a627-eb50-4218-9fae-2681417a6aba
# ╟─e432060e-1a3c-4804-bb00-18fcc4fb50a5
# ╠═10138dfe-3c81-4273-9fc7-765c23a8afa2
# ╠═c8453242-276c-4816-8723-2b69424e6084
# ╠═eea9917a-5d84-4c18-bb27-c33761f10cde
# ╠═0562a901-2f5c-4c12-b074-840e6298e6be
# ╠═7424f8e0-b6db-4e99-b655-42fd17ef43c8
# ╠═3465eeba-1b9a-4728-8a2f-98ec0acee4b4
# ╠═e29f2f01-1bb7-4545-b2c1-708ec278187f
# ╠═ba2c6b54-b594-43eb-a29f-24d8cb8f4349
# ╠═492e15cc-e83e-46d1-a912-0a209d2837c8
# ╠═cf8c4190-a439-416f-a1b0-94ec95831a2f
# ╠═a9828d2f-0624-4c20-a7c2-ba4ed507bf87
# ╠═5ca34394-3a13-4cd5-965a-cd87dca95711
# ╠═e4f8b637-240f-4c2e-a891-9ed3ac54712f
# ╠═1951feae-88c3-41d4-978a-6b9903ea37f0
# ╟─879d0c93-617e-4be7-9bae-e73cb5380b6c
# ╠═34471754-7722-4dc6-873b-a85893b05326
# ╠═de8cd0d1-5889-4b70-a01f-e913a6c54a89
# ╠═b6bc6d8e-45d2-4ab6-87ee-60eaef5d42af
# ╠═c3761a65-e8dd-43c5-a02f-84ddb437a21d
# ╠═7c94177c-796c-41e0-bf8e-8b8a6af23b51
# ╠═691c3299-c738-4f40-89cc-fa23a095e8f1
# ╠═6044befd-5691-4487-b026-36ab74f77038
# ╠═23ae0570-4ac0-4a47-b75b-ca58f0a2ed93
# ╠═988735ee-9c03-4015-8c7f-544c41de869e
# ╠═6073fa3d-004f-4b2d-a918-7ddfb47157fd
# ╠═b2e0d795-60c9-4cdf-871b-557ff9c5412b
# ╠═f1eb54e4-7a6a-4a31-88c8-4070602c7884
# ╠═7f75e93e-7ec1-4cce-aef6-5f04ea4eeefd
# ╠═cef8c186-f93e-4dca-8418-999327a25afd
# ╟─565161a7-aa79-4f3a-947f-f23826da3850
# ╠═7a26d100-9b5e-4cd0-82b0-1978082805f7
# ╠═b9acd0d1-1cb5-4bf4-ba9d-c8035c5af28e
# ╟─64dcc216-57dc-4e40-b71b-371884d01e30
# ╠═ef01e908-8b5c-419f-99cc-a28791d1aa3a
# ╠═38bab1fd-4246-49c5-8f91-95045f4162c8
# ╠═354dac55-fef0-4f49-931d-5ad97ba7194f
# ╠═3a248255-dcad-410c-9f1b-f5a1fef3af1f
# ╠═fa2dc924-a72f-46ea-9979-7ffe55e23ed3
# ╠═707d739d-e3e7-4db4-881a-0253bb1210bb
# ╠═00e7e06d-634f-4bc3-adef-e4aa8f06d35c
# ╟─25c40d98-d67b-4bc6-994b-c343f3fb06fe
# ╠═a9cbbb3c-ec38-4742-8595-e4953c0a2c97
# ╠═bec82953-b272-4c90-96be-a4147e211c28
# ╠═0c362908-f122-4c9b-8ba1-909262cd80b4
