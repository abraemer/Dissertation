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

# ╔═╡ 914f946b-2a6a-4006-9928-d427638a8c75
md"""
# ETH
"""

# ╔═╡ 6cded40a-265c-4079-b8e6-a06fd847339f
#N=17
N=13

# ╔═╡ 95f064ef-e3d5-4bf1-8200-8c8bfeec3516
binomial(13,7)

# ╔═╡ 7e8f7f04-7db7-41bd-b256-5f54d1218ff5
md"""
## Random positions
"""

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

# ╔═╡ d4f6675a-298b-4df4-9f7b-20cc18a86677
md"""
## Energies
"""

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
## EEV
"""

# ╔═╡ 3321be84-53da-4d37-82c9-2ea1398b301b
# O = sparse(NN(Chain(N))*ZZ())
O = let couplings = zeros(N,N)
	couplings[5,4] = 1
	#sector*sparse(couplings*ZZ())*sector'
	#sector*sparse(ones(N,N)/binomial(N,2)*XX())*sector'
	#sector*sparse(ones(N)/N*X())*sector'
	sector*sparse([zeros(N÷2);0;0;-0.5;zeros(N-3-N÷2)]*Z())*sector'
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

# ╔═╡ 35d9cf54-fabf-4810-85a7-1978eddfba38
md"""
## EON
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
## Ensembles
"""

# ╔═╡ 9f95d12f-0933-49d9-a031-1e720033697e
md"""
### Diagonal Ensemble
"""

# ╔═╡ f475a627-eb50-4218-9fae-2681417a6aba
begin
	diag_ensemble_ordered = dot(EEV_ordered, EON_ordered)
	diag_ensemble_disordered = dot(EEV_disordered, EON_disordered)
end

# ╔═╡ e432060e-1a3c-4804-bb00-18fcc4fb50a5
md"""
### Microcanonical ensemble
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
binned_EON_disordered = [sum(EON_disordered[s:e-1]; init=0.0) for (s,e) in zip(energy_bins_disordered,energy_bins_disordered[2:end])]

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
### Canonical Ensemble
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
## Time evolution
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
## Plot
"""

# ╔═╡ 7a26d100-9b5e-4cd0-82b0-1978082805f7
with_theme(theme(; height=2.5, width=3)) do
	fig = Figure()

	ax0 = Axis(fig[0,1:2]; title="Positions")
	scatter!(ax0, vec(pos1), 1.5ones(N); alpha=0.9, markersize=15, label="weak disorder")
	scatter!(ax0, vec(pos2), ones(N); alpha=0.9, markersize=15, label="strong disorder")
	ylims!(ax0, 0.25,1.75)
	xlims!(ax0, 0, N+0.5)
	hidedecorations!(ax0)
	axislegend(ax0; position=:cb, orientation=:horizontal, framevisible=false, margin=(0,0,0,0), padding=0)
	Label(fig[0,1,TopLeft()],"(a)"; tellwidth=false, tellheight=false)
	
	ax1 = Axis(fig[1,1]; 
		ylabel = L"$\langle k|O|k\rangle$", 
		xlabel = L"Energy (normalized) $E_k$",
		title = "Eigenstate expectation")
	scatter!(ax1, eigenvals_ordered_normalized,
		EEV_ordered; markersize=5)
	scatter!(ax1, eigenvals_disordered_normalized,
		EEV_disordered; markersize=5)
	Label(fig[1,1,TopLeft()],"(b)"; tellwidth=false, tellheight=false)
	ax1.yticks = -0.5:0.5:0.5
	
	ax2 = Axis(fig[1,2];
		ylabel = L"|c_k|^2",
		xlabel = L"Energy (normalized) $E_k$",
		title = "Eigenstate occupation")
	vspan!(ax2, E0_ordered_normalized-sqrt(varE_ordered_normalized), E0_ordered_normalized + sqrt(varE_ordered_normalized); color=(Makie.wong_colors()[1],0.2))
	vspan!(ax2, E0_disordered_normalized-sqrt(varE_disordered_normalized), E0_disordered_normalized + sqrt(varE_disordered_normalized); color=(Makie.wong_colors()[2],0.2))
	# scatter!(ax2, eigenvals_ordered_normalized, EON_ordered; label="EON")
	# scatter!(ax2, eigenvals_disordered_normalized, EON_disordered; label="EON")
	scatter!(ax2, energy_bin_values[2:end], binned_EON_ordered; markersize=5)
	scatter!(ax2, energy_bin_values[2:end], binned_EON_disordered; markersize=5)
	# ax2.xlabel = L"Eigenstate energy $E_k$"
	# ax2.ylabel = L"Eigenstate occupation $|\langle k| \psi \rangle|^2$"
	# vlines!(ax2, [E0_ordered_normalized])
	# vlines!(ax2, [E0_disordered_normalized])
	Label(fig[1,2,TopLeft()],"(c)"; tellwidth=false, tellheight=false)

	tspan = range(0,40;length=201)
	ax3 = Axis(fig[2,1:2]; 
		xlabel=L"Time t [$1/J$]",
		ylabel=L"\langle\! \hat{O}(t) \rangle",
		title="Time trace")
	l1 = lines!(ax3, tspan, Ot_ordered; label="ordered")
	l2 = lines!(ax3, tspan, Ot_disordered; label="disordered")
	hlines!(ax3, [diag_ensemble_ordered]; color=l1.color, linestyle=:dash)
	hlines!(ax3, [diag_ensemble_disordered]; color=l2.color, linestyle=:dash)
	hlines!(ax3, [microcanonical_ensemble2_ordered.expvals[end]]; color=l1.color, linestyle=:dot)
	hlines!(ax3, [microcanonical_ensemble2_disordered.expvals[end]]; color=l2.color, linestyle=:dot)
	xlims!(ax3, -1, 48)
	ax3.yticks = -0.5:0.25:0.5
	text!(ax3, 40, diag_ensemble_ordered+0.02; color=l1.color, text=L"\langle\! O\rangle_{diag} = \langle\! O\rangle_{mc}", align=(:left,:bottom))

	text!(ax3, 44, microcanonical_ensemble2_disordered.expvals[end]-0.02; color=l2.color, text=L"\langle\! O\rangle_{mc}", align=(:center,:top))
	text!(ax3, 44, diag_ensemble_disordered-0.02; color=l2.color, text=L"\langle\! O\rangle_{diag}", align=(:center,:top))
	
	Label(fig[2,1,TopLeft()],"(d)"; tellwidth=false, tellheight=false)
	
	rowsize!(fig.layout, 0, Relative(0.2))
	fig
end |> save_and_display("ETH-thermalization", "part1")

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

# ╔═╡ f9119509-6e24-48ee-a1fb-71b970afc11b
md"""
# Pair timecrystal
"""

# ╔═╡ d1fd0411-db74-450f-bc08-21719117f99a
function random_couplings(N; β=1)
	# β = d/α <= 1
	return (1 ./ rand(Xoshiro(2), N) .- 1) .^ (1/β)
end

# ╔═╡ c6c227e9-0958-4844-8295-9df5055c8cb7
function single_pair_analytic(;t_wait, ϕ, J=1)
	a = J*t_wait/4
	b = ϕ
	c = acos(cos(a)*cos(b))
	f = sin(a)*cos(b)/sin(c)
	n -> 0.5*(cos(a*n)*cos(c*n) + sin(a*n)*sin(c*n)*f)
	#n -> 2/π*sin(n*(1-cos(ϕ))*π/2)/(n*(1-cos(ϕ)))
end

# ╔═╡ b7c8eb86-1182-4325-8930-14ef52cca045
function single_pair_analytic_approx(;t_wait, ϕ, J=1)
	# valid for ϕ<J*twait
	a = pi/2-abs(pi/2-mod(J*t_wait/4,2pi/2))
	b = ϕ
	#c = acos(cos(a)*cos(b))
	
	#f=1
	#f = (a < b) ? sqrt(a^2/(a^2+b^2)) : 1-b^2/2/sin(a)^2
	#f = sin(a)*cos(b)/sin(c)
	
	# Δac = b^2/2*(a-π/2) # not precise enough
	# Δac = a-c
	#(a<b) && return n->0.5*cos(n*b)*cos(n*a) # more precise: b+a^2/2*cot(b)
	n -> 0.5*cos(n*b^2/2*cot(a))
end

# ╔═╡ c2a2d9f4-9e30-40f6-9c84-d982298087b8
pair_couplings = sort!(random_couplings(10000))

# ╔═╡ 28716ac8-d601-4007-a0d0-1b8d28976091
let fig = Figure()
	ax = Axis(fig[1,1])
	t_wait = 10
	ϵ = 0.02
	couplings = filter(J->π > J*t_wait/4>ϵ*π, pair_couplings)
	hist(couplings)
end
	

# ╔═╡ 1ff8c58e-c681-40b8-aa5c-78ea1cf54e91
with_theme(theme(; height=2, width=3)) do
	ncycles = 1000
	t_wait = 10
	ϵ = 0.02
	ϕ = π*(1-ϵ)
	interaction_dom_couplings = filter(J->J*t_wait/4>ϵ*π, pair_couplings)
	κ = length(interaction_dom_couplings)/length(pair_couplings)

	times = 0:ncycles
	
	fig = Figure()

	ax1 = Axis(fig[1,1])
	hist!(ax1, log10.(pair_couplings); bins=30, normalization=:probability)
	vlines!(ax1, [log10(4*ϵ*π/t_wait)]; color=Makie.wong_colors()[2],
		label=L"J\tau \geq 4\phi")
	vspan!(ax1, log10(4*ϵ*π/t_wait), 4; color=(Makie.wong_colors()[2], 0.3))
	vlines!(ax1, [log10(4/t_wait)]; color=:gray, linestyle=:dash,
		label=L"4/\tau")
	ax1.xlabel = L"pair coupling $J$"
	ax1.ylabel = L"$P_{pair}(b_i \leq J < b_{i+1})$"
	ax1.xticks = -3:1:3
	ax1.xtickformat = xs -> [L"10^{%$(round(Int, x))}" for x in xs]
	ax1.title = "Pair coupling distribution"
	xlims!(ax1, -3.5,3.5)
	axislegend(ax1; position=:rt)
	
	ax2 = Axis(fig[1,2]; xticks=([0,π/6,π/3,π/2], ["0", L"\pi/6", L"\pi/3", L"\pi/2"]))#["0",L"\frac{\pi}{6}",L"\frac{\pi}{3}",L"\frac{\pi}{2}"]))
	hist!(ax2, mod.(interaction_dom_couplings .* t_wait ./ 4, pi/2); 
		bins=30, normalization=:pdf,
		color=Makie.wong_colors()[2])
	hlines!(ax2, [2/π], color=:gray, linestyle=:dot, label=L"\frac{2}{\pi}")
	text!(ax2, 1.4, 2/π; text=L"\frac{2}{\pi}", color=:gray)
	ylims!(ax2, -0.05, 0.9)
	ax2.yticks = 0:0.2:0.8
	ax2.xlabel = L"\bar{J}\ \mathrm{mod}\ \pi/2"
	ax2.title = "Phase-wrapped interactions"
	ax2.ylabel = L"\mathrm{P}(\bar{J}\ \mathrm{mod}\ \pi/2\ |\ \bar{J} > \phi)"
	
	ax = Axis(fig[2,1:2])
	lines!(ax, times, 
		mean(J->single_pair_analytic(;t_wait, ϕ=π*(1-ϵ), J).(times), pair_couplings), color=:grey, alpha=0.5, label="exact")
	lines!(ax, times, 
		mean(J->single_pair_analytic_approx(;t_wait, ϕ=ϵ*π, J).(times), pair_couplings), label="approx. (stat. average)")
	lines!(ax, times, @.(0.5*κ*exp(-times/2*ϵ^2*pi^2)); linestyle=:dash, linewidth=2, label=L"\kappa/2\ \exp(-n\epsilon^2/2)")
	ylims!(ax, -0.025, 0.525)
	ax.yticks = 0:0.1:0.5
	ax.xticks = 0:200:ncycles
	ax.xlabel = L"Cycle $n$"
	ax.ylabel = L"\langle M_z(n) \rangle"
	axislegend(ax; position=:rt)
	text!(ax, 400, 0.4; text=L"d=\alpha,\ \epsilon=%$(round(Int,ϵ*100))%,\ \tau=%$(t_wait)J_{med}", align=(:center,:bottom))
	
	Label(fig[1,1,TopLeft()], "(a)")
	Label(fig[1,2,TopLeft()], "(b)")
	Label(fig[2,1,TopLeft()], "(c)")
	rowgap!(fig.layout, 0)
	fig
end |> save_and_display("pair-model-timecrystal", "part2")

# ╔═╡ 25898d74-976f-4da3-87f8-0490e666428e
H = NN(Chain(9))*Hopp()+ 0.5 .* (-1).^(1:9)*Z()

# ╔═╡ 97a08149-0c82-4beb-825c-034c89e07a29
ψ_neel = let v = zeros(2^9)
	v[1+sum(k->2^k, 1:2:8)] = 1
	v
end

# ╔═╡ bd8c8d2c-47de-407a-be5f-19d0fdb3193c
O2 = -sparse(NN(Chain(9))*ZZ()/16)

# ╔═╡ ed00bbac-74e7-4faf-bfe6-ba3920656861
dot(ψ_neel, O2, ψ_neel)

# ╔═╡ 9e8e1e17-ad21-40f3-a651-e24e84fd9a04
E,U = eigen!(Hermitian(Matrix(H)))

# ╔═╡ a537a230-4991-4c49-a818-8f296b9dad52
ψt(t; E=E,U=U) = U*(Diagonal(cis.(-E .*t))*(U'ψ_neel))

# ╔═╡ 2bd60e33-6c8d-40f7-bb39-f0585bd85aed
lines(range(0,200; length=101), t->real(dot(ψt(t), O2, ψt(t))))

# ╔═╡ 4249337d-691c-4da4-915a-aa4ad52eefde
drive = Hermitian(Matrix(sparse(ones(9)*X())))

# ╔═╡ 82344207-1345-47e5-bf5b-349147e38d36
T=0.001

# ╔═╡ d44cd737-afe6-4033-96c2-c30020c1ead2
U_F = U'*Diagonal(cis.(-E*T))*U*cis(drive*T)

# ╔═╡ 51ad9170-49b0-4b2e-8639-96aa0415e061
E_2, U_2 = eigen(U_F)

# ╔═╡ fe13e6dd-4427-4b2d-b2e1-29d43040bffb
lines(range(0,20; length=101), t->real(dot(ψt(-im*t;E=E_2,U=U_2), O2, ψt(-im*t;E=E_2,U=U_2))))

# ╔═╡ 5cae33c3-f816-49c7-8646-8cc86c87d7a4
md"""
# Gas in Box with divider
"""

# ╔═╡ 9d209080-e01f-4102-a590-84cad949ebe1
with_theme(theme(; height=1, width=3)) do
	fig = Figure()
	ax1 = Axis(fig[1,1]; title=L"t=0")
	hidedecorations!(ax1)
	hidespines!(ax1)
	box = Rect2(0,0,2,1)
	poly!(ax1, box; strokecolor=:black, color=(:white,0.0), strokewidth=2)
	lines!(ax1, [1,1], [0,1]; color=:gray, linestyle=:dot)
	Nparticles = 40
	particle_size = 10
	mindist = 0.03
	rng = Xoshiro(3)
	positions1 = zeros(Nparticles, 2)
	at = 1
	while at <= Nparticles
		newx,newy = rand(rng, 2)
		mindist <= newx <= 1-mindist || continue
		mindist <= newy <= 1-mindist || continue
		positions1[at, :] = [newx,newy]
		at += 1
		for i in 1:at-2
			hypot(newx-positions1[i,1], newy-positions1[i,2]) > 2mindist || (at -= 1; break)
		end
	end
	scatter!(ax1, positions1[:,1], positions1[:,2], markersize=particle_size)

	ax2 = Axis(fig[1,2]; title=L"t\gg1")
	hidedecorations!(ax2)
	hidespines!(ax2)
	poly!(ax2, box; strokecolor=:black, color=(:white,0.0), strokewidth=2)
	positions2 = zeros(Nparticles, 2)
	at = 1
	while at <= Nparticles
		newx,newy = rand(rng, 2)
		newx *= 2
		mindist <= newx <= 2-mindist || continue
		mindist <= newy <= 1-mindist || continue
		positions2[at, :] = [newx,newy]
		at += 1
		for i in 1:at-2
			hypot(newx-positions2[i,1], newy-positions2[i,2]) > 2mindist || (at -= 1; break)
		end
	end
	scatter!(ax2, positions2[:,1], positions2[:,2], markersize=particle_size)

	Label(fig[1,1,TopLeft()]; text="(a)")
	Label(fig[1,2,TopLeft()]; text="(b)")
	
	fig
end |> save_and_display("gas-in-box", "part1")

# ╔═╡ c3bfe179-595d-4faa-8884-edec9c802b24
md"""
# Density - Coupling distribution
"""

# ╔═╡ a729535b-5e7f-419c-b058-702b81701d11
int = PowerLaw(1)

# ╔═╡ 79b8e13f-bcd5-422c-bbca-73db891ad946
geom_density1 = Blockaded(SpinModels.Box(12,[1,1]); blockade=0.025, retries=-1)

# ╔═╡ fbeb058a-811f-4f73-bbe1-49e4e78d478d
geom_density2 = Blockaded(SpinModels.Box(12,[1,1]); blockade=0.1, retries=-1)

# ╔═╡ 58fb9ff8-a186-40a1-a46c-8d73e4504ff3
geom_density3 = Blockaded(SpinModels.Box(12,[1,1]); blockade=0.2 , retries=-1)

# ╔═╡ bfee42e8-07db-495b-8148-6228ab59276b
begin
	rng1 = Xoshiro(1)
	Js1 = mapreduce(_->inv.(vec(maximum(interaction_matrix(int, geom_density1; rng=rng1); dims=1))), vcat, 1:1000)
end

# ╔═╡ 7ff461fb-522f-4457-9dcf-534fcf21ce2a
begin
	rng2 = Xoshiro(2)
	Js2 = mapreduce(_->inv.(vec(maximum(interaction_matrix(int, geom_density2; rng=rng2); dims=1))), vcat, 1:1000)
end

# ╔═╡ 8cfd627b-c0a0-4800-8bc6-b7043928e49c
begin
	rng3 = Xoshiro(3)
	Js3 = mapreduce(_->inv.(vec(maximum(interaction_matrix(int, geom_density3; rng=rng3); dims=1))), vcat, 1:100)
end

# ╔═╡ 816ef5b5-16e9-4a51-8b57-4945feb1c94e
sqrt(2)/pi

# ╔═╡ e1890d15-eeeb-42a2-a37d-f7064f046496
pi/2/sqrt(3)

# ╔═╡ d194dce5-934c-4f1b-9b8b-25fa7d7469e4
a0 =  √(1/12π)

# ╔═╡ a96dcfbe-3895-426b-af1b-1470cd5c42c9
with_theme(theme(; height=2, width=3)) do
	fig = Figure()

	# column headers
	Label(fig[0,3]; text=L"r_b=0.025", tellwidth=false, halign=:center)
	Label(fig[0,2]; text=L"r_b=0.1", tellwidth=false, halign=:center)
	Label(fig[0,1]; text=L"r_b=0.2", tellwidth=false, halign=:center)

	# row labels
	Label(fig[1,0]; text=L"Sample configuration$$", tellheight=false, valign=:center, rotation=π/2)
	Label(fig[2,0]; text=L"NN distances $P(r_{NN})$", tellheight=false, valign=:center, rotation=π/2)
	
	# sample configurations
	ax1 = Axis(fig[1,3]; aspect=1)
	ax2 = Axis(fig[1,2]; aspect=1)
	ax3 = Axis(fig[1,1]; aspect=1)
	hidedecorations!.((ax1,ax2,ax3))
	
	pos1 = positions(geom_density1; rng=Xoshiro(2))
	pos2 = positions(geom_density2; rng=Xoshiro(1))
	pos3 = positions(geom_density3; rng=Xoshiro(3))

	scatter!(ax1, pos1[1,:], pos1[2,:]; markersize=5)
	scatter!(ax2, pos2[1,:], pos2[2,:]; markersize=20)
	scatter!(ax3, pos3[2,:], pos3[1,:]; markersize=40)

	# coupling distributions
	ax32 = Axis(fig[2,1]; xlabel=L"r_{NN}", xticks=0:0.2:0.6)
	ax22 = Axis(fig[2,2]; xlabel=L"r_{NN}", xticks=0:0.2:0.6)
	ax12 = Axis(fig[2,3]; xlabel=L"r_{NN}", xticks=0:0.2:0.6)

	hist!(ax12, Js1; normalization=:pdf, bins=range(0,0.6;length=41))
	hist!(ax22, Js2; normalization=:pdf, bins=range(0,0.6;length=41))
	hist!(ax32, Js3; normalization=:pdf, bins=range(0,0.6;length=41))
	xlims!.((ax12,ax22,ax32), -0.03, 0.63)
	ax32.yticks=0:4:12

	# arrow underneath
	# ax_arrow = Axis(fig[3,1:3])
	# hidedecorations!(ax_arrow)
	# hidespines!(ax_arrow)
	# arrows!(ax_arrow, [0.05], [0], [0.85], [0]; linewidth=4,arrowsize=21, color=:gray)
	# xlims!(ax_arrow, 0, 1)
	# vlines!(ax_arrow, [0.15, 0.505, 0.86]; ymin=0.5, color=:gray)
	# Label(fig[3,0], L"W=\frac{a_0}{r_B}"; tellwidth=false, tellheight=false)
	# text!(ax_arrow, 0.86, 0; text=L"%$(round(a0/0.025;sigdigits=2))", align=(:center, :top))
	
	# minor layout improvements
	rowgap!(fig.layout, Fixed(5))
	colgap!(fig.layout, Fixed(20))
	rowgap!(fig.layout, 1, Fixed(0))
	# rowsize!(fig.layout, 3, Relative(0.1))

	# labels
	Label(fig[1,1,TopLeft()]; text="(a)", tellwidth=false, tellheight=false, alignmode=Mixed(;left=18,bottom=-11))
	Label(fig[1,2,TopLeft()]; text="(b)", tellwidth=false, tellheight=false, alignmode=Mixed(;left=11,bottom=-11))
	Label(fig[1,3,TopLeft()]; text="(c)", tellwidth=false, tellheight=false, alignmode=Mixed(;left=12,bottom=-11))

	Label(fig[2,1,TopLeft()]; text="(d)", tellwidth=false, tellheight=false, alignmode=Mixed(;left=20,bottom=0))
	Label(fig[2,2,TopLeft()]; text="(e)", tellwidth=false, tellheight=false, alignmode=Mixed(;left=12,bottom=0))
	Label(fig[2,3,TopLeft()]; text="(f)", tellwidth=false, tellheight=false, alignmode=Mixed(;left=13,bottom=0))
	fig
end |> save_and_display("disorder-in-experiment", "part1")

# ╔═╡ Cell order:
# ╠═fe94fe12-3e9a-11ef-3a89-7f23e9876f8e
# ╠═be9a4606-ac35-45c2-b86f-04a886e7c6ff
# ╠═a2341a96-f390-4045-a327-53c422074516
# ╟─914f946b-2a6a-4006-9928-d427638a8c75
# ╠═6cded40a-265c-4079-b8e6-a06fd847339f
# ╠═95f064ef-e3d5-4bf1-8200-8c8bfeec3516
# ╠═3ad05f78-e58f-460a-97a2-3aaaeb827468
# ╟─7e8f7f04-7db7-41bd-b256-5f54d1218ff5
# ╠═3fefa6cb-1416-44ff-95f0-61d4fda577ce
# ╠═9a654b24-2a37-47b7-af75-e44a1ab36c16
# ╠═d4f6675a-298b-4df4-9f7b-20cc18a86677
# ╠═a486a8cf-190a-4ec7-8a53-cb736b6791f5
# ╠═a4a21eda-eab2-4697-9dc2-18318ac23892
# ╠═361bf4c3-93b5-4d38-82a5-63348bb37cdc
# ╟─623ebbaa-690c-4251-8397-51c7343d2268
# ╠═3321be84-53da-4d37-82c9-2ea1398b301b
# ╠═7f0514fe-3567-4f3e-8342-646822d74fa4
# ╠═8f3a4e8b-9d92-4101-9611-c76dfe72aa4a
# ╠═9966a21f-090b-4d06-855b-f9bffa50bf6a
# ╟─35d9cf54-fabf-4810-85a7-1978eddfba38
# ╠═fdfdac0a-e91e-405b-aef5-55089c6b080f
# ╠═ffb54136-d1d0-46fa-a126-9b239e08a49d
# ╠═bdd01911-03cb-488e-b180-52e282820901
# ╠═de9d9482-9395-4426-83d6-240f35bfcb53
# ╠═6b1eb731-cf69-4e8d-96df-77a81b404896
# ╠═bb545fe2-0274-4000-8cba-9972e5fb3d19
# ╟─a3fb9d7f-d360-4bb8-8cc7-ee63b40911be
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
# ╟─6044befd-5691-4487-b026-36ab74f77038
# ╠═23ae0570-4ac0-4a47-b75b-ca58f0a2ed93
# ╠═988735ee-9c03-4015-8c7f-544c41de869e
# ╠═6073fa3d-004f-4b2d-a918-7ddfb47157fd
# ╠═b2e0d795-60c9-4cdf-871b-557ff9c5412b
# ╠═f1eb54e4-7a6a-4a31-88c8-4070602c7884
# ╠═7f75e93e-7ec1-4cce-aef6-5f04ea4eeefd
# ╠═cef8c186-f93e-4dca-8418-999327a25afd
# ╟─565161a7-aa79-4f3a-947f-f23826da3850
# ╠═7a26d100-9b5e-4cd0-82b0-1978082805f7
# ╟─25c40d98-d67b-4bc6-994b-c343f3fb06fe
# ╠═a9cbbb3c-ec38-4742-8595-e4953c0a2c97
# ╠═0c362908-f122-4c9b-8ba1-909262cd80b4
# ╟─f9119509-6e24-48ee-a1fb-71b970afc11b
# ╠═d1fd0411-db74-450f-bc08-21719117f99a
# ╠═c6c227e9-0958-4844-8295-9df5055c8cb7
# ╠═b7c8eb86-1182-4325-8930-14ef52cca045
# ╠═c2a2d9f4-9e30-40f6-9c84-d982298087b8
# ╠═28716ac8-d601-4007-a0d0-1b8d28976091
# ╠═1ff8c58e-c681-40b8-aa5c-78ea1cf54e91
# ╠═25898d74-976f-4da3-87f8-0490e666428e
# ╠═97a08149-0c82-4beb-825c-034c89e07a29
# ╠═bd8c8d2c-47de-407a-be5f-19d0fdb3193c
# ╠═ed00bbac-74e7-4faf-bfe6-ba3920656861
# ╠═9e8e1e17-ad21-40f3-a651-e24e84fd9a04
# ╠═a537a230-4991-4c49-a818-8f296b9dad52
# ╠═2bd60e33-6c8d-40f7-bb39-f0585bd85aed
# ╠═4249337d-691c-4da4-915a-aa4ad52eefde
# ╠═82344207-1345-47e5-bf5b-349147e38d36
# ╠═d44cd737-afe6-4033-96c2-c30020c1ead2
# ╠═51ad9170-49b0-4b2e-8639-96aa0415e061
# ╠═fe13e6dd-4427-4b2d-b2e1-29d43040bffb
# ╟─5cae33c3-f816-49c7-8646-8cc86c87d7a4
# ╠═9d209080-e01f-4102-a590-84cad949ebe1
# ╟─c3bfe179-595d-4faa-8884-edec9c802b24
# ╠═a729535b-5e7f-419c-b058-702b81701d11
# ╠═79b8e13f-bcd5-422c-bbca-73db891ad946
# ╠═fbeb058a-811f-4f73-bbe1-49e4e78d478d
# ╠═58fb9ff8-a186-40a1-a46c-8d73e4504ff3
# ╠═bfee42e8-07db-495b-8148-6228ab59276b
# ╠═7ff461fb-522f-4457-9dcf-534fcf21ce2a
# ╠═8cfd627b-c0a0-4800-8bc6-b7043928e49c
# ╠═816ef5b5-16e9-4a51-8b57-4945feb1c94e
# ╠═e1890d15-eeeb-42a2-a37d-f7064f046496
# ╠═d194dce5-934c-4f1b-9b8b-25fa7d7469e4
# ╠═a96dcfbe-3895-426b-af1b-1470cd5c42c9
