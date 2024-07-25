### A Pluto.jl notebook ###
# v0.19.43

using Markdown
using InteractiveUtils

# ╔═╡ ddcb6926-3c44-11ef-3c7b-654dec1e65df
begin
	import Pkg
	Pkg.activate(".")
	using CairoMakie, PlutoUI, JLD2, Integrals
	include("helpers.jl")
end

# ╔═╡ f03b3dae-b7e3-4b33-82fa-088960d23e15
begin
	function nn_coupling_cdf(; d, α, λ)
		return J->nn_coupling_cdf(J; d, α, λ)
	end
	function nn_coupling_cdf(J; d, α, λ)
		β = d/α
		return @. 0.5*exp(-λ^d * abs(J)^(-β))
	end
end

# ╔═╡ 39734bc6-e4e6-406c-9b7e-151699187b7e
begin
	function nn_coupling_pdf(; d, α, λ)
		return J->nn_coupling_pdf(J; d, α, λ)
	end
	function nn_coupling_pdf(J; d, α, λ)
		β = d/α
		return @. 0.5*β*λ^d * abs(J)^(-β-1)*exp(-λ^d * abs(J)^(-β))
	end
end

# ╔═╡ 8c36cedd-94b8-4761-ab73-fe00e6ed6310
begin
	function pair_coupling_cdf(; d, α, λ)
		return J->pair_coupling_cdf(J; d, α, λ)
	end
	function pair_coupling_cdf(J; d, α, λ)
		β = d/α
		return @. 0.5/(1+λ^d*abs(J)^(-β))
	end
end

# ╔═╡ ac404ac2-6351-4601-8fd9-a277a58ac0e9
begin
	function pair_coupling_pdf(; d, α, λ)
		return J->pair_coupling_cdf(J; d, α, λ)
	end
	function pair_coupling_pdf(J; d, α, λ)
		β = d/α
		return @. 0.5β*λ^d*abs(J)^(-β-1)/(1+λ^d*abs(J)^(-β))^2
	end
end

# ╔═╡ a36490fd-4eef-400d-a2d0-8565923769c1
color_pairs = RGBf(0x08/255,0x08/255,0x64/255)

# ╔═╡ 58a01ffa-091e-410e-8c90-de159f0686f1
color_nn = RGBf(0xBF/255,0x40/255,0x40/255)

# ╔═╡ ec0c2a98-6c05-45dd-85e8-ec648ea6381f
function cusp_plot!(ax; d,α,λ,Ωspace = range(-5,5;length=201))
	lines!(ax, Ωspace, nn_coupling_cdf(;d,α,λ); label=L"M_{NN}", color=color_nn)
	lines!(ax, Ωspace, pair_coupling_cdf(;d,α,λ); label=L"M_{pair}", color=color_pairs)
	ax.xticks = -4:2:4
	ax.xlabel = L"Field $\Omega$ [$\lambda^\alpha$]"
	ax.ylabel=L"Magnetization $\langle M(\Omega)\rangle$"
end

# ╔═╡ 8032340d-b92e-4301-857c-c5f9cb86304e
200/1.5

# ╔═╡ 461879ec-c7c5-452c-b649-2b352f216a2f
with_theme(theme(;height=1, width=2)) do
	let λ = 1, d=3,
		fig = Figure()
		
		ax1 = Axis(fig[1,1])
		cusp_plot!(ax1; d,α=d,λ)
		ax1.title=L"\alpha=d=%$d"
		Label(fig[1,1, Top()], "(a)"; tellwidth=false,tellheight=true, alignmode=Outside(1))
		
		ax2 = Axis(fig[1,2])
		cusp_plot!(ax2; d,α=2d,λ)
		ax2.title=L"\alpha=2d=%$(2d)"
		ax2.ylabelvisible = false
		ax2.yticklabelsvisible = false
		Label(fig[1,2, Top()], "(b)"; tellwidth=false,tellheight=true, alignmode=Outside(1))
		
		linkyaxes!(ax1,ax2)
		axislegend(ax1; position=:lb)
		fig
	end
end |> save_and_display("analytical_cusp", "part1")

# ╔═╡ 66d57495-4670-4a57-8c6e-c7c433b7cd17
with_theme(theme(;height=1, width=2)) do
	let λ = 1, d=3, ir = 0.2, # inset range
		fig = Figure()

		ax1 = Axis(fig[1,1])
		cusp_plot!(ax1; d,α=d,λ)
		ax1.title=L"\alpha=d=%$d"
			
		ax2 = Axis(fig[1,2])
		ax2.title=L"\alpha=2d=%$d"

		ax_in = Axis(fig; bbox=BBox(290,400,150,190))
		cusp_plot!(ax_in; d,α=2d,λ, Ωspace=range(-ir,ir;length=201))
		ax_in.xticks = -0.2:0.2:0.2
		ax_in.xminorticks = -0.2:0.05:0.2
		ax_in.xminorticksvisible = true
		ax_in.xminortickalign = 1
		ax_in.xgridvisible = false
		ax_in.xlabelvisible = false
		ax_in.ylabelvisible = false
		ax_in.yticks = -0.0:0.1:0.2
		ax_in.yminorticks = 0:0.025:0.2
		ax_in.yminorticksvisible = true
		ax_in.yminortickalign = 1
		ax_in.ygridvisible = false
		#translate!(ax_in.scene, 0,0,10)
		translate!(ax_in.blockscene, 0,0,100)
		coords = lift(ax_in.xaxis.attributes.limits, ax_in.yaxis.attributes.limits) do xlims, ylims
			x1,x2 = xlims
			y1,y2 = ylims
			Point2f[(x1,y1),(x2,y1),(x2,y2),(x1,y2)]
		end
		poly!(ax2, coords; color=(:grey,0.7), strokecolor=:grey, strokewidth=1)
		cusp_plot!(ax2; d,α=2d,λ)

		xlims!(ax2, -5.1, 5.1)
		ylims!(ax2, -0.02, 0.55)
		ax2.ylabelvisible = false
		ax2.yticklabelsvisible = false
		ax1.yticks = 0.0:0.1:0.5
		ax2.yticks = 0.0:0.1:0.5
		
		linkyaxes!(ax1,ax2)
		axislegend(ax1; position=:lb)

		Label(fig[1,2, Top()], "(b)"; tellwidth=false,tellheight=true, alignmode=Outside(1))
		Label(fig[1,1, Top()], "(a)"; tellwidth=false,tellheight=true, alignmode=Outside(1))
		fig
	end
end  |> save_and_display("analytical_cusp_inset", "part1")

# ╔═╡ 40288a5e-1cdb-43d1-9a4d-c391c7fc5ddf
with_theme(theme(;height=1, width=2)) do
	let λ = 1, d=3, ir = 0.2, # inset range
		fig = Figure()

		ax1 = Axis(fig[1,1])
		cusp_plot!(ax1; d,α=1,λ, Ωspace=range(-2,2;length=201))
		ax1.title=L"\alpha=1, d=%$d"
			
		ax2 = Axis(fig[1,2])
		ax2.title=L"\alpha=2, d=%$d"

		ax_in = Axis(fig; bbox=BBox(290,400,150,190))
		cusp_plot!(ax_in; d,α=2,λ, Ωspace=range(-ir,ir;length=201))
		ax_in.xticks = -0.2:0.2:0.2
		ax_in.xminorticks = -0.2:0.05:0.2
		ax_in.xminorticksvisible = true
		ax_in.xminortickalign = 1
		ax_in.xgridvisible = false
		ax_in.xlabelvisible = false
		ax_in.ylabelvisible = false
		ax_in.yticks = -0.0:0.1:0.2
		ax_in.yminorticks = 0:0.025:0.2
		ax_in.yminorticksvisible = true
		ax_in.yminortickalign = 1
		ax_in.ygridvisible = false
		#translate!(ax_in.scene, 0,0,10)
		translate!(ax_in.blockscene, 0,0,100)
		coords = lift(ax_in.xaxis.attributes.limits, ax_in.yaxis.attributes.limits) do xlims, ylims
			x1,x2 = xlims
			y1,y2 = ylims
			Point2f[(x1,y1),(x2,y1),(x2,y2),(x1,y2)]
		end
		poly!(ax2, coords; color=(:grey,0.7), strokecolor=:grey, strokewidth=1)
		cusp_plot!(ax2; d,α=2,λ, Ωspace=range(-2,2;length=201))

		xlims!(ax2, -2.2, 2.2)
		ylims!(ax2, -0.02, 0.55)
		ax2.ylabelvisible = false
		ax2.yticklabelsvisible = false
		ax1.yticks = 0.0:0.1:0.5
		ax2.yticks = 0.0:0.1:0.5
		ax1.xticks = -2:1:2
		ax2.xticks = -2:1:2
		
		linkyaxes!(ax1,ax2)
		axislegend(ax1; position=:ct)

		Label(fig[1,2, Top()], "(b)"; tellwidth=false,tellheight=true, alignmode=Outside(1))
		Label(fig[1,1, Top()], "(a)"; tellwidth=false,tellheight=true, alignmode=Outside(1))
		fig
	end
end  |> save_and_display("analytical_cusp_alpha_small", "part1")

# ╔═╡ 5c354f0f-51b3-4b7d-a4fd-bb3fc15f3126
lorentz(ΩJ) = ΩJ^2/(1+ΩJ^2)

# ╔═╡ f75bc96a-a32d-42dd-8c1f-624065826108
theta(ΩJ) = abs2(ΩJ) > 1

# ╔═╡ 2a93c508-ca64-40dc-a797-6920f160d031
f,a,l = lines(range(-2,2;length=201), lorentz); lines!(a, range(-2,2;length=201), theta); f

# ╔═╡ ccb3dce2-4fe1-4786-87ae-40be9d7c413d
function theta_pair(Ω; d=3,α=3,λ=1, cutoff = 1e-10)
	function integrand(J, _)
		ΩJ = abs(Ω)/J
		ΩJ < cutoff && return 0.0
		return pair_coupling_pdf(J;d,α,λ)
	end
	abs(Ω) < cutoff && return 0
	sol = solve(
		IntegralProblem(integrand, (0.0,abs(Ω))),
		QuadGKJL();
		reltol=1e-4,abstol=1e-4)
	return sol.u
end

# ╔═╡ cbd42f68-bdcf-4113-86e8-6b2e2e6da5cc
function lorentz_pair(Ω; d=3,α=3,λ=1, cutoff = 1e-6)
	function integrand!(J, p)
		ΩJ = abs(Ω)/J
		ΩJ < cutoff && return 0.0
		return lorentz(ΩJ) * pair_coupling_pdf(J;d,α,λ)
	end
	sol = solve(
		IntegralProblem(integrand!, (0.0,Inf)),
		QuadGKJL();
		reltol=1e-5,abstol=1e-5)
	return sol.u
end

# ╔═╡ 4233a7a8-e425-4def-ad78-ac4fd5f8f7a4
function lorentz_nn(Ω; d=3,α=3,λ=1, cutoff = 1e-6)
	function integrand!(J, p)
		ΩJ = abs(Ω)/J
		ΩJ < cutoff && return 0.0
		return lorentz(ΩJ) * nn_coupling_pdf(J;d,α,λ)
	end
	sol = solve(
		IntegralProblem(integrand!, (0.0,Inf)),
		QuadGKJL();
		reltol=1e-5,abstol=1e-5)
	return sol.u
end

# ╔═╡ 8c439262-765d-47a0-8a5d-bf5cd133a883
lorentz_pair(0.5), theta_pair(0.5)

# ╔═╡ c61638c6-753f-4afa-b9a5-0cff1c218999
let
	f = Figure()
	a = Axis(f[1,1])
	lines!(a, range(-4,4;length=201), pair_coupling_cdf(;d=3,λ=1,α=3))
	lines!(a, range(-4,4;length=201), lorentz_pair)
	lines!(a, range(-4,4;length=201), nn_coupling_cdf(;d=3,λ=1,α=3))
	lines!(a, range(-4,4;length=201), lorentz_nn)
	f
end

# ╔═╡ 26d8c657-397d-4c0b-839f-b1e1721e1e05
with_theme(theme(;height=1, width=2)) do
	let λ = 1, d=3, ir = 0.2, # inset range
		fig = Figure()

		ax1 = Axis(fig[1,1])
		lines!(ax1, range(0,3;length=201), x->theta(x)/2; label=L"Heaviside $\Theta$", color=:grey, linestyle=:dash)
		lines!(ax1, range(0,3;length=201), x->lorentz(x)/2; label="Lorentzian", color=:grey);
		ax1.xlabel=L"|\Omega/J|"
		ax1.ylabel="Magnetization"
		ax1.title="Activation function"
		ax1.titlefont = :regular
			
		ax2 = Axis(fig[1,2])
		ax2.title=L"\alpha=d=%$d"
		ax2.xlabel=L"Field $\Omega$ [$\lambda^\alpha$]"
		lines!(ax2, range(-1,1;length=201), nn_coupling_cdf(;d,α=d,λ=1), linestyle=:dash, color=color_nn)
		lines!(ax2, range(-1,1;length=201), x->lorentz_nn(x;d,α=d,λ=1), color=color_nn, label=L"M_{NN}")
		lines!(ax2, range(-1,1;length=201), pair_coupling_cdf(;d,α=d,λ=1), linestyle=:dash, color=color_pairs)
		lines!(ax2, range(-1,1;length=201), x->lorentz_pair(x;d,α=d,λ=1), color=color_pairs, label=L"M_{pair}")

		# ax_in = Axis(fig; bbox=BBox(290,400,150,190))
		# ax_in.xticks = -0.2:0.2:0.2
		# ax_in.xminorticks = -0.2:0.05:0.2
		# ax_in.xminorticksvisible = true
		# ax_in.xminortickalign = 1
		# ax_in.xgridvisible = false
		# ax_in.xlabelvisible = false
		# ax_in.ylabelvisible = false
		# ax_in.yticks = -0.0:0.1:0.2
		# ax_in.yminorticks = 0:0.025:0.2
		# ax_in.yminorticksvisible = true
		# ax_in.yminortickalign = 1
		# ax_in.ygridvisible = false
		# #translate!(ax_in.scene, 0,0,10)
		# translate!(ax_in.blockscene, 0,0,100)
		# coords = lift(ax_in.xaxis.attributes.limits, 
		# 		ax_in.yaxis.attributes.limits) do xlims, ylims
		# 	x1,x2 = xlims
		# 	y1,y2 = ylims
		# 	Point2f[(x1,y1),(x2,y1),(x2,y2),(x1,y2)]
		# end
		# poly!(ax2, coords; color=(:grey,0.7), strokecolor=:grey, strokewidth=1)
		# lines!(ax_in, range(-ir,ir;length=201), x->theta_pair(x;d,α=d,λ=1))
		# lines!(ax_in, range(-ir,ir;length=201), x->lorentz_pair(x;d,α=d,λ=1))
		# lines!(ax_in, range(-ir,ir;length=201), x->lorentz_nn(x;d,α=d,λ=1))

		# xlims!(ax2, -5.1, 5.1)
		# ylims!(ax2, -0.02, 0.65)
		ax2.ylabelvisible = false
		#ax2.yticklabelsvisible = false
		ax1.yticks = 0.0:0.1:0.5
		ax2.yticks = 0.0:0.1:0.5
		
		#linkyaxes!(ax1,ax2)

		axislegend(ax1; position=:rb)
		axislegend(ax2; position=:ct, padding=(6,6,3,3), rowgap=3)
		#Legend(fig[2,:], ax1; orientation=:horizontal)

		Label(fig[1,2, Top()], "(b)"; tellwidth=false,tellheight=true, alignmode=Outside(1))
		Label(fig[1,1, Top()], "(a)"; tellwidth=false,tellheight=true, alignmode=Outside(0,0,4,0))
		fig
	end
end |> save_and_display("analytical_cusp_comparison_lorentz", "part1")

# ╔═╡ Cell order:
# ╠═ddcb6926-3c44-11ef-3c7b-654dec1e65df
# ╠═f03b3dae-b7e3-4b33-82fa-088960d23e15
# ╠═39734bc6-e4e6-406c-9b7e-151699187b7e
# ╠═8c36cedd-94b8-4761-ab73-fe00e6ed6310
# ╠═ac404ac2-6351-4601-8fd9-a277a58ac0e9
# ╠═a36490fd-4eef-400d-a2d0-8565923769c1
# ╠═58a01ffa-091e-410e-8c90-de159f0686f1
# ╠═ec0c2a98-6c05-45dd-85e8-ec648ea6381f
# ╠═8032340d-b92e-4301-857c-c5f9cb86304e
# ╠═461879ec-c7c5-452c-b649-2b352f216a2f
# ╠═66d57495-4670-4a57-8c6e-c7c433b7cd17
# ╠═40288a5e-1cdb-43d1-9a4d-c391c7fc5ddf
# ╠═5c354f0f-51b3-4b7d-a4fd-bb3fc15f3126
# ╠═f75bc96a-a32d-42dd-8c1f-624065826108
# ╠═2a93c508-ca64-40dc-a797-6920f160d031
# ╠═ccb3dce2-4fe1-4786-87ae-40be9d7c413d
# ╠═cbd42f68-bdcf-4113-86e8-6b2e2e6da5cc
# ╠═4233a7a8-e425-4def-ad78-ac4fd5f8f7a4
# ╠═8c439262-765d-47a0-8a5d-bf5cd133a883
# ╠═c61638c6-753f-4afa-b9a5-0cff1c218999
# ╠═26d8c657-397d-4c0b-839f-b1e1721e1e05
