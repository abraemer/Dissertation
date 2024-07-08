### A Pluto.jl notebook ###
# v0.19.43

using Markdown
using InteractiveUtils

# ╔═╡ ddcb6926-3c44-11ef-3c7b-654dec1e65df
begin
	import Pkg
	Pkg.activate(".")
	using CairoMakie, PlutoUI, JLD2
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

# ╔═╡ a36490fd-4eef-400d-a2d0-8565923769c1
color_pairs = RGBf(0x08/255,0x08/255,0x64/255)

# ╔═╡ 58a01ffa-091e-410e-8c90-de159f0686f1
color_nn = RGBf(0xBF/255,0x40/255,0x40/255)

# ╔═╡ ec0c2a98-6c05-45dd-85e8-ec648ea6381f
function cusp_plot!(ax; d,α,λ)
	Ωspace = range(-5,5;length=201)
	lines!(ax, Ωspace, nn_coupling_cdf(;d,α,λ); label="NN", color=color_nn)
	lines!(ax, Ωspace, pair_coupling_cdf(;d,α,λ); label="pair", color=color_pairs)
	ax.xticks = -4:2:4
	ax.xlabel = L"Field $\Omega$ [$\lambda^\alpha$]"
	ax.ylabel=L"Magnetization $\langle M_x(\Omega)\rangle$"
end

# ╔═╡ 8032340d-b92e-4301-857c-c5f9cb86304e
200/1.5

# ╔═╡ 461879ec-c7c5-452c-b649-2b352f216a2f
with_theme(theme(;height=1, width=2)) do
	let λ = 1, d=3, α=3,
		fig = Figure(),
		ax1 = Axis(fig[1,1])
		cusp_plot!(ax1; d,α,λ)
		ax1.title=L"\alpha=d=3"
		Label(fig[1,1, Top()], "(a)"; tellwidth=false,tellheight=true, alignmode=Outside(1))
		
		ax2 = Axis(fig[1,2])
		cusp_plot!(ax2; d,α=2α,λ)
		ax2.title=L"\alpha=2d=6"
		ax2.ylabelvisible = false
		ax2.yticklabelsvisible = false
		Label(fig[1,2, Top()], "(b)"; tellwidth=false,tellheight=true, alignmode=Outside(1))
		
		linkyaxes!(ax1,ax2)
		axislegend(ax1; position=:lb)
		fig
	end
end |> save_and_display("analytical_cusp", "part1")

# ╔═╡ 24d49333-72f1-4a53-87b2-03e1fd165f71
with_theme(theme(;height=0, scale=2)) do
	let λ = 1, d=3, α=3,
		fig = Figure(;size=(50+250, 50+200*2), pt_per_unit=1),
		ax1 = Axis(fig[1,1])
		cusp_plot!(ax1; d,α,λ)
		ax1.title=L"\alpha=d=3"
		Label(fig[1,1, Top()], "(a)"; tellwidth=false,tellheight=true, alignmode=Outside(1))
		ax1.xlabelvisible=false
		ax1.xticklabelsvisible=false
		
		ax2 = Axis(fig[2,1])
		cusp_plot!(ax2; d,α=2α,λ)
		ax2.title=L"\alpha=2d=6"
		# ax2.ylabelvisible = false
		# ax2.yticklabelsvisible = false
		Label(fig[2,1, Top()], "(b)"; tellwidth=false,tellheight=true, alignmode=Outside(1))
		
		linkyaxes!(ax1,ax2)
		linkxaxes!(ax1,ax2)
		axislegend(ax1; position=:lb)
		fig
	end
end

# ╔═╡ Cell order:
# ╠═ddcb6926-3c44-11ef-3c7b-654dec1e65df
# ╠═f03b3dae-b7e3-4b33-82fa-088960d23e15
# ╠═8c36cedd-94b8-4761-ab73-fe00e6ed6310
# ╠═a36490fd-4eef-400d-a2d0-8565923769c1
# ╠═58a01ffa-091e-410e-8c90-de159f0686f1
# ╠═ec0c2a98-6c05-45dd-85e8-ec648ea6381f
# ╠═8032340d-b92e-4301-857c-c5f9cb86304e
# ╠═461879ec-c7c5-452c-b649-2b352f216a2f
# ╠═24d49333-72f1-4a53-87b2-03e1fd165f71
