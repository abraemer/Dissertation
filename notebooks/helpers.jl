# Figure themeing
HEIGHT = 250
SCALE = 2
bgcolor = RGBf(1,1,1)

const THEME = merge(theme_latexfonts(),
	Theme(fontsize=9*SCALE, size=(246*SCALE,HEIGHT*SCALE), pt_per_unit=1/SCALE,
		figure_padding=(1,7,1,1),
		Axis=(; xtickalign=1, ytickalign=1, xscale=identity,
			backgroundcolor=bgcolor),
		Legend=(; backgroundcolor=bgcolor),
		Label=(; font=:bold,
			halign=:left, valign=:top)))

# save figures
function save_and_display(name, folder="")
	return fig -> save_and_display(name, fig, folder)
end
function save_and_display(name, fig, folder)
	mkpath(joinpath("../gfx", folder))
	mkpath(joinpath("../gfx", folder, "png"))
	Makie.save(joinpath("../gfx", folder, name*".pdf"), fig)
	Makie.save(joinpath("../gfx", folder, "png", name*".png"), fig)
	fig
end

# loading/saving data

function load(name; folder="")
    path = joinpath("data/", folder)
	name = endswith(name, ".jld2") ? name : name*".jld2"
	fullpath = joinpath(path, name)
    return jldopen(fullpath, "r")["data"]
end

function load_or_generate(f, name; folder="")
	path = joinpath("data/", folder)
	name = endswith(name, ".jld2") ? name : name*".jld2"
	isdir(path) || mkpath(path)
	fullpath = joinpath(path, name)
	isfile(fullpath) && return jldopen(fullpath, "r")["data"]
	data = f()
	jldsave(fullpath; data)
	return data
end
macro load_or_generate(expr)
	expr.head != :(=) && error("@load_or_generate is only applicable to simple assignments of the form `a = b`")
	var, rhs = expr.args
	#eqsign == := 
	return quote
		$(esc(var)) = $(esc(:load_or_generate))($(string(var)); folder=basename(replace(@__FILE__, r"\.jl#==#.*" => ""))) do
			return $(esc(rhs))
		end
	end
end
md"""Setup save to $(joinpath(homedir(), "results/notebooks", NOTEBOOKNAME))"""
