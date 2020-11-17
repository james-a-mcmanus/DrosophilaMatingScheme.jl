function path_exists(start_genotypes; depth=1)


	d = IdDict{Genotype,Parents}()

	sizehint!(d, 500 * (4* depth)) # account for length of start_genotypes

	for original in start_genotypes
		d[original] = origin()
	end


	for _ = 1:depth
		
		dkeys = keys(copy(d))

		for key in dkeys
			print(".")
			cross_to_dict!(key, d)
		end

	end
	return d
end

function cross_to_dict!(g1, d)

	for key in keys(d)
		cross(d, g1, key)
	end
end

function add_children_to_dict(parents, children, d)

	for child in children
		if !haskey(d, child) 
			d[child] = parents
		end
	end
	return d
end


function f1(g1, g2, target)

	children = cross(g1, g2)
	target ∈ children ? true : children
end


Base.hash(g::Genotype, h::UInt) = hash(string(g), h)
Base.hash(gs::Tuple{Genotype, Genotype}, h::UInt) = hash(string(gs[1]) * "x" * string(gs[2]),h)