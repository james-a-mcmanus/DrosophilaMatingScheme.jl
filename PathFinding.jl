function path_exists(start_genotypes; depth=2)
	
	#maximum_arraysize = calc_max_permutations(start_genotypes)

	out_array = Vector{Family}(undef, 500 * (length(start_genotypes)*depth))
	i = 1

	for g in start_genotypes
		out_array[i] = Family(g, origin())
		i+=1
	end
	
	for layer in 1:depth

		new_gens = unique(cross(start_genotypes));
		out_array[i:i+length(new_gens)-1] .= new_gens
		i = i+length(new_gens)

	end

	return out_array
end

# okay we can't preallocate all possible combinations. This leads to ballooning. e.g. 3 start genotypes and 3 layers = 28Gb of data.


function f1(g1, g2, target)

	children = cross(g1, g2)
	target ∈ children ? true : children
end

function calc_max_permutations(n, depth)

	(n * 256)^depth	

end