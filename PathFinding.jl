function path_exists(start_genotypes; depth=2)
	
	out_array = Vector{Genotype}(undef, 500 * (length(start_genotypes)*depth))
	i = 1

	for g in start_genotypes
		out_array[i] = g
		i+=1
	end
	
	for layer in 1:depth

		new_gens = unique(cross(start_genotypes));

		f1_number = length(new_gens)

		# if our output_array gets too big
		if f1_number + i > length(out_array)
			new_array = Vector{eltype(out_array)}(undef, 2 * length(out_array))
			new_array[1:i] .= out_array[1:i]
			out_array = new_array
		end

		out_array[i:i+length(new_gens)-1] .= new_gens
		i = i+length(new_gens)

		unique_genotypes = unique(out_array[1:i-1])
		i = 1 + length(unique_genotypes)
		out_array[1:i-1] .= unique_genotypes

	end

	return out_array[1:i-1]
end



function f1(g1, g2, target)

	children = cross(g1, g2)
	target ∈ children ? true : children
end

function calc_max_permutations(n, depth)

	(n * 256)^depth	

end