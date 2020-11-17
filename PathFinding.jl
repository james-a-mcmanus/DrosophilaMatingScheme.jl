function path_exists(start_genotypes, target; depth=2)
	
	#maximum_arraysize = calc_max_permutations(start_genotypes)

	for g in available_genotypes
		Family(g, origin())
	end
	
	for layer in 1:depth

		available_genotypes = vcat(available_genotypes, unique(cross(available_genotypes)))

	end

	return available_genotypes
end

# okay we can't preallocate all possible combinations. This leads to ballooning. e.g. 3 start genotypes and 3 layers = 28Gb of data.


function f1(g1, g2, target)

	children = cross(g1, g2)
	target ∈ children ? true : children
end

function calc_max_permutations(n, depth)

	(n * 256)^depth	

end