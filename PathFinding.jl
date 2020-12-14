function all_combinations(start_genotypes; depth=1)


	d = IdDict{Genotype,Parents}()

	sizehint!(d, 500 * (4* depth)) # account for length of start_genotypes

	# add start genotypes to dict.
	for original in start_genotypes
		d[original] = origin()
	end

	for _ = 1:depth
		
		dkeys = keys(copy(d)) # get available genotypes to cross.

		for key in dkeys
			print(".")
			cross_to_dict!(key, d)
		end

	end
	return d
end

function find_path(start_genotypes, target; depth=1)

	d = IdDict{Genotype,Parents}()

	sizehint!(d, 500 * (4* depth)) # account for length of start_genotypes

	# add start genotypes to dict.
	for original in start_genotypes
		d[original] = origin()
	end

	for current_depth = 1:depth
		
		dkeys = keys(copy(d)) # get available genotypes to cross.

		for key in dkeys
			print(".")
			cross_to_dict!(key, d)
			@infiltrate
			if target in keys(d)
				print("Found the Genotype!")
				return d#target_to_origin(d, target, current_depth)
			end
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

function target_to_origin(crossing_tree, target, depth)

	parents_found = false



	while !parents_found

		# find the parents
		parents = crossing_tree[target]

		# then find the parents' parents.

		# add it to the final crossing dict.
	end
end

function find_origin(dict, target)

	if !is_origin(dict[target].mum)
		print(dict[target].mum)
		find_origin(dict, dict[target].mum)	
	end
	if !is_origin(dict[target].dad)
		print(dict[target].dad)
		find_origin(dict, dict[target].dad)
	end

	return 
end


Base.hash(g::Genotype, h::UInt) = hash(string(g), h)
Base.hash(ch::Chromosome, h::UInt) = hash(ch.genes)
Base.hash(gs::Tuple{Genotype, Genotype}, h::UInt) = hash(string(gs[1]) * "x" * string(gs[2]),h)