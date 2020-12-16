import Base: show, length, isequal, parse, hash, string, unique, ==, length, unique

include("Structs.jl")
include("PathFinding.jl")

"""
Parse strings into genetic types
"""
function Base.parse(::Type{Chromosome}, n, s::AbstractString)

	isempty(s) && return wildtype(Chromosome, n)

	allele_names = (lstrip.(split(s, ",")))
	alleles = Vector{Allele}(undef, length(allele_names))
	for (i, al) in enumerate(allele_names)
		lethality = al[end] == '!' 
		alleles[i] = Allele(al[1:end-lethality], lethality)
	end

	Chromosome(n, alleles) # splatting is slow, avoid parsing if performance critical
end
function Base.parse(::Type{Couple}, n, s::AbstractString)

	isempty(s) && return wildtype(Couple, n)

	chrom_strings = split(s, "/")

	if length(chrom_strings) == 1
		return Couple(Chromosome(n, chrom_strings[1]))
	elseif length(chrom_strings) == 2
		return Couple(Chromosome(n, chrom_strings[1]), Chromosome(n, chrom_strings[2]))
	else
		error("More than 1 divisor in a chromosome?")
	end
end
function Base.parse(::Type{Genotype}, s::String)

	s = replace(s, " " => "")
	chromosomes = split(s, ";")

	length(chromosomes) != 4 && error("There should be 3 semicolons in a genotype definition.")

	Genotype(Couple(1, chromosomes[1]), Couple(2, chromosomes[2]), Couple(3, chromosomes[3]), Couple(4, chromosomes[4]))
end


"""
Change representation to string. 
"""
Base.string(a::Allele) = a.name
function Base.string(c::Chromosome)

	out = ""

	for g in c.genes
	
		out = out * string(g) * ", "

	end

	return out[1:end-2]
end
Base.string(p::Couple) = string(p.c1) * "/" * string(p.c2)
function Base.string(g::Genotype)
	is_origin(g) ? "" : string(g.p1) * " ; " * string(g.p2) * " ; " * string(g.p3) * " ; " * string(g.p4)
end

"""
Return all combinations of a genetic cross
"""
function cross(d::AbstractDict{Genotype,Parents}, g1::Genotype, g2::Genotype)

	chr1 = cross(g1.p1, g2.p1)
	chr2 = cross(g1.p2, g2.p2)
	chr3 = cross(g1.p3, g2.p3)
	chr4 = cross(g1.p4, g2.p4)	
	permute_chromosomes(d, chr1, chr2, chr3, chr4, Parents(g1, g2))
end
function cross(g1::Genotype, g2::Genotype)::Vector{Genotype}

	chr1 = cross(g1.p1, g2.p1)
	chr2 = cross(g1.p2, g2.p2)
	chr3 = cross(g1.p3, g2.p3)
	chr4 = cross(g1.p4, g2.p4)

	return unique(permute_chromosomes(chr1, chr2, chr3, chr4))
end
function cross(p1::Couple{N}, p2::Couple{N}) where {N}

	if homozygous(p1) && homozygous(p2)
		return [Couple(p1.c1, p2.c1)]
	elseif homozygous(p1)
		return [Couple(p1.c1, p2.c1), Couple(p1.c1, p2.c2)]
	elseif homozygous(p2)
		return [Couple(p1.c1, p2.c1), Couple(p1.c2, p2.c1)]
	else
		return [Couple(p1.c1, p2.c1), Couple(p1.c1, p2.c2), Couple(p1.c2, p2.c1), Couple(p1.c2, p2.c2)]
	end
end
cross(s1::AbstractString, s2::AbstractString) = cross(Genotype(s1), Genotype(s2))
function cross(g1::Genotype, genotypes::Vector{Genotype})

	out = Genotype[]

	for g2 in genotypes
		out = vcat(out, cross(g1, g2))
	end
	return out
end
cross(genotypes::Vector{Genotype}, g1::Genotype) = cross(g1, genotypes)
function cross(genotypes::Vector{Genotype})

	out = Genotype[]

	for g1 in genotypes
		out = vcat(out, cross(g1, genotypes))
	end
	return out
end





function cross!(genotypes::Vector{Genotype}, i)

	b = copy(i)

	for g1 in genotypes[1:b-1]
	
		new_genotypes = cross(g1, genotypes[1:i-1])
		new_genotypes = unique(new_genotypes)
		new_genotypes = unique_genotypes(genotypes[i:i-1], new_genotypes)

		if i + length(new_genotypes) > length(genotypes)
			@infiltrate
			new_array = Vector{Genotype}(undef, 10 * length(genotypes))
			new_array[1:i-1] .= genotypes[1:i-1]
			genotypes = new_array
		end

		genotypes[i:i+length(new_genotypes)-1] .= new_genotypes
		i += length(new_genotypes)
	end
	return (genotypes, i)
end

function cross!(genotypes::Vector{Family}, i)

	b = copy(i)

	for g1 in genotypes[1:b-1]
	
		new_genotypes = cross(g1, genotypes[1:i-1])
		@infiltrate
		new_genotypes = unique(new_genotypes)
		@infiltrate
		new_genotypes = unique_genotypes(genotypes[1:i-1], new_genotypes)

		if i + length(new_genotypes) > length(genotypes)
			new_array = Vector{Family}(undef, 10 * length(genotypes))
			new_array[1:i-1] .= genotypes[1:i-1]
			genotypes = new_array
		end

		genotypes[i:i+length(new_genotypes)-1] .= new_genotypes
		i += length(new_genotypes)
	end
	return (genotypes, i)

end

cross(f1::Family, families::Vector{Family}) = cross(families, f1)
function cross(families::Vector{Family}, f1::Family)

	# we can massively improve this by predicting how big the output is going to be.
	out = Family[]
	# cross every family in "families" to "f1"
	
	for f2 in families

		@infiltrate
		new_genotypes = cross(f2.child, f1.child)
		out_families = Vector{Family}(undef, length(new_genotypes))
		for i in 1:length(new_genotypes)
			out_families[i] = Family(new_genotypes[i], Parents(f1.child, f2.child))
		end

		out = vcat(out, out_families)
	end
	return out
end



function unique_genotypes(old_genotypes, new_genotypes)

	uniques = Vector{eltype(old_genotypes)}(undef, length(new_genotypes))
	i = 1
	for g in new_genotypes
		if !(g in old_genotypes)
			uniques[i] = g
			i+=1
		end

	end
	return uniques[1:i-1]
end
function Base.unique(fms::Vector{Family})

	out = Vector{Family}(undef, length(fms))
	i=1
	for fam in fms
		if !(fam.child in children(out[1:i-1]))
			out[i] = fam
			i+=1
		end
	end
	return out[1:i-1]
end

# only considers families unique if they have different children (parents don't matter)
function unique_genotypes(old_genotypes::Vector{Family}, new_genotypes::Vector{Family})


	old_children = children(old_genotypes)
	uniques = Vector{eltype(old_genotypes)}(undef, length(new_genotypes))
	i = 1
	for g in new_genotypes	
		if !(g.child in old_children)
			uniques[i] = g
			i+=1
		end
	end
	return uniques[1:i-1]
end
"""
Takes tuples of chromosomes, combines them into all possible genotypes.
"""
#function permute_chromosomes(cpl1::NTuple{4, Couple{1}}, cpl2::NTuple{4, Couple{2}}, cpl3::NTuple{4, Couple{3}}, cpl4::NTuple{4, Couple{4}})
function permute_chromosomes(cpl1, cpl2, cpl3, cpl4)	
	out = Vector{Genotype}(undef, 256)
	i = 1
	for pr1 in cpl1
		for pr2 in cpl2
			for pr3 in cpl3
				for pr4 in cpl4
					g = Genotype(pr1, pr2, pr3, pr4)
					if !islethal(g)
						out[i] = g
						i+=1
					end
				end
			end
		end
	end
	return out[1:i-1]
end
function permute_chromosomes(d::AbstractDict{Genotype, Parents}, cpl1, cpl2, cpl3, cpl4, parents::Parents)

	for pr1 in cpl1
		for pr2 in cpl2
			for pr3 in cpl3
				for pr4 in cpl4
					g = Genotype(pr1, pr2, pr3, pr4)
					if !haskey(d, g)  && !islethal(g)
						d[g] = parents;
					end
				end
			end
		end
	end
end


"""
Custom printing of genetic types.
"""
longprint(al::Allele) = print(al.name)
function longprint(ch::Chromosome)
	for i in ch.genes
		show(i)
		print(", ")
	end
end
function longprint(cpl::Couple)

	show(cpl.c1)
	print("\n")
	println(Char(8212)^max(namelength(cpl.c1),namelength(cpl.c2)))
	show(cpl.c2)
end
function longprint(gn::Genotype)

	spacing = 3

	print("\n\n")

	print_of_pair(gn.p1, gn.p1.c1, spacing)
	print_of_pair(gn.p2, gn.p2.c1, spacing)
	print_of_pair(gn.p3, gn.p3.c1, spacing)
	print_of_pair(gn.p4, gn.p4.c1, spacing, semi=false)

	# Print the dividing line
	print("\n")
	print(Char(8212)^namelength(gn.p1), " "^(2*spacing + 1))
	print(Char(8212)^namelength(gn.p2), " "^(2*spacing + 1))
	print(Char(8212)^namelength(gn.p3), " "^(2*spacing + 1))
	print(Char(8212)^namelength(gn.p4), " "^(2*spacing + 1))
	print("\n")

	print_of_pair(gn.p1, gn.p1.c2, spacing)
	print_of_pair(gn.p2, gn.p2.c2, spacing)
	print_of_pair(gn.p3, gn.p3.c2, spacing)
	print_of_pair(gn.p4, gn.p4.c2, spacing, semi=false)
	print("\n\n")
end


Base.show(io::IO, ::MIME"text/plain", gn::Genotype) = show(io, gn)
function Base.show(io::IO, gn::Genotype)

	print("[" * string(gn) * "]")
end
Base.show(io::IO, ::MIME"text/plain", fm::Family) = show(io, fm)
function Base.show(io::IO, fm::Family)

	show(fm.parents)
	print(" = ")
	show(fm.child)
end
function Base.show(io::IO, pr::Parents)
	show(pr.mum)
	print(" x ")
	show(pr.dad)
end
Base.show(io::IO, ::MIME"text/plain", vfm::Vector{Family}) = show(io, vfm)
function Base.show(io::IO, vfm::Vector{Family})

	for fm in vfm
		show(fm)
		println()
	end
end


""" 
Equivalence Functions
"""
Base.isequal(a1::Allele, a2::Allele) = a1.name == a2.name
Base.isequal(chr1::Chromosome{N}, chr2::Chromosome{N}) where {N} = issetequal(chr1.genes, chr2.genes)
Base.isequal(p1::Couple{N}, p2::Couple{N}) where {N} = (p1.c1 == p2.c2 && p1.c2 == p2.c1) || (p1.c1 == p2.c1 && p1.c2 == p2.c2)
Base.isequal(g1::Genotype, g2::Genotype) = isequal(g1.p1, g2.p1) && isequal(g1.p2, g2.p2) && isequal(g1.p3, g2.p3) && isequal(g1.p4, g2.p4)


(==)(a1::Allele, a2::Allele) = isequal(a1,a2)
(==)(chr1::Chromosome{N}, chr2::Chromosome{N}) where {N} = isequal(chr1, chr2)
(==)(p1::Couple{N}, p2::Couple{N}) where {N} = isequal(p1,p2)
(==)(g1::Genotype, g2::Genotype) = isequal(g1,g2)


"""
Lethal()
Find out if a genotype is lethal.
"""
islethal(g::Genotype) = islethal(g.p1) || islethal(g.p2) || islethal(g.p3) || islethal(g.p4)
#islethal(p::Couple) = any(islethal.(same_genes(p)))
islethal(a::Allele) = a.lethality

@inline function islethal(p::Couple)

	for ch1_gene in p.c1.genes 
		(ch1_gene.lethality && ch1_gene in p.c2.genes) && return true
	end

	for ch2_gene in p.c2.genes
		(ch2_gene.lethality && ch2_gene in p.c1.genes) && return true
	end
	return false
end

"""
same_genes()
Returns which alleles are the same between two chromosomes.
"""
same_genes(p::Couple) = intersect(p.c1.genes, p.c2. genes)
#p.c1.genes[findall(in(p.c1.genes),p.c2.genes)]


"""
Some helper functions.
"""
function example1()
	Genotype("sn/FM7a!;sp/CyO!;ser/TM3-sb!;")
end
function example2()
	Genotype("w[+]; If/CyO! ; MKRS!/TM6b! ;")
end
function example3()
	Genotype("w[+]; If/CyO! ; mCD8_GFP /TM6b! ;")
end



"""
How long a component will be when printed
"""
namelength(a::Allele) = length(a.name)
namelength(ch::Chromosome) = sum([namelength(g) + 2 for g in ch.genes])
namelength(pr::Couple) = max(namelength(pr.c1), namelength(pr.c2)) -1
namelength(gn::Genotype) = namelength(gn.p1) + namelength(gn.p2) + namelength(gn.p3) + namelength(gn.p4)


"""
Given a couple, print out a given chromosome in that couple. 
"""
function print_of_pair(couple, chrom, spacing; semi=true) 
	show(chrom)
	print(" "^(spacing+namelength(couple) - namelength(chrom)))
	semi && print(";")
	print(" "^spacing)
end



"""
given genetic info, is it wildtype
"""
function wildtype(::Type{Chromosome}, n)
	return Chromosome(n, (Allele("+"),))
end
function wildtype(::Type{Couple}, n)
	return Couple(Chromosome(n, (Allele("+"),)))
end
function wildtype(::Type{Genotype})

	Genotype(wildtype(Couple, 1), wildtype(Couple, 2), wildtype(Couple, 3), wildtype(Couple, 4))
end



function origin()

	Parents(
		Genotype(
			Chromosome(1, (Allele(""))),
			Chromosome(2, (Allele(""))),
			Chromosome(3, (Allele(""))),
			Chromosome(4, (Allele("")))),

		Genotype(
			Chromosome(1, (Allele(""))),
			Chromosome(2, (Allele(""))),
			Chromosome(3, (Allele(""))),
			Chromosome(4, (Allele("")))))
end

"""
Empty Genotype 
"""
function origin()

	Parents(
		Genotype(
			Chromosome(1, (Allele("", not_lethal, recessive))),
			Chromosome(2, (Allele("", not_lethal, recessive))),
			Chromosome(3, (Allele("", not_lethal, recessive))),
			Chromosome(4, (Allele("", not_lethal, recessive)))),


		Genotype(
			Chromosome(1, (Allele("", not_lethal, recessive))),
			Chromosome(2, (Allele("", not_lethal, recessive))),
			Chromosome(3, (Allele("", not_lethal, recessive))),
			Chromosome(4, (Allele("", not_lethal, recessive)))))
end


"""
Chceks if genotype is origin
"""
function is_origin(g::Genotype)
	isempty(rand(g.p1.c1.genes).name)
end


"""
return set of chromosomes in genetic type
"""
chromosomes(p::Couple) = Set([p.c1, p.c2])
chromosomes(g::Genotype) = Set([g.p1.c1, g.p1.c2, g.p2.c1, g.p2.c2, g.p3.c1, g.p3.c2, g.p4.c1, g.p4.c2])

function children(gs::Vector{Family})

	out = Vector{Genotype}(undef, length(gs))
	for (i, g) in enumerate(gs)
		out[i] = g.child
	end
	return out
end

"""
is pair of chromosomes homozygous.
"""
homozygous(p1::Couple) = p1.c1 == p1.c2

