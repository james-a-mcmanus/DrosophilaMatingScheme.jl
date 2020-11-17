import Base: show, length, isequal, parse, hash, string

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
		alleles[i] = Allele(al) 
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
Base.string(g::Genotype) = string(g.p1) * " ; " * string(g.p2) * " ; " * string(g.p3) * " ; " * string(g.p4)


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

	return (Couple(p1.c1, p2.c1), Couple(p1.c1, p2.c2), Couple(p1.c2, p2.c1), Couple(p1.c2, p2.c2))
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


"""
Takes tuples of chromosomes, combines them into all possible genotypes.
"""
# This isnt' very efficient. If we insist on there being 4 posible pairs, then we can't optimise by e.g. reducing the pairs when they are homozygous etc. But it is type-stable.
function permute_chromosomes(cpl1::NTuple{4, Couple{1}}, cpl2::NTuple{4, Couple{2}}, cpl3::NTuple{4, Couple{3}}, cpl4::NTuple{4, Couple{4}})
	out = Vector{Genotype}(undef, 256)
	i = 1
	for pr1 in cpl1
		for pr2 in cpl2
			for pr3 in cpl3
				for pr4 in cpl4
					out[i] = Genotype(pr1, pr2, pr3, pr4)
					i += 1
				end
			end
		end
	end
	return out
end


function permute_chromosomes(d::AbstractDict{Genotype, Parents}, cpl1::NTuple{4, Couple{1}}, cpl2::NTuple{4, Couple{2}}, cpl3::NTuple{4, Couple{3}}, cpl4::NTuple{4, Couple{4}}, parents::Parents)

	for pr1 in cpl1
		for pr2 in cpl2
			for pr3 in cpl3
				for pr4 in cpl4
					
					g = Genotype(pr1, pr2, pr3, pr4)
					if !haskey(d, g) 
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
Base.show(io::IO, al::Allele) = print(al.name)
Base.show(io::IO, ::MIME"text/plain", al::Allele) = print(al.name)

Base.show(io::IO, ::MIME"text/plain", ch::Chromosome) = show(io, ch)
function Base.show(io::IO, ch::Chromosome)
	for i = 1:length(ch.genes)-1
		show(ch.genes[i])
		print(", ")
	end
	show(ch.genes[length(ch.genes)])
end

Base.show(io::IO, ::MIME"text/plain", cpl::Couple) = show(io, cpl)
function Base.show(io::IO, cpl::Couple)

	show(cpl.c1)
	print("\n")
	println(Char(8212)^max(namelength(cpl.c1),namelength(cpl.c2)))
	show(cpl.c2)
end

Base.show(io::IO, ::MIME"text/plain", gn::Genotype) = show(io, gn)
function Base.show(io::IO, gn::Genotype)

	spacing = 3

	print("\n\n")

	print_of_pair(gn.p1, gn.p1.c1, spacing)
	print_of_pair(gn.p2, gn.p2.c1, spacing)
	print_of_pair(gn.p3, gn.p3.c1, spacing)
	print_of_pair(gn.p4, gn.p4.c1, spacing, semi=false)

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

# Don't know why this has to be overloaded, otherwise it prints it 3 times???
Base.show(io::IO, ::MIME"text/plain", gn::Vector{Genotype}) = show(io, gn)
function Base.show(io::IO, gn::Vector{Genotype})

	for gt in gn
		show(gt)
	end
end



""" 
Equivalence Functions
"""
Base.isequal(a1::Allele, a2::Allele) = a1.name == a2.name
Base.isequal(chr1::Chromosome{N}, chr2::Chromosome{N}) where {N} = Set(chr1.genes) == Set(chr2.genes)
Base.isequal(p1::Couple{N}, p2::Couple{N}) where {N} = Set(p1.c1, p1.c2) == Set(p2.c1, p2.c2)


"""
Some helper functions.
"""
function example1()

	Genotype("sn/FM7a;sp/CyO;ser/TM3-sb;")
end

function example2()
	Genotype("w[+]; If/CyO ; MKRS/TM6b ;")
end

namelength(a::Allele) = length(a.name)
namelength(ch::Chromosome) = sum([namelength(g) + 2 for g in ch.genes]) - 2
namelength(pr::Couple) = max(namelength(pr.c1), namelength(pr.c2))
namelength(gn::Genotype) = namelength(gn.p1) + namelength(gn.p2) + namelength(gn.p3) + namelength(gn.p4)

function print_of_pair(couple, chrom, spacing; semi=true) 
	show(chrom)
	print(" "^(spacing+namelength(couple) - namelength(chrom)))
	semi && print(";")
	print(" "^spacing)
end

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