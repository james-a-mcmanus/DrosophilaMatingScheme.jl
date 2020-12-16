struct Allele
	name::String
end

struct CValue{N}
end

struct Chromosome{N,i}
	chromosome::CValue{N}	
	genes::NTuple{i, Allele}
end

struct Couple{N}
	c1::Chromosome{N}
	c2::Chromosome{N}
end

struct Genotype
	p1::Couple{1}
	p2::Couple{2}
	p3::Couple{3}
	p4::Couple{4}
end

struct Parents
	mum::Genotype
	dad::Genotype
end

struct Family
	child::Genotype
	parents::Parents
end


CValue(x) = CValue{x}()

Chromosome(chromosome_number::Int, alleles::NTuple{i, Allele}) where {i} = Chromosome(CValue(chromosome_number), alleles)
Chromosome(chromosome_number::Int, alleles::Vector{Allele}) = Chromosome(CValue(chromosome_number), tuple(alleles...))
Chromosome(chromosome_number::Int, alleles::Allele) where {i} = Chromosome(CValue(chromosome_number), (alleles,))
Chromosome(chromosome_number::Int, s::AbstractString) = parse(Chromosome, chromosome_number, s)

Couple(c1::Chromosome{N}) where {N} = Couple(c1, c1);
Couple(c2::Couple) = c2
Couple(i::Int, s::AbstractString) = parse(Couple, i, s)

Genotype(a,b,c,d) = Genotype(Couple(a), Couple(b), Couple(c), Couple(d))
Genotype(s::AbstractString) = parse(Genotype, s)