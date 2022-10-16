#=
After downloading and processing 1000 Genomes genotype array data with `01-download-and-prep.sh`, 
this Julia code can be used to prepare LD reference panels by ancestry groups as follows: 
- keep SNPs with minor allele frequency (MAF) > 0.05
- keep SNPs with polymorphic genotypes (i.e. remove SNPs that are heterozygous for all samples)
=#

dir = "/u/project/gandalm/shared/refGenomes/1kg"
cd(dir)

using Pkg
Pkg.activate(".")
Pkg.add(["SnpArrays", "CSV", "DataFrames", "Glob"])
using Downloads, SnpArrays, CSV, DataFrames, Glob

isdir("data") || mkdir("data")

@info "Download sample metadata"
url = joinpath("https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a", 
    "sample_info/integrated_call_samples_v3.20130502.ALL.panel")
meta = basename(url)
isfile("data/$(meta)") || Downloads.download(url, "data/$(meta)")

@info "Load sample metadata"
meta = CSV.read("data/$(meta)", DataFrame, header = false, skipto = 2)
rename!(meta, [:sample, :pop, :super, :gender])
ancestries = unique(meta.super)
for ancestry in ancestries
    df = filter(row -> row.super == ancestry, meta)
    println(ancestry, ": ", size(df, 1))
end

@info "Load 1000 Genomes data"
kgp = SnpData("kgp")

for ancestry in ancestries
    @info "Processing $(ancestry) samples"
    isfile("data/kgp.$(ancestry).maf0.05.geno.bed") ? continue : nothing

    @info "Subset to $(ancestry) samples"
    sample = meta.sample[meta.super .== "$(ancestry)"]
    rowinds = findall(in(sample), kgp.person_info.iid)
    SnpArrays.filter(kgp, rowinds, trues(size(kgp)[2]); des = "data/kgp.$(ancestry)")

    kgp_ancestry = SnpData("data/kgp.$(ancestry)")
    @info "# of samples: $(size(kgp_ancestry)[1]), # of SNPs: $(size(kgp_ancestry)[2])"

    @info "Subset to SNPs w/ MAF > 0.05"
    colinds = SnpArrays.filter(kgp_ancestry.snparray; min_maf = 0.05)[2]
    SnpArrays.filter(kgp_ancestry, trues(size(kgp_ancestry)[1]), colinds; des = "data/kgp.$(ancestry).maf0.05")   
    kgp_ancestry = SnpData("data/kgp.$(ancestry).maf0.05")
    @info "# of samples: $(size(kgp_ancestry)[1]), # of SNPs: $(size(kgp_ancestry)[2])"

    @info "Subset to SNPs w/ polymorphic genotypes"
    c = counts(kgp_ancestry.snparray, dims = 1)
    colinds = Int[]
    for j in 1:size(kgp_ancestry)[2]
        size(kgp_ancestry)[1] in c[:, j] ? nothing : push!(colinds, j)
    end
    SnpArrays.filter(kgp_ancestry, trues(size(kgp_ancestry)[1]), colinds; des = "data/kgp.$(ancestry).maf0.05.geno")
    kgp_ancestry = SnpData("data/kgp.$(ancestry).maf0.05.geno")
    @info "# of samples: $(size(kgp_ancestry)[1]), # of SNPs: $(size(kgp_ancestry)[2])"
    println("")
end

for ancestry in ancestries
    for plink in ["bed", "bim", "fam"]
        rm("data/kgp.$(ancestry).$(plink)")
        rm("data/kgp.$(ancestry).maf0.05.$(plink)")
    end
end