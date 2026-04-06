using XAM
import GenomicFeatures.eachoverlap
import BioGenerics.header

"""
Convenience struct that holds the BAM filename, a reader for the file,
and the current record of the reader (if reading record-by-record).
Also includes access to convenience functions `close!()` `fetch()`, `next!()`, and `getValidBX()`.
"""
mutable struct BamReader
    file::String
    reader::XAM.BAM.Reader{IOStream}
    record::BAM.Record
    contigs::Dict{String,Int}
end

function BamReader(infile::String)
    rdr = open(BAM.Reader, infile, index = infile * ".bai")
    d = Dict{String,Int}()
    for i in findall(header(rdr), "SQ")
        d[i["SN"]] = parse(Int64, i["LN"]) 
    end
    BamReader(infile, rdr, BAM.Record(), d)
end

# range(start=0, step=0.01, length=2^10)

"""
Mutating reader of the BAM file in a `BamReader` that will overwrite
the existing `.record` field with the next record in the BAM reader.
"""
function next!(bamObj::BamReader)
    empty!(bamObj.record)
    if eof(bamObj.reader)
        return 
    end
    read!(bamObj.reader, bamObj.record)
end

"""
Close the underlying BAM.Reader
"""
function close!(bamObj::BamReader)
    close(bamObj.reader)
end

"""
Wraps `GenomicFeatures.eachoverlap` for a `BamReader` object,
returning a iterator of `BAM.record` over the requested interval.
"""
function fetch(bamObj::BamReader, chrom::String, pos::UnitRange{T} where T<:Core.Real)::XAM.BAM.OverlapIterator{IOStream}
    return eachoverlap(bamObj.reader, chrom, pos)
end

"""
Return the BX barcode of a `BAM.Record` if it has a BX tag and valid VX tag,
otherwise return `nothing`. Will still return `nothing` if there is a BX tag
without a VX tag.
"""
function getValidBX(rec::BAM.Record)::Union{String, Nothing}
    if haskey(rec, "VX")
        if rec["VX"]::UInt8 == 0x01
            if haskey(rec, "BX")
                return rec["BX"]::String
            end
        end
    end
    return nothing
end
