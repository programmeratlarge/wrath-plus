using XAM
import GenomicFeatures.eachoverlap

"""
Convenience struct that holds the BAM filename, a reader for the file,
and the current record of the reader (if reading record-by-record).
Also includes access to convenience functions `fetch()`, `next!()`, and `getValidBX()`.
"""
mutable struct BamReader
    file::String
    reader::XAM.BAM.Reader{IOStream}
    record::BAM.Record

    BamReader(infile) = new(infile, open(BAM.Reader, infile, index = infile * ".bai"), BAM.Record())
end


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
Wraps `GenomicFeatures.eachoverlap` for a `BamReader` object,
returning a iterator of `BAM.record` over the requested interval.
"""
function fetch(bamObj::BamReader, chrom::String, pos::UnitRange{T} where T<:Core.Real)::XAM.BAM.OverlapIterator{IOStream}
    return GenomicFeatures.eachoverlap(bamObj.reader, chrom, pos)
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
