module TOGTypst

export typst

using PNGFiles, Typst_jll
using TOGMatrix: MATRIX_CACHE, rgbamatrix

const TYPST_TEMPLATE(content) = """
#set page(width: auto, height: auto, margin: (top: 5pt, bottom: 5pt, left: 5pt, right: 5pt))
#set text(size: 10pt)
$content
"""

function typst_to_matrix(; typst_code, dpi)
    h = hash(typst_code)
    haskey(MATRIX_CACHE, h) && return MATRIX_CACHE[h]
    Typst_jll.typst() do exe
        cmd = pipeline(
            `$exe compile - --format png --ppi $(dpi) -`,
            stdin=IOBuffer(TYPST_TEMPLATE(typst_code))
        )
        png_bytes = read(cmd)
        IOBuffer(png_bytes) |> PNGFiles.load
    end
end

typst(; typst_code, dpi) = rgbamatrix(typst_to_matrix(typst_code=typst_code, dpi=dpi))

end
