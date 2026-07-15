module TOGImagine

using FileIO, ColorTypes
using LoopOS: Peripheral
using TOGXAI: generateimage
import Base: take!

struct ImagineImage <: Peripheral end
# const IMAGINEIMAGE ?
take!(::Type{ImagineImage}, prompt::String) = generateimage(prompt=prompt)

# struct ImagineVideo <: Peripheral end
# function take!(::Type{ImagineVideo}, prompt::String)
#     bytes = generatevideo(prompt=prompt)
#     img = load(IOBuffer(bytes))
#     RGBA{N0f8}.(img)
# end

end
