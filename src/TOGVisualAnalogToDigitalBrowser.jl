module TOGVisualAnalogToDigitalBrowser

export take!

using Serialization
using JSON3, Base64, PNGFiles, ColorTypes, FixedPointNumbers
using LoopOS: Peripheral, listen
using TOGBroadcastBrowser: BroadcastBrowser
import Base: take!, put!

struct VisualAnalogToDigitalBrowser <: Peripheral
    # channel::Channel{Matrix{ColorTypes.RGBA{FixedPointNumbers.N0f8}}}
    channel::Channel{String}
end
# const VISUALANALOGTODIGITALBROWSER = VisualAnalogToDigitalBrowser(Channel{Matrix{ColorTypes.RGBA{FixedPointNumbers.N0f8}}}())
const VISUALANALOGTODIGITALBROWSER = VisualAnalogToDigitalBrowser(Channel{String}())
# put!(::VisualAnalogToDigitalBrowser, img::Matrix{ColorTypes.RGBA{FixedPointNumbers.N0f8}}) = begin
put!(::VisualAnalogToDigitalBrowser, description::String) = begin
  @info "TOGVisualAnalogToDigitalBrowser.jl, put!"
  put!(VISUALANALOGTODIGITALBROWSER.channel, description)
end
"""
See the real world via the browser webcam using.
example: `whatisee = take!(TOGVisualAnalogToDigitalBrowser.VisualAnalogToDigitalBrowser)`
"""
take!(::Type{VisualAnalogToDigitalBrowser}) = _take!(VISUALANALOGTODIGITALBROWSER)
# take!(::VisualAnalogToDigitalBrowser) = take!(VISUALANALOGTODIGITALBROWSER.channel)
# function take!(::Type{VisualAnalogToDigitalBrowser})
function _take!(::VisualAnalogToDigitalBrowser)
  @info "TOGVisualAnalogToDigitalBrowser.jl, take!"
  put!(BroadcastBrowser, JSLOOK)
  take!(VISUALANALOGTODIGITALBROWSER.channel)
end
function webcam(bytes)
  @info "TOGVisualAnalogToDigitalBrowser.jl, webcam"
  # serialize("bytes", bytes)
  json = JSON3.read(String(bytes))
  b64 = split(json.image, "base64,")[2]
  description = imagetotext(b64=b64)
  # img = PNGFiles.load(IOBuffer(base64decode(b64)))
  # PNGFiles.save("img.png",img)
  # put!(VISUALANALOGTODIGITALBROWSER, img)
  put!(VISUALANALOGTODIGITALBROWSER, description)
end

const JSLOOK = """
(() => {
  navigator.mediaDevices.getUserMedia({ video: true })
    .then(stream => {
      const video = document.createElement('video')
      video.srcObject = stream
      video.onloadedmetadata = () => {
        video.play()
        const canvas = document.createElement('canvas')
        canvas.width = video.videoWidth
        canvas.height = video.videoHeight
        const ctx = canvas.getContext('2d')
        ctx.drawImage(video, 0, 0)
        const dataUrl = canvas.toDataURL('image/png')
        stream.getTracks().forEach(t => t.stop())
        fetch('/webcam', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ image: dataUrl })
        }).then(r => r.text()).then(console.log).catch(console.error)
      }
    })
    .catch(console.error)
})()
"""

end
