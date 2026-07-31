module TOGVisualAnalogToDigitalBrowser

export take!

using JSON3, Base64, PNGFiles, ColorTypes, FixedPointNumbers
using LoopOS: Peripheral, listen
using TOGBroadcastBrowser: BroadcastBrowser
using TOGXAI: imageb64totext
import Base: take!, put!

struct VisualAnalogToDigitalBrowser <: Peripheral
    channel::Channel{String}
end
const VISUALANALOGTODIGITALBROWSER = VisualAnalogToDigitalBrowser(Channel{String}())
put!(::VisualAnalogToDigitalBrowser, description::String) = put!(VISUALANALOGTODIGITALBROWSER.channel, description)
"""
See the real world via the browser webcam using.
example: `whatisee = take!(TOGVisualAnalogToDigitalBrowser.VisualAnalogToDigitalBrowser)`
"""
take!(::Type{VisualAnalogToDigitalBrowser}) = _take!(VISUALANALOGTODIGITALBROWSER)
function _take!(::VisualAnalogToDigitalBrowser)
  # @info "TOGVisualAnalogToDigitalBrowser.jl, take!"
  put!(BroadcastBrowser, JSLOOK)
  take!(VISUALANALOGTODIGITALBROWSER.channel)
end
function webcam(bytes)
  # @info "TOGVisualAnalogToDigitalBrowser.jl, webcam"
  json = JSON3.read(String(bytes))
  b64 = split(json.image, "base64,")[2]
  description = imageb64totext(b64=b64)
  put!(VISUALANALOGTODIGITALBROWSER, description)
end

const JSLOOK = raw"""
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
        fetch('${window.BASE}/webcam', {
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
