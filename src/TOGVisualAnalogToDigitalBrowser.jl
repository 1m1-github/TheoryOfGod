module TOGVisualAnalogToDigitalBrowser

using JSON3, Base64, PNGFiles, ColorTypes, FixedPointNumbers
using LoopOS: Peripheral, listen
using TOGBroadcastBrowser: BroadcastBrowser
import Base: take!, put!

struct VisualAnalogToDigitalBrowser <: Peripheral
    channel::Channel{Matrix{ColorTypes.RGBA{FixedPointNumbers.N0f8}}}
end
const VISUALANALOGTODIGITALBROWSER = VisualAnalogToDigitalBrowser(Channel{Matrix{ColorTypes.RGBA{FixedPointNumbers.N0f8}}}())
put!(::VisualAnalogToDigitalBrowser, img) = put!(VISUALANALOGTODIGITALBROWSER.channel, img)
take!(::VisualAnalogToDigitalBrowser) = take!(VISUALANALOGTODIGITALBROWSER.channel)
function take!(::Type{VisualAnalogToDigitalBrowser})
    put!(BroadcastBrowser, JSLOOK)
    take!(VISUALANALOGTODIGITALBROWSER)
end
function webcam(bytes)
    json = JSON3.read(String(bytes))
    b64 = split(json.image, "base64,")[2]
    img = PNGFiles.load(IOBuffer(base64decode(b64)))
    put!(VISUALANALOGTODIGITALBROWSER, img)
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
