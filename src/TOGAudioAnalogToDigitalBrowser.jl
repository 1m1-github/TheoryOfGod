"""
http://localhost:yourport
"""
module TOGAudioAnalogToDigitalBrowser

using HTTP, JSON3
using LoopOS
using LoopOS: Peripheral, listen
using TOGBroadcastBrowser: BroadcastBrowser
using TOGXAI
using TOGXAI: sttws
import Base: take!, put!

struct AudioAnalogToDigitalBrowser <: Peripheral
  channel::Channel{String}
end
take!(::AudioAnalogToDigitalBrowser) = begin
  @info "TOGAudioAnalogToDigitalBrowser.jl, take!"
  take!(AUDIOANALOGTODIGITALBROWSER.channel)
end
put!(::AudioAnalogToDigitalBrowser, message::String) = begin
  @info "TOGAudioAnalogToDigitalBrowser.jl, put!", message
  put!(AUDIOANALOGTODIGITALBROWSER.channel, message)
end
const AUDIOANALOGTODIGITALBROWSER = AudioAnalogToDigitalBrowser(Channel{String}())

const WSTASK = Ref{Task}()
const WS = Ref{HTTP.WebSockets.WebSocket}()

function awaken()
  @info "TOGAudioAnalogToDigitalBrowser.jl, awaken"
  put!(BroadcastBrowser, JSLISTEN)

  url = "wss://" * TOGXAI.URL * "/stt?interim_results=true&diarizetrue&filler_words=true&smart_turn=0.8&smart_turn_timeout=3000&language=en"
  headers = HTTP.Headers()
  push!(headers, "Authorization" => """Bearer $(ENV["XAI_API_KEY"])""")
  WS[] = HTTP.WebSockets.open(url, headers=headers)
  @info "wss", WS[]
  WSTASK[] = @async for msg in WS[]
        # @info "WS[] msg", msg
        # try
            output(JSON3.read(msg))
        # catch e 
            # @info e
            # sprint(showerror, e, catch_backtrace())
        # end
    end

  # WSTASK[], WS[] = sttws(output=output)
  listen(AUDIOANALOGTODIGITALBROWSER)
end
WSSTARTED = false
function output(data)
  # @info "TOGAudioAnalogToDigitalBrowser.jl, output", data
  if data["type"] == "transcript.created"
    global WSSTARTED
    WSSTARTED = true
    return
  end
  msg = Dict(
    :type => data["type"],
    :text => data["text"],
    :is_final => data["is_final"],
    :speech_final => data["speech_final"],
    :end_of_turn_confidence => data["end_of_turn_confidence"],
  )
  isempty(msg[:text]) && return
  # haskey(data, "is_final") && push!(msg, data["is_final"])
  # haskey(data, "speech_final") && push!(msg, data["speech_final"])
  global AUDIOANALOGTODIGITALBROWSER
  put!(AUDIOANALOGTODIGITALBROWSER, string(msg))
end

ws(msg) = begin
  # @info "TOGAudioAnalogToDigitalBrowser.jl, ws", msg
  # @info "ws(msg)", WS, WSTASK, length(msg)
  global WSSTARTED
  WSSTARTED && HTTP.WebSockets.send(WS[], msg)
end

const JSLISTEN = """
(async()=>{
window.AUDIOSTREAM=await navigator.mediaDevices.getUserMedia({audio:true})
window.AUDIOCTX=new AudioContext({sampleRate:16000})
window.AUDIOSOURCE=window.AUDIOCTX.createMediaStreamSource(window.AUDIOSTREAM)
const workletCode = `
    class PCMProcessor extends AudioWorkletProcessor {
      process(inputs) {
        const input = inputs[0];
        if (input && 0 < input.length) {
          const data = input[0];
          const pcm = new Int16Array(data.length);
          for (let i = 0; i < data.length; i++) {
            pcm[i] = Math.max(-32768, Math.min(32767, data[i] * 32767 | 0));
          }
          this.port.postMessage(pcm.buffer, [pcm.buffer]);
        }
        return true;
      }
    }
    registerProcessor('pcm-processor', PCMProcessor);
  `
  const blob = new Blob([workletCode], { type: 'application/javascript' });
  await window.AUDIOCTX.audioWorklet.addModule(URL.createObjectURL(blob));
  window.WORKLETNODE = new AudioWorkletNode(window.AUDIOCTX, 'pcm-processor');
  window.AUDIOSOURCE.connect(window.WORKLETNODE);
  window.WORKLETNODE.port.onmessage = (event) => {fetch('/audio', {method: 'POST',body: event.data,headers: { 'Content-Type': 'application/octet-stream' }}).catch(() => {});};
})()
"""
# const JSLISTEN = """
# (async()=>{
# window.AUDIOSTREAM=await navigator.mediaDevices.getUserMedia({audio:true})
# window.AUDIOCTX=new AudioContext({sampleRate:16000})
# window.AUDIOSOURCE=window.AUDIOCTX.createMediaStreamSource(window.AUDIOSTREAM)
# window.WS.binaryType='arraybuffer'
# const workletCode = `
#     class PCMProcessor extends AudioWorkletProcessor {
#       process(inputs) {
#         const input = inputs[0];
#         if (input && 0 < input.length) {
#           const data = input[0];
#           const pcm = new Int16Array(data.length);
#           for (let i = 0; i < data.length; i++) {
#             pcm[i] = Math.max(-32768, Math.min(32767, data[i] * 32767 | 0));
#           }
#           this.port.postMessage(pcm.buffer, [pcm.buffer]);
#         }
#         return true;
#       }
#     }
#     registerProcessor('pcm-processor', PCMProcessor);
#   `
#   const blob = new Blob([workletCode], { type: 'application/javascript' });
#   await window.AUDIOCTX.audioWorklet.addModule(URL.createObjectURL(blob));
#   window.WORKLETNODE = new AudioWorkletNode(window.AUDIOCTX, 'pcm-processor');
#   window.AUDIOSOURCE.connect(window.WORKLETNODE);
#   window.WORKLETNODE.port.onmessage = (event) => {if (window.WS.readyState === 1) window.WS.send(event.data)}
# })()
# """
const JSIGNORE = """
if(window.AUDIOSTREAM)window.AUDIOSTREAM.getTracks().forEach(t=>t.stop())
if(window.AUDIOCTX)window.AUDIOCTX.close()
"""

end
