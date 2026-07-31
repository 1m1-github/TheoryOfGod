"""
http://localhost:yourport
"""
module TOGAudioAnalogToDigitalBrowser

export JSLISTEN, JSIGNORE

# using HTTP, JSON3
using LoopOS
using LoopOS: Peripheral, listen
using TOGBroadcastBrowser: BroadcastBrowser
# using TOGXAI
# using TOGXAI: sttws
using TOGSTT: transcribe
# import TOGState: state
import Base: take!, put!

struct AudioAnalogToDigitalBrowser <: Peripheral
  channel::Channel{String}
end
take!(::AudioAnalogToDigitalBrowser) = take!(AUDIOANALOGTODIGITALBROWSER.channel)
put!(::AudioAnalogToDigitalBrowser, message::String) = put!(AUDIOANALOGTODIGITALBROWSER.channel, message)
# state(::AudioAnalogToDigitalBrowser) = "" # todo does anything?
const AUDIOANALOGTODIGITALBROWSER = AudioAnalogToDigitalBrowser(Channel{String}())

# const WSTASK = Ref{Task}()
# const WS = Ref{HTTP.WebSockets.WebSocket}()

function awaken()
  # @info "TOGAudioAnalogToDigitalBrowser.jl, awaken"
  put!(BroadcastBrowser, JSLISTEN)
  # WSTASK[], WS[] = sttws(output=output, interim_results=false)
  listen(AUDIOANALOGTODIGITALBROWSER)
end
WSSTARTED = false
# function output(data)
#   # @info "TOGAudioAnalogToDigitalBrowser.jl, output", data
#   if data["type"] == "transcript.created"
#     global WSSTARTED
#     WSSTARTED = true
#     return
#   end
#   data["speech_final"] || return
#   # msg = Dict(
#   #   :type => data["type"],
#   #   :text => data["text"],
#   #   :is_final => data["is_final"],
#   #   :speech_final => data["speech_final"],
#   #   :end_of_turn_confidence => data["end_of_turn_confidence"],
#   # )
#   isempty(data["text"]) && return
#   global AUDIOANALOGTODIGITALBROWSER
#   # put!(AUDIOANALOGTODIGITALBROWSER, string(msg))
#   put!(AUDIOANALOGTODIGITALBROWSER, data["text"])
# end
# audio(msg) = begin
#   global WSSTARTED
#   WSSTARTED && HTTP.WebSockets.send(WS[], msg)
# end
audio(data) = put!(AUDIOANALOGTODIGITALBROWSER, transcribe(data))

"Run this to start listening to audio on the browser. On by default."
const JSLISTEN = raw"""
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
  window.WORKLETNODE.port.onmessage = (event) => {fetch(`${window.BASE}/audio``, {method: 'POST',body: event.data,headers: { 'Content-Type': 'application/octet-stream' }}).catch(() => {});};
})()
"""

"Run this on the browser to stop audio recording on the browser"
const JSIGNORE = raw"""
if(window.AUDIOSTREAM)window.AUDIOSTREAM.getTracks().forEach(t=>t.stop())
if(window.AUDIOCTX)window.AUDIOCTX.close()
"""

end
