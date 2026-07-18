"""
http://localhost:yourport
"""
module TOGAudioAnalogToDigitalBrowser

using HTTP
using LoopOS
using LoopOS: Peripheral, listen
using TOGBroadcastBrowser: BroadcastBrowser
using TOGXAI: sttws
import Base: take!, put!

struct AudioAnalogToDigitalBrowser <: Peripheral
  channel::Channel{String}
end
take!(::AudioAnalogToDigitalBrowser) = take!(AUDIOANALOGTODIGITALBROWSER.channel)
put!(::AudioAnalogToDigitalBrowser, message::String) = put!(AUDIOANALOGTODIGITALBROWSER.channel, message)
const AUDIOANALOGTODIGITALBROWSER = AudioAnalogToDigitalBrowser(Channel{String}())

const WSTASK = Ref{Task}()
const WS = Ref{HTTP.WebSockets.WebSocket}()

function awaken()
  @info "TOGAudioAnalogToDigitalBrowser.jl, awaken"
  WSTASK[], WS[] = sttws(output=output)
  put!(BroadcastBrowser, JSLISTEN)
  listen(AUDIOANALOGTODIGITALBROWSER)
end

function output(data)
  @info "TOGAudioAnalogToDigitalBrowser.jl, output"
  msg = string(Dict(
    :type => data["type"],
    :text => data["text"],
    :is_final => data["is_final"],
    :speech_final => data["speech_final"],
    :end_of_turn_confidence => data["end_of_turn_confidence"],
  ))
  put!(AUDIOANALOGTODIGITALBROWSER, msg)
end

ws(msg) = begin
  @info "TOGAudioAnalogToDigitalBrowser.jl, ws"
  HTTP.WebSockets.send(WS[], msg)
end

const JSLISTEN = """
(async()=>{
window.AUDIOSTREAM=await navigator.mediaDevices.getUserMedia({audio:true})
window.AUDIOCTX=new AudioContext({sampleRate:16000})
window.AUDIOSOURCE=window.AUDIOCTX.createMediaStreamSource(window.AUDIOSTREAM)
window.WS.binaryType='arraybuffer'
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
  window.WORKLETNODE.port.onmessage = (event) => {if (window.WS.readyState === 1) window.WS.send(event.data)}
})()
"""
const JSIGNORE = """
if(window.AUDIOSTREAM)window.AUDIOSTREAM.getTracks().forEach(t=>t.stop())
if(window.AUDIOCTX)window.AUDIOCTX.close()
"""

end
