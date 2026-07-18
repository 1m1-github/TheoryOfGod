module TOGTextToAudioBrowser

export put!

using HTTP
using LoopOS: Peripheral
using TOGBroadcastBrowser: BroadcastBrowser
using TOGXAI: tts
using Base64
import Base: put!

struct TextToAudioBrowser <: Peripheral end

"""
Speak audio using the browser
"""
function put!(::Type{TextToAudioBrowser}, message::String)
    @info "TOGTextToAudioBrowser.jl, put!"
    audio = tts(message)
    mime = "audio/mpeg"
    b64 = base64encode(audio)
    js = """(()=>{const u=Uint8Array.from(atob('$b64'),c=>c.charCodeAt(0));const o=URL.createObjectURL(new Blob([u],{type:'$mime'}));const a=new Audio(o);a.play().catch(console.error)})()"""
    put!(BroadcastBrowser, js)
end

end
