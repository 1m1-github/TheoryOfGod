module TOGXAI

using HTTP, JSON3, Base64, FileIO, ImageIO, Images, ColorTypes
using TOGLogging: LOGS

const URL = "api.x.ai/v1"
const PREVIOUS_RESPONSE_ID = Ref{String}("")
const MAX_USD_IN_TICKS = 25 * 10^10

intelligence(complexity::Number) =
    if complexity ≤ 0.5
        (a, b, c, d) -> intelligence(a, b, c, d, "grok-build-0.1")
    else
        (a, b, c, d) -> intelligence(a, b, c, d, "grok-4.5")
    end

function intelligence(
    input_system,
    input_user,
    max_output_tokens,
    temperature,
    model,
)
    @info "TOGXAI.jl, intelligence"
    headers = [
        "Authorization" => """Bearer $(ENV["XAI_API_KEY"])""",
        "Content-Type" => "application/json",
    ]
    input = [
        Dict("role" => "system", "content" => input_system),
        Dict("role" => "user", "content" => input_user),
    ]
    body = Dict(
        "model" => model,
        "input" => input,
        # "reasoning" => Dict("effort" => "high"),
        "max_output_tokens" => max_output_tokens,
    )
    if !isempty(PREVIOUS_RESPONSE_ID[])
        body["previous_response_id"] = PREVIOUS_RESPONSE_ID[]
    end

    response = HTTP.post("https://" * URL * "/responses", headers, JSON3.write(body))
    response_body = String(response.body)
    result = JSON3.parse(response_body)

    # DEBUG
    ts = time()
    write(joinpath(LOGS, "latest-response.txt"), response_body)
    write(joinpath(LOGS, "$ts-response.txt"), response_body)
    # @info result
    # DEBUG

    PREVIOUS_RESPONSE_ID[] = result["id"]
    result["output"][2]["content"][1]["text"], ΔEnery(result, model)
end

ΔEnery(result, model) = result["usage"]["cost_in_usd_ticks"] / MAX_USD_IN_TICKS # depent on modelintelligence

function imageb64totext(;b64,mime="image/png",detail="high")
    data_uri = "data:$mime;base64,$b64"
    headers = [
        "Authorization" => """Bearer $(ENV["XAI_API_KEY"])""",
        "Content-Type" => "application/json",
    ]
    image = Dict(
        "type" => "input_image",
        "image_url" => data_uri,
        "detail" => detail,
    )
    input = [
        Dict("role" => "system", "content" => "This is what you see. Describe it for yourself."),
        Dict("role" => "user", "content" => [image])
    ]
    body = Dict(
        "model" => "grok-4.5",
        "input" => input,
    )

    response = HTTP.post("https://" * URL * "/responses", headers, JSON3.write(body))
    response_body = String(response.body)
    result = JSON3.parse(response_body)

    result["output"][2]["content"][1]["text"]
end

# """
# """
# stt(mp3filepath) = string(open(mp3filepath, "r") do io
#     headers = ["Authorization" => """Bearer $(ENV["XAI_API_KEY"])"""]
#     form = HTTP.Form([
#         "language" => "en",
#         "format" => "true",
#         "diarize" => true,
#         "file" => HTTP.Multipart(mp3filepath, io, "audio/mp3")
#     ])
#     response = HTTP.post(URL * "/stt", headers, form)
#     response_body = String(response.body)
#     result = JSON3.parse(response_body)
#     result["text"]
# end)

function sttws(;
    output,
    interim_results=true,
    diarize=true,
    filler_words=true,
    smart_turn=0.8,
    smart_turn_timeout=3000,
    language="en"
)
    # @info "TOGXAI.jl, sttws"
    url = "wss://" * URL * "/stt?interim_results=$(interim_results)&diarize=$(diarize)&filler_words=$(filler_words)&smart_turn=$(smart_turn)&smart_turn_timeout=$(smart_turn_timeout)&language=$(language)"
    headers = HTTP.Headers()
    push!(headers, "Authorization" => """Bearer $(ENV["XAI_API_KEY"])""")
    wss = HTTP.WebSockets.open(url, headers=headers)
    task = @async for msg in wss output(JSON3.read(msg)) end
    task, wss
end

"""
Use inline speech tags per X AI Speech API as needed.
Examples:
inline tags: [pause], [long-pause], [hum-tune], [laugh], [chuckle], [giggle], [cry], [tsk], [tongue-click], [lip-smack], [breath], [inhale], [exhale], [sigh]
wrapping tags: <soft>, <whisper>, <loud>, <build-intensity>, <decrease-intensity>, <higher-pitch>, <lower-pitch>, <slow>, <fast>, <sing-song>, <singing>, <laugh-speak>, <emphasis>
https://docs.x.ai/developers/model-capabilities/audio/text-to-speech#speech-tags
""" # todo maybe add all?
function tts(text)
    @info "TOGXAI.jl, tts"
    headers = [
        "Authorization" => """Bearer $(ENV["XAI_API_KEY"])""",
        "Content-Type" => "application/json",
    ]
    body = Dict(
        "text" => text,
        "language" => "en",
        "voice_id" => "ara",
        "text_normalization" => true,
        "output_format" => Dict(
            "codec" => "mp3",
            "sample_rate" => 48000,
            "bit_rate" => 32000,
        )
    )
    response = HTTP.post("https://" * URL * "/tts", headers, JSON3.write(body))
    response.body
end

function generateimage(; prompt, aspect_ratio="auto", resolution="1k")
    @info "TOGXAI.jl, generateimage"
    headers = [
        "Authorization" => """Bearer $(ENV["XAI_API_KEY"])""",
        "Content-Type" => "application/json",
    ]
    body = Dict(
        "prompt" => prompt,
        "model" => "grok-imagine-image-quality",
        "aspect_ratio" => aspect_ratio,
        "resolution" => resolution,
        "response_format" => "b64_json",
    )
    response = HTTP.post("https://" * URL * "/images/generations", headers, JSON3.write(body))
    json = JSON3.parse(String(response.body))
    b64 = json[:data][1][:b64_json]
    bytes = base64decode(b64)
    load(IOBuffer(bytes))
end

# function generatevideo(;prompt, aspect_ratio="auto", resolution="1k")
#     headers = [
#         "Authorization" => """Bearer $(ENV["XAI_API_KEY"])""",
#         "Content-Type" => "application/json",
#     ]
#     body = Dict(
#         "prompt" => prompt,
#         "model" => "grok-imagine-image-quality",
#         "aspect_ratio" => aspect_ratio,
#         "resolution" => resolution,
#         # "response_format" => "b64_json",
#     )
#     response = HTTP.post("https://" * URL * "/images/generations", headers, JSON3.write(body))
#     response.body
# end

end
