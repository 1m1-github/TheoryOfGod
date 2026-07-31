module TOGIntelligenceLocal

using HTTP, JSON3
using TOGLogging: LOGS

const URL = "http://127.0.0.1:8080/v1/messages"

intelligence(complexity) = (a, b, c, d) -> intelligence(a, b, c, d, "DavidAU/Qwen3.6-27B-Fable-Fusion-711-Uncensored-Heretic-NM-DAU-NEO-MAX-MTP-GGUF:Q4_K_M")

function intelligence(
    input_system,
    input_user,
    max_output_tokens,
    temperature,
    model,
)
    headers = [
        "Content-Type" => "application/json",
    ]
    messages = [
        Dict("role" => "system", "content" => input_system),
        Dict("role" => "user", "content" => input_user),
    ]
    body = Dict(
        "model" => model,
        "messages" => messages,
        "temperature" => temperature,
        "max_tokens" => max_output_tokens,
    )
    body_string = JSON3.write(body)
    
    response = HTTP.post(URL, headers, body_string)
    response_body = String(response.body)
    
    # DEBUG
    ts = time()
    write(joinpath(LOGS, "latest-response.txt"), response_body)
    write(joinpath(LOGS, "$ts-response.txt"), response_body)
    # @info result
    # DEBUG

    result = JSON3.parse(response_body)
    result[:content][2][:text], 0.0
end

# const MAX_USD_IN_TICKS = 1000 * 10^10
# ΔEnery(result, model) = result["usage"]["cost_in_usd_ticks"] / MAX_USD_IN_TICKS

end