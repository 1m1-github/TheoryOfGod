module TOGIntelligenceLocal

using HTTP, JSON3

const URL = "http://localhost:8888/v1/messages"

function intelligence(
    model,
    input_system,
    input_user,
    max_output_tokens,
    temperature,
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
    
    response = HTTP.post(url, headers, body_string)
    response_body = String(response.body)
    result = JSON3.parse(response_body)
    result["choices"][1]["message"]["content"], ΔEnery(result, model)
end

const MAX_USD_IN_TICKS = 25 * 10^10
ΔEnery(result, model) = result["usage"]["cost_in_usd_ticks"] / MAX_USD_IN_TICKS

end