module TOGAnthropic

using HTTP, JSON3

const URL = "https://api.anthropic.com/v1/messages"

function intelligence(
    model,
    input_system,
    input_user,
    max_output_tokens,
    temperature,
)
    headers = [
        "x-api-key" => ENV["ANTHROPIC_API_KEY"],
        "anthropic-version" => "2023-06-01",
        "Content-Type" => "application/json"
    ]
    system = [Dict("type" => "text", "text" => input_system, "cache_control" => Dict("type" => "ephemeral"))]
    messages = [Dict("role" => "user", "content" => input_user)]
    body = Dict(
        "model" => model,
        "system" => system,
        "messages" => messages,
        "temperature" => temperature,
        "max_tokens" => max_output_tokens,
    )
    body_string = JSON3.write(body)
    response = HTTP.post(url, headers, body_string)
    response_body = String(response.body)
    result = JSON3.parse(response_body)
    result["content"][1]["text"], ΔEnery(result, model)
end

const MAX_USD = 25
const MAX_CUMULATIVE_CACHED_READ_TOKENS = MAX_USD / 0.5 * 1e6
const MAX_CUMULATIVE_CACHED_WRITE_TOKENS = MAX_USD / 6.25 * 1e6
const MAX_CUMULATIVE_READ_TOKENS = MAX_USD / 5 * 1e6
const MAX_CUMULATIVE_WRITE_TOKENS = MAX_USD / 25 * 1e6
function ΔEnery(result, model)
    ΔE = result["usage"]["cache_read_input_tokens"] / MAX_CUMULATIVE_CACHED_READ_TOKENS
    ΔE += result["usage"]["cache_creation_input_tokens"] / MAX_CUMULATIVE_CACHED_WRITE_TOKENS
    # result["usage"]["ephemeral_5m_input_tokens"] / MAX_CUMULATIVE_CACHED_READ_BITS
    # result["usage"]["ephemeral_1h_input_tokens"] / MAX_CUMULATIVE_CACHED_READ_BITS
    ΔE += result["usage"]["input_tokens"] / MAX_CUMULATIVE_READ_TOKENS
    ΔE += result["usage"]["output_tokens"] / MAX_CUMULATIVE_WRITE_TOKENS
end

end
