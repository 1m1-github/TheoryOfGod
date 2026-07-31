module TOGSTT

using Whisper, Suppressor

const SILENCE_THRESHOLD = 1e-6
const WHISPER_FILENAME = "/Users/1m1/Downloads/aos/transcription/ggml-large-v3.bin"
# todo init on awaken only
const WHISPER_CONTEXT = @suppress Whisper.whisper_init_from_file(WHISPER_FILENAME)
const WHISPER_PARAMS = @suppress Whisper.whisper_full_default_params(Whisper.LibWhisper.WHISPER_SAMPLING_GREEDY)

function transcribe(data)
    @info "transcribe", typeof(data)
    @suppress begin
        result = ""
        all(abs.(data) .< SILENCE_THRESHOLD) && return result
        Whisper.whisper_full_parallel(WHISPER_CONTEXT, WHISPER_PARAMS, data, length(data), 1)
        n_segments = Whisper.whisper_full_n_segments(WHISPER_CONTEXT)
        for i in 0:n_segments-1
            txt = Whisper.whisper_full_get_segment_text(WHISPER_CONTEXT, i)
            result *= unsafe_string(txt)
        end
        @info "transcribe", result
        result
    end
end

end
