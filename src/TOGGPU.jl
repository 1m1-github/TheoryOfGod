module TOGGPU

using KernelAbstractions
# using Metal
# const GPU_BACKEND = MetalBackend()
const GPU_BACKEND = CPU()
const GPU_BACKEND_WORKGROUPSIZE = 2^8

end
