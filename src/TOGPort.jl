module TOGPort

using Sockets

function openport(hint=8888)
    port, server = listenany(hint)
    close(server)
    Int(port)
end

end
