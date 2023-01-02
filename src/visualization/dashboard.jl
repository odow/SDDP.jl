#  Copyright (c) 2017-23, Oscar Dowson and SDDP.jl contributors.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

function launch_websocket(server_to_client::Channel{Log})
    host, port = HTTP.ip"127.0.0.1", 8000
    server = HTTP.Sockets.listen(host, port)
    HTTP.WebSockets.listen!(host, port; server = server) do ws
        for msg in ws
            break  # Wait for the first message, then exit the loop
        end
        while true
            if !isready(server_to_client)
                sleep(0.01)
                continue
            end
            new_log = take!(server_to_client)
            data = Dict(
                "iteration" => new_log.iteration,
                "bound" => new_log.bound,
                "simulation" => new_log.simulation_value,
                "time" => new_log.time,
            )
            HTTP.send(ws, JSON.json(data))
        end
        return
    end
    return server
end

function launch_dashboard()
    server_to_client = Channel{Log}(typemax(Int))
    server = launch_websocket(server_to_client)
    launch_file(joinpath(@__DIR__, "dashboard.html"))
    return (log::Union{Log,Nothing}, close_flag::Bool) -> begin
        if close_flag
            close(server)
        elseif log !== nothing
            put!(server_to_client, log)
        end
        return
    end
end
