# Copyright (c) 2017-26: Oscar Dowson and SDDP.jl contributors.
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v2.0. If a copy of the MPL was not distributed with this file, You can obtain
# one at http://mozilla.org/MPL/2.0/.

function launch_websocket(server_to_client::Channel{Log})
    return HTTP.serve!("127.0.0.1", 8000) do req::HTTP.Request
        return HTTP.sse_stream(200) do stream
            for log in server_to_client
                data = Dict(
                    "iteration" => log.iteration,
                    "bound" => log.bound,
                    "simulation" => log.simulation_value,
                    "time" => log.time,
                    "solves" => log.total_solves,
                )
                event = HTTP.SSEEvent(JSON.json(data); event = "iteration")
                write(stream, event)
            end
        end
    end
end

function launch_dashboard()
    server_to_client = Channel{Log}(typemax(Int))
    server = launch_websocket(server_to_client)
    launch_file(joinpath(@__DIR__, "dashboard.html"))
    function dashboard_callback(log::Union{Log,Nothing}, close_flag::Bool)
        if close_flag
            HTTP.forceclose(server)
        elseif log !== nothing
            put!(server_to_client, log)
        end
        return
    end
    return dashboard_callback
end
