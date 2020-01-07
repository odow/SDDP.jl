#  Copyright 2017-20, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

function launch_websocket(server_to_client::Channel{Log}, client_to_server::Channel{Bool})
    @async HTTP.WebSockets.listen("127.0.0.1", UInt16(8001)) do ws
        while !eof(ws)
            if isopen(server_to_client) && isopen(client_to_server)
                new_log = take!(server_to_client)
                write(
                    ws,
                    JSON.json(Dict(
                        "iteration" => new_log.iteration,
                        "bound" => new_log.bound,
                        "simulation" => new_log.simulation_value,
                        "time" => new_log.time,
                    )),
                )
                # Now, wait for Plotly to return a signal that it has updated
                # the plot.
                ready = String(readavailable(ws))
                put!(client_to_server, ready == "ready")
            else
                close(ws)
            end
        end
    end
end

function launch_dashboard()
    # We set up two channels:
    #  1) a channel from the server (Julia) to the client (JS)
    #  2) a channel from the client (JS) to the server (Julia)
    server_to_client = Channel{Log}(typemax(Int))
    client_to_server = Channel{Bool}(typemax(Int))
    # Launch the websocket, and wait for the task to start before launching the
    # HTML file to avoid the connection failing on the client-side.
    task = launch_websocket(server_to_client, client_to_server)
    while !istaskstarted(task)
        sleep(0.1)
    end
    launch_file(joinpath(@__DIR__, "dashboard.html"))
    # Return the plotting callback.
    return (log::Union{Log,Nothing}, close_flag::Bool) -> begin
        if close_flag
            # We've received the signal to close the websocket because SDDP has
            # terminated.
            close(server_to_client)
            close(client_to_server)
        elseif log !== nothing
            # Send the latest log to the client.
            put!(server_to_client, log)
            # Wait for the signal that Plotly has updated the plot.
            take!(client_to_server)
        end
        return
    end
end
