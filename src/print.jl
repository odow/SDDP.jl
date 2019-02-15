# ========================== Begin MIT Licensed Code ========================= #
# The following a modified version of that found at
#
# Humanize.jl    https://github.com/IainNZ/Humanize.jl
# as it was at commit be2c55008b501e17ed13c0a9aa791d40214385ea
#
# Copyright (c) 2017 Oscar Dowson, Iain Dunning, Julian Gehring
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software. # THE SOFTWARE IS PROVIDED "AS
# IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
# LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE
# AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
# CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

function humanize(value::Integer)
    return Printf.@sprintf("% 9d", value)
end

function humanize(value::Real)
    SUFFIX = [" ", "K", "M", "G", "T", "P", "E", "Z", "Y"]
    BASE = 1000.0
    bytes = abs(float(value)) # O.D. abs value
    unit = BASE
    suffix = SUFFIX[1]
    for (i, s) in enumerate(SUFFIX)
        unit = BASE ^ i
        suffix = s
        if bytes < unit
            break
        end
    end
    normalized_value = sign(value) * BASE * bytes / unit
    return Printf.@sprintf("% 8.3f%s", normalized_value, suffix)
end
# =========================== End MIT Licensed Code ========================== #

# ======================== Begin MPL2.0 Licensed Code ======================== #
#  Copyright 2018, Oscar Dowson. This Source Code Form is subject to the terms
#  of the Mozilla Public License, v.2.0. If a copy of the MPL was not
#  distributed with this file, You can obtain one at
#  http://mozilla.org/MPL/2.0/.

function print_banner(io=stdout)
    println(io, "———————————————————————————————————————————————————————————————————————————————")
    println(io, "                        Kokako - © Oscar Dowson, 2018-19.                      ")
    println(io, "———————————————————————————————————————————————————————————————————————————————")
    println(io, " Iteration | Simulation |      Bound |   Time (s)")
    println(io, "———————————————————————————————————————————————————————————————————————————————")
end

function print_iteration(log::Log, io=stdout)
    line = string(" ", humanize(log.iteration), " |  ",
        humanize(log.simulation_value), " |  ", humanize(log.bound), " |  ",
        humanize(log.time))
    println(io, line)
end
