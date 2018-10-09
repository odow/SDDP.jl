# The following a modified version of that found at
#
# Humanize.jl    https://github.com/IainNZ/Humanize.jl
# as it was at commit be2c55008b501e17ed13c0a9aa791d40214385ea
#
# Copyright (c) 2017 Oscar Dowson, Iain Dunning, Julian Gehring
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# as a future note, we do this because generated functions are currently not GC'd

# O.D. 2016 renamed
const suffix = [" ", "K", "M", "G", "T", "P", "E", "Z", "Y"]
# O.D. fix base
const base   = 1000.0

humanize51f(v,s) = @sprintf("% 5.1f%s", v, s)
humanize52f(v,s) = @sprintf("% 5.2f%s", v, s)
humanize83f(v,s) = @sprintf("% 8.3f%s", v, s)
humanize5d(v,s)  = @sprintf("% 5d%s", v, s)

function humanize(value::Number, fmt_str::String="5.1f")
    if fmt_str == "5.1f"
        return humanize(value, humanize51f)
    elseif fmt_str == "5.2f"
        return humanize(value, humanize52f)
    elseif fmt_str == "8.3f"
        return humanize(value, humanize83f)
    elseif fmt_str == "5d"
        return humanize(value, humanize5d)
    end
    error("Format string $fmt_str not intialised.")
end
# O.D. 2017 renamed. drop style option
function humanize(value::Number, fmt_str::Function)
    bytes   = abs(float(value)) # O.D. abs value
    unit    = base
    s       = suffix[1]
    for (i,s) in enumerate(suffix)
        unit = base ^ (i)
        bytes < unit && break
    end
    # O.D. add sign
    return fmt_str(sign(value)*base * bytes / unit, s)::String
end
