#  Copyright 2017-19, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

# In this file we install some unregistered dependencies. This is a hack until
# such time that Julia gets its act together.

# For some reason we need to add this here.
push!(LOAD_PATH, "@stdlib")

const UNREGISTERED_URLS = [
    "https://github.com/odow/MathOptFormat.jl"
]

import Pkg
for url in UNREGISTERED_URLS
    Pkg.add(Pkg.PackageSpec(url=url))
end
