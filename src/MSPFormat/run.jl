#  Copyright (c) 2023, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using SDDP
import HiGHS

include("MSPFormat.jl")
import .MSPFormat

problem = "Simple_Hydrothermal"
# problem = "GasStorage"

# ==============================================================================

model = MSPFormat.read_from_file(problem)
SDDP.write_to_file(model, "$problem.sof.json"; validation_scenarios = 1)
model, _ = SDDP.read_from_file("$problem.sof.json")
JuMP.set_optimizer(model, HiGHS.Optimizer)
SDDP.train(model; iteration_limit = 10)

# ==============================================================================

MSPFormat.convert_to_stochoptformat(
    "$problem.problem.json",
    "$problem.lattice.json",
    "$problem.sof.json",
)
model, _ = SDDP.read_from_file("$problem.sof.json")
JuMP.set_optimizer(model, HiGHS.Optimizer)
SDDP.train(model; iteration_limit = 10)
