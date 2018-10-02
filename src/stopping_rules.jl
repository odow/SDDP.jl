abstract type AbstractStoppingRules end

struct Statistical end
stopping_rule_status(::Statistical) = :statistical

struct BoundStalling end
stopping_rule_status(::BoundStalling) = :bound_stalling

function convergence_test(graph::PolicyGraph,
                          stopping_rules::Vector{AbstractStoppingRules})
    for stopping_rule in stopping_rules
        if check_stopping_rule(graph)
            return true, stopping_rule_status(stopping_rule)
        end
    end
    return false, :not_solved
end
