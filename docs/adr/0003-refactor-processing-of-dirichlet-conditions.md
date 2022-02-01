# 3. Refactor processing of Dirichlet conditions

Date: 2022-02-01

## Status

Accepted

## Context

Until now Dirichlet boundary conditions were handled at the same time
as other degrees of freedom introducing a lot of `if` switches
and increasing code obfuscation.

Moreover considering the coupling between several numerical schemes
would progressively lead to more and more complex code with many additional
conditionals (`if`, `then`, `else` nested clauses).


## Decision

The processing of boundary conditions in the assembly of the Jacobian
are dissociated from the processing of other degrees of freedom.
Boundary nodes are first handled as ordinary nodes.
Then the rows corresponding to Dirichlet nodes are enforced to
the equations solved for boundary conditions.

## Consequences

The main gain is that the code is clearer and it will be more
straightforward to implement new boundary conditions.
The `if` clauses spared a bit of computing work and memory,
but the overhead associated with the new algorithm
is at the most very light (to be confirmed on more tests)
and we might even gain time in some configurations with few
Dirichlet boundary nodes.
