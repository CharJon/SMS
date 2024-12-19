# KL

## phase1OptimizationSubRoutine

Given a partition of the graph, this subroutine tries to improve it:

Improvement is done by computing Gain Pair, a pair of nodes in different partitions, which improve the solution value,
if swapped.
For this purpose the D-value of each node is computed, this value is the difference in internal cost vs external cost,
i.e. how much cost the nodes prevents from other nodes within its current partition, vs how much in currently incurs
from in neighbours in the other one.
We are looking at cost, as the KL heuristic tries to minimize the cost of a cut, by multiplying the weight with -1, we
get a max-cut heuristic.

Using the D-values, a pair of nodes, the Gain Pair, whose gain is maximized is selected. After selecting it, the
D-values are updated, to reflect the swap (however the swap itself is not yet executed). This process is repeated until
a number of Gain Pairs is found (usually as many as the number of nodes in the smaller partition).

Then, the index k is computed, for which the first k found Gain Pairs maximize the gain. This k pairs are then swapped.

As long the there is a k such that there is a positive total gain, the procedure is repeated, starting at calculating
the D-values again.

## repartitionOptimization

Using an already optimized partition (by the Phase 1 algorithm), create 2 subgraphs, one for each parition. Then
calculate a KL Partition on both graphs, using the Phase 1 algorithm.
This results in two new solution, each with a 0 and 1 partition. Recombine those into one new partition, and optimize
using Phase 1 again.

Currently implemented: (1,1) with (0,0) and (1,0) with (0,1). Where (subgraph_index, partition_side).

## phase1Optimization(McSolution &)

Use a given McSolution to start optimizing from, initial partition will be the McSolution. With no new nodes.

## phase1Optimization()

Put all nodes into one partition and fill the other one with n dummy nodes, then optimize. Dummy nodes have weight zero
to all other nodes.

# Markdown Math Examples

$`u \notin V`$

```math
u \in V
```
