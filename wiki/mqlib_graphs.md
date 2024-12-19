# MQLib graph structure

MQLib works using the two instances data type, one for MaxCut and one for QUBO instances.
The MaxCut one can be constructed via providing an edgelist:

```c++
// Constructor (provide graph by providing edge list through vector of tuples)
  MaxCutInstance(const std::vector<Instance::InstanceTuple>& edgeList,
                 int dimension);
```

The type InstanceTuple is defined as:
```c++
typedef std::pair<std::pair<int, int>, double> InstanceTuple;
```

Where the pair of integers are the nodes belonging to the edge and the double value is its weight.
**Node labels start at 1**.

The dimension parameter is the number of nodes, as given in the first line of mc-files.


## Transforming Networkit graphs to MQLib

This should bedoable in two steps:

1. Iterate over the edges and create the edgelist vector
2. With the vector and the number of nodes to create the MQLib instances by using the constructor
