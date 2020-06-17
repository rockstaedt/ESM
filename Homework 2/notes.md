## Task 2

### 1.

It is useful to exclude the input data from the source code because then one can easily adapt the source code to another example by simply changing the input data. In addition, the code becomes more clean since there are only one read in statement instead of large declarations of input vectors.


### 2.

The function "dictzip" creates a dictionary using a dataframe and a pair value as arguments. With this function it is easy to create a dictionary which can be then used to look up certain values like the marginal costs of the power plants. Without the function one would always need to index the input dataframe. Thus the code would be more complicated.

As all timeseries include values for all given zones it would be handy to create a dictionary to look up for example the load timeseries for a specific zone. This is done by using the function "coldict" which takes a dataframe as an argument.

### 3.

A limited net capacity splits the market into zones. Otherwise we would have one market with one merit order curve since all electricity could be distributed.


### 4.

Extend all input data (technologies, zones, net transfer capacity, timeseries for heat, pv, wind and load) by adding all relevant information for zone 3:

- technologies
  - 
