## Task 2

### 1.

It is useful to exclude the input data from the source code because then one can easily adapt the source code to another example by simply changing the input data. In addition, the code becomes more clean since there are only read in statements instead of large declarations of input vectors and matrices.


### 2.

The function "dictzip" creates a dictionary using a dataframe and a pair value as arguments. With this function it is easy to create a dictionary which can be then used to look up certain values like the marginal costs of the power plants. Without the function one would always need to index the corresponding input dataframe. Thus the code would be more complicated.

As all timeseries include values for all given zones it would be handy to create a dictionary to look up for example the load timeseries for a specific zone. This is done by using the function "coldict" which takes a dataframe as an argument.

### 3.

A limited net capacity splits the market into zones. Otherwise we would have one market with one merit order curve since all electricity could be distributed through the lines.


### 4.

In order to implement a case with three zones one would just need to extend all input data by simply adding all relevant information for zone 3 which are:

1. If there is a new technology one would need to add this technology to the input data set technologies.
2. In the input data set zones one would need to declare which technologies are available in zone 3 along with the information about the generation or the storage capacity.
3. Adding a third zone means that there are two more lines connecting zone 3 with zone 1 and 2. Therefore these lines need to be added to the input data set lines along with the information about the net capacity.
4. As a last step, one would have to determine the timeseries for heat, load, pv and wind in zone 3 and add these information to the corresponding input data sets.
