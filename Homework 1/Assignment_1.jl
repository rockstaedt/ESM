using DataFrames
using CSV

# creating a DataFrames

data = DataFrame()

# creating arrays for each column
technology = ["Nuclear","Lignite","HardCoal","NaturalGas"]
FuelCost   = [3,6.21,10.6,31.08]
efficiency = [0.33,0.42,0.42,0.59]
varCost    = [10,6,6,2]
emission        = [0,.399,.337,.201]

#As no CO2 price was given, we used the price from eex dated 26.05.20
CO2Cost    = 21.71

#Filling the DataFrame
data = DataFrame(technology=technology,FuelCost=FuelCost,efficiency=𝜂,OM=varCost,emissionFactor=λ)


#Not necessary steps, just to get used to csv :)
CSV.write("data.csv",data)
powerplants= CSV.read("data.csv")

#Calculating marginal cost by applying the formula on the columns elementwise
newcol = (FuelCost./efficiency)+CO2Cost*emission.+varCost

#attaching the results as column to already existing DataFrame
insertcols!(powerplants,:MC =>newcol)

powerplants

CSV.write("results.csv",powerplants)
