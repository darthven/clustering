module App

include("lib.jl")

using .ClusterizationModule

# Reading data from .csv file
dataset = datasetFromCSV("creditcard.csv")

# Getting values for dataset by column names
println("Enter 1st column name, please")
col1 = readline()
println("Enter 2nd column name, please")
col2 = readline()
println("Enter 3rd column name, please")
col3 = readline()
values = datasetValues(dataset, col1, col2, col3)

# Uncomment if need hardcoded column names
# values = datasetValues(dataset, "V1", "V3", "V6")

# Defining count of the clusters
println("Enter clusters count, please")
clustersCount = parse(Int64, readline())

# Uncomment if need hardcoded clusters' count
# clustersCount = 6

# Running k-means clustering algorithm
kmeans(clustersCount, values)

end
