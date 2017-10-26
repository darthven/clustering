module App

include("lib.jl")

using .ClusterizationModule

# Reading data from .csv file
dataset = datasetFromCSV("creditcard.csv")

# Getting values for dataset by column names
values = datasetValues(dataset, "V1", "V3", "V6")

# Defining count of the clusters
clustersCount = 6

# Running k-means clustering algorithm
kmeans(clustersCount, values)

end