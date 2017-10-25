module App

include("lib.jl")
using .ClusterizationModule

# Reading data from .csv file
fileData = getDataSetFromCSV("creditcard.csv")

# Getting values from data as the array of coordinates 
xs = getDataFromColumnByName(fileData, "\"V1\"")
ys = getDataFromColumnByName(fileData, "\"V3\"")
# Uncomment this for only 'X' and 'Y' axis
zs = zeros(getDataFromColumnByName(fileData, "\"V6\""))
# zs = getDataFromColumnByName(fileData, "\"V6\"")

# Creating final dataset based on arrays of coordinates
dataset = linkData(xs, ys, zs)

# Creating points according to the following dataset
points = determinePoints(dataset)

# Defining count of the clusters
clustersCount = 6

# Creating inital clusters
clusters = defineInitialClusters(clustersCount, points)

counter = 1;
# Do clustering until clusters' centers will change
while(!kmeansIteration(clusters, points))  
    println("==============================")
    counter += 1
end
println(`Iterations: $counter`)

end