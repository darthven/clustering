module App

include("lib.jl")
using .ClusterizationModule

# Reading dataset from .csv file
fileData = getDataSetFromCSV("creditcard.csv")


xs = getDataFromColumnByName(fileData, "\"V1\"")
ys = getDataFromColumnByName(fileData, "\"V3\"")
zs = getDataFromColumnByName(fileData, "\"V6\"")

print(linkData(xs, ys, zs))

dataset = linkData(xs, ys, zs)
points = determinePoints(dataset)
println(getPointsProjections(points, "x"))


# Removing headers from dataset
# shift!(fileData)
# Converting string values to floats of the dataset
# dataset = convertToFloatData(fileData, 200, 10)
# println(dataset)

end