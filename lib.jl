#=
    Module that provides clusterization algorithm
    for the custom datasets, defined in .csv format files
=#
module ClusterizationModule

using DataFrames

#=
    Type of the point in the Cartesian coordinate system
=#
type Point{T <: Real}
    x::T
    y::T
    z::T
end

#=
    Type of the cluster with center and points,
    which belongs to its set
=#
type Cluster 
    center::Point
    points::Array{Point, 1}
end

#=
    Function that returns array of arrays
    of string values by reading of .csv file
=#
function datasetFromCSV(filename::String, customNastrings...)::DataFrame   
    readtable(filename, nastrings=["NA", "na", "n/a", "missing", customNastrings...])       
end    

#=
    Function that links columns of the dataset into 2-d matrix
=#
function datasetValues(dataset::DataFrame, columnNames...)::Matrix{Float64}    
    selectedData::DataFrame = dataset[:, map(Symbol, [columnNames...])] 
    values = Matrix{Float64}(nrow(selectedData), ncol(selectedData))
    # Excluding headers and rows' numbers
    for i::Int64 = 1 : nrow(selectedData)
        for j::Int64 = 1 : ncol(selectedData)
            values[i, j] = selectedData[i, j]
        end
    end
    values
end

#=
    Function that creates array of the points based on the chosen dataset,
    represented in 2-d matrix
=#
function determinePoints(datasetValues::Matrix{Float64})::Array{Point, 1}
    len::Int64 = size(datasetValues, 1)
    points = Array{Point, 1}(len)  
    for i::Int64 = 1 : len
        point = Point(datasetValues[i, 1], datasetValues[i, 2], datasetValues[i, 3])      
        points[i] = point      
    end  
    points
end

#=
    Function that creates array of points' projections
    based on the selected axis ("x", "y", "z")
=#
function pointsProjections(points::Array{Point, 1}, axis::String)::Array{Float64, 1}
    coords = Array{Float64, 1}(0)
    for p::Point in points
        if(lowercase(axis) == "x")
            push!(coords, p.x)
        elseif(lowercase(axis) == "y")
            push!(coords, p.y)
        elseif(lowercase(axis) == "z")
            push!(coords, p.z)
        else    
            push!(coords, 0.0)      
        end
    end
    coords
end

#=
    Function that returns point by random index
    from the array of points   
=#
function randomPoint(points::Array{Point, 1})::Point
    index::Int64 = rand(1: size(points, 1))
    points[index];
end

#=
    Function that creates initial clusters   
=#
function defineInitialClusters(clustersCount::Int64, points::Array{Point, 1})::Array{Cluster, 1}
    clusters = Array{Cluster, 1}(clustersCount)    
    for i::Int64 = 1 : clustersCount
        center::Point = randomPoint(points)
        cluster = Cluster(center, [center])
        clusters[i] = cluster
    end
    clusters
end

#=
    Function that returns distance between two points 
=#
function distance(point1::Point, point2::Point)::Float64   
    norm([point1.x, point1.y, point1.z] - [point2.x, point2.y, point2.z])    
end

#=
    Function that returns distances' matrix based on distances between
    clusters' centers and all the points
=#
function distanceMatrix(points::Array{Point, 1}, clusters::Array{Cluster, 1})::Matrix{Float64}
    distanceMatrix = Matrix{Float64}(length(points), length(clusters))
    for i::Int64 = 1 : length(points)
        for j::Int64 = 1 : length(clusters)
            distanceMatrix[i, j]::Float64 = distance(points[i], clusters[j].center)
        end
    end    
    distanceMatrix
end

#=
    Function that compares two points by their coordinates  
=#
function comparePointsCoordinates(point1::Point, point2::Point)::Bool
   point1.x == point2.x && point1.y == point2.y && point1.z == point2.z
end

#=
    Function that returns new center of the cluster as the center of mass   
=#
function centerMassOfCluster(cluster::Cluster)::Point   
    x::Float64 = reduce(+, pointsProjections(cluster.points, "x")) / length(cluster.points)
    y::Float64 = reduce(+, pointsProjections(cluster.points, "y")) / length(cluster.points)
    z::Float64 = reduce(+, pointsProjections(cluster.points, "z")) / length(cluster.points)    
    Point(x, y, z)    
end

#=
    Function for displaying clusters
=#
function displayCluster(clusters::Array{Cluster, 1})::Void
    # Display clustering iteration
    for i::Int64 = 1 : length(clusters)
        println("CLUSTER $i:")
        println(clusters[i])
        println("SIZE: $(length(clusters[i].points))")
    end
end

#=
    Function that makes one iteration of the clustering 
=#
function kmeansIteration(clusters::Array{Cluster, 1}, points::Array{Point, 1})::Bool 
    status::Bool  = false;     
    dtMatrix::Matrix{Float64} = distanceMatrix(points, clusters)      
    for i::Int64 = 1 : size(dtMatrix, 1)
        minIndex::Int64 = indmin(dtMatrix[i, 1 : size(dtMatrix, 2)])  
        if(!(points[i] in clusters[minIndex].points))            
            push!(clusters[minIndex].points, points[i])    
        end
        for j::Int64 = 1 : size(dtMatrix, 2)
            if(j != minIndex && (points[i] in clusters[minIndex].points))                                 
                clusters[j].points = filter!(pnt -> pnt != points[i], clusters[j].points)                     
            end
        end      
    end
    for cl::Cluster in clusters
        newCenter::Point = centerMassOfCluster(cl)
        if(comparePointsCoordinates(cl.center, newCenter))
           status = true; 
        else
           cl.center = newCenter
        end
    end
    # displayClustersInfo(clusters)
    status   
end

#=
    Function that makes the whole clustering by the k-means algorithm
=#
function kmeans(clustersCount::Int64, datasetValues::Matrix{Float64})::Array{Cluster, 1}
    iteration::Int64 = 0;
    # Creating points based on dataset
    points::Array{Point, 1} = determinePoints(datasetValues)
    # Creating inital clusters
    clusters::Array{Cluster, 1} = defineInitialClusters(clustersCount, points)
    # Do clustering until clusters' centers will change
    while(!kmeansIteration(clusters, points))  
        iteration += 1
    end
    println("Iterations: $iteration") 
    return clusters
end

export Point, Cluster, datasetFromCSV, datasetValues, kmeans

end



