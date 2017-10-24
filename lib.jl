#=
    Module that provides clusterization algorithm
    for the custom datasets, defined in .csv format files
=#
module ClusterizationModule

#=
    Type of the point in the Cartesian coordinate system
=#
type Point{T<:Real}
    x::T
    y::T
    z::T
end

#=
    Type of the cluster with center and points,
    which belongs to its set
=#
struct Cluster
    points::Array{Point, 1}
    center::Point
end

#=
    Function that returns array of arrays
    of string values by reading of .csv file
=#
function getDataSetFromCSV(filename::String)::Array{Array{String, 1}, 1}
    arrayOfData::Array{Array{String, 1}, 1} = []
    open(filename) do file   
        for line::String in eachline(file)
          push!(arrayOfData, split(line, ','))  
        end    
    end   
    arrayOfData
end    

#=
    Function that converts array of arrays
    of string values to the 2-dimensional array
    of parsed values to float format 
=#
function convertToFloatData(array::Array{Array{String, 1}, 1},
     rowNumber::Int64, colNumber::Int64)::Array{Float64, 2}
    arrayOfFloats = Array{Float64, 2}(rowNumber, colNumber);   
    for i::Int64 = 1 : rowNumber
        for j::Int64 = 1 : colNumber
            arrayOfFloats[i, j] = parse(Float64, array[i][j]) 
        end
    end    
    arrayOfFloats
end

#=
    Function that gets index of the column from dataset
    by its name (if it contains it)    
=#
function getColumnIndexByName(array::Array{Array{String, 1}, 1},
     colName::String)::Int64 
    for i::Int64 = 1 : size(array[1], 1)
        if(array[1][i] == colName)           
            return i
        end
    end   
    return -1
end

#=
    Function that returns column's data by its index in dataset    
=#
function getDataFromColumnByIndex(array::Array{Array{String, 1}, 1},
     index::Int64)::Array{Float64, 1} 
    dataFromColumn = Array{Float64, 1}(0);
    for i::Int64 = 2 : size(array, 1)
        push!(dataFromColumn, parse(Float64, array[i][index])) 
    end   
    dataFromColumn
end

#=
    Function that returns column's data by its name in dataset    
=#
function getDataFromColumnByName(array::Array{Array{String, 1}, 1},
     colName::String)::Array{Float64, 1}       
    dataFromColumn = Array{Float64, 1}(0);
    colIndex::Int64 = getColumnIndexByName(array, colName) 
    for i::Int64 = 2 : size(array, 1)          
        push!(dataFromColumn, parse(Float64, array[i][colIndex])) 
    end   
    dataFromColumn
end

#=
    Function that links columns of the dataset into 2-d matrix
=#
function linkData(arrays...)::Matrix{Float64}
    result = Matrix{Float64}(size(arrays[1], 1), length(arrays))
    for i::Int64 = 1 : size(arrays[1], 1)
        for j::Int64 = 1 : length(arrays)
            result[i, j] = arrays[j][i]
        end
    end    
    result
end

#=
    Function that creates array of the points based on the chosen dataset,
    represented in 2-d matrix
=#
function determinePoints(dataset::Matrix{Float64})::Array{Point, 1}
    len::Int64 = size(dataset, 1)
    points = Array{Point, 1}(len)  
    for i::Int64 = 1 : len
        point = Point(dataset[i, 1], dataset[i, 2], dataset[i, 3])      
        points[i] = point      
    end  
    points
end

#=
    Function that creates array of points' projections
    based on the selected axis ("x", "y", "z")
=#
function getPointsProjections(points::Array{Point, 1},
     axis::String)::Array{Float64, 1}
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

export Point, Cluster, getDataSetFromCSV, convertToFloatData, getColumnIndexByName,
    getDataFromColumnByIndex, getDataFromColumnByName, linkData, determinePoints, getPointsProjections
end



