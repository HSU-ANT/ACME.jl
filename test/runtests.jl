using ACME
using Base.Test

# Pathological cases for topomat:
# two nodes, one loop branch (short-circuited) -> voltage==0, current arbitrary
@test ACME.topomat(spzeros(Int, 2, 1)) == (sparse([1]), sparse([]))
# two nodes, one branch between them -> voltage arbitrary, current==0
@test ACME.topomat(sparse([1,2], [1,1], [1,-1])) == (sparse([]), sparse([1]))
