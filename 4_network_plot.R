library(RTN)
library(networkD3)

load("./expression_analysis/rdata_files/network/g34_rtn.RData")

g <- rtn::tni.graph(rtni)

net <- igraph_to_networkD3()