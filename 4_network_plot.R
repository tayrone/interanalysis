# This script plots a regulatory network using a rtni object.
# Obtaining the final relaxed network may be extremelly laborious. On the 
# other hand, the set of relax() instructions on this document makes it much
# easier, since its parameter values and iteraing structure were set 
# after several upsetting but illuminating experiences. Finally, this should
# be run at desktop, and not on a remote machine. You'll need your graphical 
# interface to manipulate the network.

library(RTN)
library(RedeR)
library(data.table)
library(RColorBrewer)
library(classInt)
library(gdata)

#load("g34_rtn.RData")

tfs <- rtni@regulatoryElements
tree <- tni.graph(rtni, tnet = "dpi", gtype = "amapDend", tfs = tfs)
tna <- tna.get(rtna, what = "mra")

#---- Defining network MRs ----

index <- which(tfs %in% tna$Regulon)
temp <- data.frame(gene = tfs[-index])
temp$is_mr <- "no"

temp <- rbind(temp, data.frame(gene = tna$Regulon, is_mr = "yes"))

#---- Initial plot ----

tal <- att.mapv(tree$g, dat = temp)
tal <- att.setv(g = tal, from = "name", to = "nodeAlias")

rdp <- RedPort()

calld(rdp)

addGraph(rdp, tal)


addLegend.color(rdp, colvec = c("#A5DD0E", "#BDBDBD"), 
                labvec = c("MR", "Not MR"), title = "node color")

#---- Relax it ----

#Let this run for a couple times. It works wonders.

for(x in 1:50){
  
  for(i in 1:1000){relax(rdp, p1 = 10, p2 = 100, p3 = 180, p4 = 180, 
                         p5 = 50, p6 = 10, p7 = 150, p8 = 200)}
  
  print("sleep!")
  
  if(x < 30)
    Sys.sleep(4) 
  Sys.sleep(8) # Initial iterations need more time in between
  
  print("wake up!")
}

for(x in 1:10){
  for(i in 1:1000){relax(rdp, p1 = 2, p2 = 100, p3 = 30, p4 = 30, p5 = 50, 
                         p6 = 10, p7 = 30, p8 = 50)}
  print("sleep!")
  Sys.sleep(1)
  print("wake up!")
}

#---- Selection of MR regulons, so it is possible to change its color ----

selectNodes(rdp, tna$Regulon)

#---- Size legend is necessary to read the amount of elements in regulon ----

scl <- round(as.numeric(tal$legNodeSize$scale))

leg <- as.character(round(as.numeric(tal$legNodeSize$legend)))

addLegend.size(rdp, sizevec = scl, labvec = leg, 
               title = "node size (Regulon size)", position = "bottomright", 
               intersp = 10, ftsize = 10, vertical = F)  


save.image("./interanalysis_files/4_net.RData")
