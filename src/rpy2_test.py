import rpy2.robjects as robjects
from rpy2.robjects.packages import importr







robjects.r("setwd('/home/vande060/Desktop/test')")
importr("edgeR")

datain = robjects.r("datain <- read.delim(file='test.csv', sep=',', row.names = 'Index')")
my_list = ['brett', 'is', 'cool']
group = robjects.r("group <- c({})".format(my_list))

dge = robjects.r("DGEList(datain, group = group, remove.zeros=TRUE)")

print(dge)












