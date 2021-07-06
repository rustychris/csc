# HEC-DSSVUE script to dump all of the DCD data to xlsx files.

import os
import java
from hec.heclib.dss import HecDss
from hec.dataTable import HecDataTableToExcel

print "Hi"

hecfile=HecDss.open("dcd_sep2016.dss",1) # must exist

#dest="excel"
dest="csv" # way faster on write and read, and slightly smaller.

if not os.path.exists(dest):
    os.makedirs(dest)

ABCs={}
count=0
for path in hecfile.getPathnameList():
    parts=path.split('/')
    parts[4]=''
    ABC="/".join( parts ) # really ABCEF
    if ABC in ABCs: continue
    ABCs[ABC]=1
    print "Path: %s"%ABC

    data=hecfile.get(path,True) # True to read full time series
    
    if dest=="excel":
        # Add Data
        datasets = java.util.Vector()
        datasets.add(data)
        # For this code, jython sees a List before a Vector
        list = []
        list.append(datasets)
        table = HecDataTableToExcel.newTable()
        target= os.path.join(dest,"export-path%07d.xlsx"%count)
        print "Writing %d values to %s"%(data.numberValues,target)
        table.createExcelFile(list,target)
    elif dest=="csv":
        target= os.path.join(dest,"export-path%07d.csv"%count)

        with open(target,'w') as fp:
            fp.write("# path=%s units=%s type=%s\n"%(ABC,data.units,data.type))
            fp.write("time,value\n")
            times = data.getTimes()
            for i,Q in enumerate(data.values):
                t = times.element(i)
                #file.write(str(t.year())+"-"+str(t.month())+"-"+str(t.day())+ " "+str(t.hour())+":"+str(t.minute())+"\n")
                fp.write(t.date(-13)+",%.2f\n"%Q)
        
    count+=1
    
hecfile.done()




