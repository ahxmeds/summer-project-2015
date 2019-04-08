#!/usr/bin/env python
import subprocess

gnuplot = subprocess.Popen(["gnuplot"],shell=False,stdin=subprocess.PIPE,)
for i in xrange(81,91):
	filename=str(i)+"/Hist_L=40_p=0."+str(i)+"0000"
	f=open(filename+".csv","r")
	x=[]
	for item in f:
		x.append((item.rstrip('\n')).split('\t'))
	f.close()
	y=[]	
	for item in x:
		z=[]
		if len(item)>0:
			for subitem in item:
				try:
					z.append(float(subitem))
				except:
					pass
			y.append(z)
	z=[]
	for item in y:
		if len(item)>0:
			z.append(item)
	filename=filename+"_c"
	f=open(filename+".csv","w")
	for item in z:
		f.write(str(item[0])+"\t"+str(item[1])+"\n")
	f.close()
	temp = 5
	for j in xrange(1,len(z)-15,16):
		print "folder=%d line=%d" %(i,j)
		gnuplot.stdin.write("set term jpeg\n")
		gnuplot.stdin.write("set output '"+filename+str(j)+".jpeg'\n")
		gnuplot.stdin.write("set yrange [0:160000]\n")
		gnuplot.stdin.write("set xrange [-0.5:15.5]\n")
		gnuplot.stdin.write("set xlabel 'Configuration No.'\n")
		gnuplot.stdin.write("set ylabel 'Counts'\n")
		gnuplot.stdin.write("unset key\n")
		gnuplot.stdin.write("set title 'T=%f'\n" % temp)
		temp-=0.01
		command='''"<(sed -n '%d,%dp' '''%(j,j+15)+filename+".csv"+''')"''' 
		gnuplot.stdin.write("plot "+command+" with boxes\n")			
			
