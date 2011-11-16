#!/usr/bin/python
print 'T\t\tT*'
for x in range(18,35):
	Ts = 5*x;
	print '%f\t%f' % (Ts, 8.617332e-5*(Ts+273.15))