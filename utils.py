from numpy import *

TOP, BOTTOM = range(2)
INSIDE, OUTSIDE = range(2)

def turn_right(vec):
	ret = vec.copy()
	ret[0] = vec[1]
	ret[1] = -vec[0]
	return ret

def turn_left(vec):
	ret = vec.copy()
	ret[0] = -vec[1]
	ret[1] = vec[0]
	return ret

def solve_quadratic(a, b, c):
	t = sqrt(b*b - 4*a*c)
	return (-b + t) / (2*a), (-b - t) / (2*a)

def rad2deg(theta):
	return (theta*360) / (2*pi)

def deg2rad(theta):
	return (theta*2*pi) / 360
