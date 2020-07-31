#!/usr/bin/env python3
from argparse import *
from math import *
import sys
import PIL.Image as Image
import PIL.ImageDraw as ImageDraw
from utils import rad2deg, deg2rad
from pdb import set_trace as brk

def distance(a, b):
	"""The distance between two points a and b"""
	ret = 0
	for i in range(2): ret += (b[i] - a[i])**2
	return sqrt(ret)

def angles(n=360*4):
	"""Split the circle into n pieces and yield the angle at each one"""
	theta = 0
	step = (2*pi) / n

	for _ in range(n):
		yield theta
		theta += step

class Ellipse:
	def __init__(self, diameters):
		self.diameters = diameters
		self.radii = [d/2 for d in diameters]
		self.perim = self._walkaround()

	def _walkaround(self):
		"""Walk around the ellipse and return an array of the perimeter at each
		angular increment."""
		ret = []

		total_dist, prev = 0, None

		for theta in angles(100000):
			point = (self.radii[0] * cos(theta), self.radii[1] * sin(theta))
			if prev:
				total_dist += distance(prev, point)
			ret.append((theta, total_dist))
			prev = point

		return ret

	def _search(self, angle, perim):
		n = len(perim)

		if n == 1:
			return perim[0][1]

		m = n // 2
		if angle < perim[m][0]:
			return self._search(angle, perim[:m])
		else:
			return self._search(angle, perim[m:])

	def circum_distance(self, theta):
		"""Given an angle theta around the ellipse, return the distance around
		the circumference at that point. We do this with a precomputed map and
		a binary search. There is better math to do this but it's harder than
		you'd think"""
		theta = fmod(theta, 2*pi)
		return self._search(theta, self.perim)

	def text(self):
		return u"{}mm, {}mm".format(*self.diameters)

class Circle(Ellipse):
	def __init__(self, diameter):
		self.diameter = diameter
		self.radius = diameter / 2

	def circum_distance(self, theta):
		return theta * self.radius

	def text(self):
		return "{}mm".format(self.diameter)

class Curve:
	def __init__(self, cut, par, wall,
			angle, rotate=False, offset=0, taper=0.0, twist=0.0):
		"""cut is diameters of cut tube, par of parent tube (or just one
		diameter if they're the same), wall is wall thickness. All units mm.
		angle is the angle between them in degrees. If rotate, display the plot
		shifted by 90 degrees around the tube. taper is the taper of the parent
		tube, expressed as a ratio of diameter change per unit length"""
		if not hasattr(cut, "__iter__"):
			cut = [cut, cut]
		else:
			self.cut = list(cut)

		self.par = par
		self.wall = wall
		self.angle = angle
		self.rotate = rotate
		self.offset = offset
		self.taper = taper
		self.twist = twist

		if cut[0] == cut[1]:
			self.ellipse = Circle(cut[0])
		else:
			self.ellipse = Ellipse(cut)

		i = float('inf')
		self.inf, self.sup = [i, i], [-i, -i]
		self.make_points()

	def make_points(self):
		cut, par, wall, angle = self.cut[:], self.par, self.wall, self.angle
		points = []

		# We use a coordinate system in which the diameter of the parent tube
		# is always 1.0 (or starts at 1.0 in the case of a taper) R is then the
		# radii of the inside of the cut tube.
		R = [0, 0]
		for i in range(2):
			cut[i] -= wall * 2
			R[i] = cut[i] / par

		offset = self.offset / par

		# The angle is the angle of the tube axes, but we want the orientation
		# of the parent tube relative to the cut tube's cross-section, so just
		# add 90 degrees.
		angle = fmod(angle + 90, 360)
		slope = tan(deg2rad(angle))
		twist = deg2rad(self.twist)

		# Radius as opposed to diameter
		par_r = par / 2

		# How much the radius changes by
		taper_r = self.taper / 2

		min_y = 0
		for theta in angles():

			# (a, b) is the point on the inside circumference of the cut tube
			# at theta
			a = R[0] * cos(theta + twist)
			b = R[1] * sin(theta + twist)

			a += offset

			# Now we're projecting the point onto the parent tube
			if a > 1 or a < -1:
				# In this case we miss the parent tube altogether, so no
				# cutting here.
				y = 1
			else:
				# phi is the angle of this point around the parent tube.
				phi = asin(a)
				r = 1.0 + b * taper_r
				y = 1 - r * cos(phi) - b * slope

			# We map the angle back onto the outside diameter, although the
			# actual curve is on the inside diameter. This implies your cut is
			# perpendicular to the tube wall.
			point = [self.ellipse.circum_distance(theta), par_r * y]
			min_y = min(min_y, point[1])
			points.append(point)

		# Now push everthing up to a couple of inches above 0 just to have a
		# bigger piece of paper to wrap
		headroom = 25.4 * 2
		if min_y < 0:
			headroom -= min_y

		for p in points:
			p[1] += headroom

			# Track the infimum and supremum of the full set of points
			for i in range(2):
				self.inf[i] = min(self.inf[i], p[i])
				self.sup[i] = max(self.sup[i], p[i])

		self.points = points

	def render_png(self, outfile, caption=None, pxi=100):
		margin = 25.4
		scale = (pxi / 25.4)
		dims = [self.sup[i] for i in range(2)]
		shape = [int(scale * (dims[i] + margin * 2)) for i in range(2)]
		im = Image.new("RGB", shape, "white")
		draw = ImageDraw.Draw(im)

		def screen_point(p):
			return [scale * (p[i] + margin) for i in range(2)]

		# Offset everything by a 3/4 turn for aesthetic reasons and so
		# it looks the same as the one on metalgeek.com
		if not self.rotate:
			offset = (3*dims[0]) / 4
		else:
			offset = 0

		prev = None
		for p in self.points:
			p = p[:]
			p[0] = fmod(p[0] + offset, dims[0])
			sp = screen_point(p)

			if prev and prev[0] < sp[0]:
				draw.line(prev + sp, "black")
			else:
				draw.point(sp, "black")

			prev = sp

		# Draw some useful guidelines
		pos, step = 0, dims[0] / 4
		for i in range(5):
			a = [pos, 0]
			b = [pos, self.sup[1]]
			coords = screen_point(a) + screen_point(b)
			draw.line(coords, "red")
			pos += step

		legend = "Cut: {}\nParent: {}mm\n" \
				"Wall: {}mm\nAngle: {:.2f}°".format(self.ellipse.text(),
						self.par, self.wall, self.angle)
		if self.offset:
			legend = "{}\nOffset: {}mm".format(legend, self.offset)
		if self.taper:
			legend = "{}\nTaper: {}mm".format(legend, self.taper)
		if caption:
			legend = "{}\n{}".format(caption, legend)
		if self.twist:
			legend = "{}\nTwist: {:.2f}°".format(legend, self.twist)

		draw.text(((margin + 3) * scale, margin * scale),
				legend, fill="black")

		im.save(outfile)
		print("Rendered to {}".format(outfile))

def parse_cut(s):
	r = [float(x.strip()) for x in s.split(',')]

	if len(r) == 1:
		return [r[0], r[0]]
	elif len(r) == 2:
		return r
	else:
		raise ArgumentTypeError("I don't understand your cut diameter(s)")

def main():
	ap = ArgumentParser(description="""
Produces templates for coping tubes, at a default resolution of 100px/in.

You create your template, print it out, then cut it out, and wrap it around the
tube. Trace the curvy edge with a pen or a pencil and cut there.

Can handle elliptical cut tubes with either ellipse axis aligned with the
parent tube.
""")

	ap.add_argument("-c", "--cut", type=str,
			default="25.4",
			help="""Cut tube diameter or two ','-separated diameters in mm if
			cut tube is elliptical. Put the bigger diameter first if it's
			oriented like a chainstay, the smaller one first if it's like a
			farm gate. Default 25.4.""")
	ap.add_argument("-p", "--parent", type=float,
			default="25.4",
			help="Parent tube diameter in mm, default 25.4.")
	ap.add_argument("-w", "--wall",
			type=float, default=0.9,
			help="Wall thickness in mm, default 0.9")
	ap.add_argument("-a", "--angle", type=float, default=90,
			help="Angle between the tubes in degrees, default 90")
	ap.add_argument("-b", "--twist", type=float, default=0,
			help="Parent tube twist in degrees")
	ap.add_argument("-r", "--resolution", type=float, default=100,
			help="Pixels per inch, default 100")
	ap.add_argument("-e", "--taper", type=float, default=0.0,
			help="Taper of parent tube as ratio of diameter change to "
			"length change")
	ap.add_argument("-o", "--outfile", default="template.png")
	ap.add_argument("-f", "--rotate", action="store_true", default=False)
	ap.add_argument("-s", "--offset", type=float, default=0.0)
	ap.add_argument("-t", "--caption")

	try:
		args = ap.parse_args()
		cut = parse_cut(args.cut)
		if args.angle < 1 or args.angle > 129:
			raise ArgumentTypeError("Angle too big or small")
	except ArgumentTypeError as error:
		print(error, file=sys.stderr)
		sys.exit(1)

	curve = Curve(cut, args.parent, args.wall,
			args.angle, args.rotate, args.offset, args.taper, args.twist)
	curve.render_png(args.outfile, args.caption, args.resolution)

if __name__ == "__main__":
	main()
