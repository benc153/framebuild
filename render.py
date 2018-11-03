from numpy import *
from utils import *
import PIL.Image as Image
import PIL.ImageDraw as ImageDraw
from pdb import set_trace as brk

class Render:
	def __init__(self, inf, sup, px_per_mm=2):
		self.margin = array([100, 100])

		self.px_per_mm = px_per_mm

		self.inf = self._project(inf)
		self.sup = self._project(sup)

		self.dims = [int(x * px_per_mm) \
				for x in self.sup - self.inf + 2*self.margin]
		self.im = Image.new("RGB", self.dims, "white")
		self.draw = ImageDraw.Draw(self.im)

	def _project(self, p):
		return array([p[0], p[1]])

	def screen_point(self, p):
		p = self._project(p)
		ret = self.px_per_mm * (p - self.inf + self.margin)
		ret = [int(x) for x in ret]
		ret[1] = self.dims[1] - ret[1]
		return tuple(ret)

	def draw_point(self, p, fill="black", margin=4):
		"""Draw a box around a point, and return the screen point"""
		p = self.screen_point(p)
		self.draw.rectangle((p[0]-margin, p[1]-margin,
			p[0]+margin, p[1]+margin), fill)
		return p

	def draw_tube(self, tube):
		bottom = tube.tube_bottom
		top = tube.tube_top

		b = self.draw_point(bottom)
		t = self.draw_point(top)
		self.draw.line(b + t, fill="red")

		radii = [tube.get_radius(end) for end in (BOTTOM, TOP)]
		for r, colour in zip(radii, ("blue", "green")):
			v = turn_left(tube.vecn) * r
			for s in (-1, 1):
				start = bottom + s * v
				end = top + s * v
				self.draw.line(
						self.screen_point(start) + self.screen_point(end),
						fill=colour)
			if radii[1] == radii[0]: break

	def draw_wheel(self, origin, diameter, tyre_height):
		origin = self._project(origin)
		r = diameter / 2
		self.circle(origin, r, "black")
		self.circle(origin, r + tyre_height, "red")

	def draw_mitre(self, mitre):
		for d in mitre.reference_points():
			t = mitre.parent
			p = t.tube_top - d * t.vecn
			w = turn_left(t.vecn)
			r = t.get_radius(TOP)
			start = p + w * r
			end = p - w * r
			self.draw.line(self.screen_point(start) +
					self.screen_point(end), fill="magenta")
		self.draw_point(mitre.corners[0], "orange")
		self.draw_point(mitre.corners[1], "lightgreen")

	def polyline(self, points, colour):
		self.draw.line(tuple([self.screen_point(p) for p in points]), colour)

	def circle(self, origin, radius, colour):
		w = array([radius, radius])
		a, b = self.screen_point(origin - w), self.screen_point(origin + w)
		e_inf = (min(a[0], b[0]), min(a[1], b[1]))
		e_sup = (max(a[0], b[0]), max(a[1], b[1]))
		self.draw.ellipse(e_inf + e_sup, outline=colour)

	def save(self, outfile):
		self.im.save(outfile)
