from argparse import *
from numpy import *
from numpy.linalg import norm
import configparser
from copecalc import Curve
from utils import *
from copy import *
from render import Render
from pdb import set_trace as brk

EPSILON = 1e-5
seterr(all='raise')

class Tube:
	def __init__(self, diameter, wall):
		self.diameter = diameter
		self.radius = diameter / 2
		self.wall = wall

	def place(self, top, bottom, tube_top=None, tube_bottom=None):
		"""top and bottom are the positions of the intersections, tube_top and
		tube_bottom are the positions of the actual ends of the tube (on the
		centreline)."""
		self.top = top
		self.bottom = bottom

		self.tube_top = tube_top if tube_top is not None else top
		self.tube_bottom = tube_bottom if tube_bottom else bottom

		self.vec = self.top - self.bottom
		self.vecn = self.vec / norm(self.vec)

	def get_end(self, end):
		if end == TOP:
			return self.top
		elif end == BOTTOM:
			return self.bottom
		assert False

	def get_diameter(self, end):
		return self.diameter

	def get_diameters(self, end):
		ret = self.get_diameter(end)
		return (ret, ret)

	def get_radius(self, end):
		return self.radius

	def get_wall(self, end):
		return self.wall

	def length(self):
		return abs(dot(self.tube_top - self.tube_bottom, self.vecn))

	def transform(self, mat):
		"""Return a transformed copy"""
		ret = deepcopy(self)

		for attr in ("top", "bottom", "tube_top", "tube_bottom", "vec", "vecn"):
			v = getattr(self, attr)
			setattr(ret, attr, dot(mat, v))

		ret.vecn /= norm(ret.vecn)
		return ret

class ExternallyButtedTube(Tube):
	def __init__(self, diameters, walls):
		self.diameters = diameters
		self.radii = [d/2 for d in diameters]
		self.walls = walls

	def get_diameter(self, end):
		return self.diameters[end]

	def get_radius(self, end):
		return self.radii[end]

	def get_wall(self, end):
		return self.walls[end]

class EllipticalTube(Tube):
	def __init__(self, diameters, wall):
		self.minor_diameter = diameters[1]
		super(EllipticalTube, self).__init__(diameters[0], wall)

	def get_diameters(self, end):
		return (self.diameter, self.minor_diameter)

class Dropout:
	def __init__(self, thickness, cs_length, ss_length):
		self.thickness = thickness
		self.cs_length = cs_length
		self.ss_length = ss_length

	def place(self, cs):
		"""cs is the chainstay"""
		self.start = cs.top
		self.end = cs.top + self.cs_length * cs.vecn
		self.vec = self.end - self.start
		self.vecn = self.vec / norm(self.vec)

class Mitre:
	def __init__(self, cut, cut_end, parent, parent_end):
		"""cut and parent are tubes, a_end and b_end tell you which ends you're
		putting the mitre between. The two tubes must be rotated to be in the
		z=0 plane."""
		self.cut = cut
		self.cut_end = cut_end

		self.parent = parent
		self.parent_end = parent_end
		self.angle = arccos(dot(self.cut.vecn, self.parent.vecn))

		# The "radius" of the mitre is the length from the centre line to the
		# top or the bottom of the cut.
		r = self.cut.get_radius(cut_end)
		self.radius = r / dot(turn_left(self.cut.vecn), self.parent.vecn)
		self._find_corners()

	def make_template(self, name, caption, rotate):
		cut = self.cut.get_diameters(self.cut_end)
		cut_wall = self.cut.get_wall(self.cut_end)

		par = self.parent.get_diameter(self.parent_end)
		curve = Curve(cut, par, cut_wall, rad2deg(self.angle), rotate)
		curve.render_png(name + ".png", caption)

	@staticmethod
	def _away_vec(tube, end):
		"""Get the normalized vector from this end of the tube to the other
		end"""
		this_end = tube.get_end(end)
		other_end = tube.get_end((1 + end) % 2)
		ret = other_end - this_end
		return ret / norm(ret)

	def _find_corners(self):
		"""The points on the cut tube that touch the parent tube"""
		par_r = self.parent.get_radius(self.parent_end)
		cut_r = self.cut.get_radius(self.cut_end)

		par_v = self._away_vec(self.parent, self.parent_end)
		cut_v = self._away_vec(self.cut, self.cut_end)

		# How far to travel down the cut tube to leave the parent
		w = turn_left(par_v)
		if dot(w, cut_v) < 0: w *= -1
		a = par_r / dot(w, cut_v)

		# How far to travel down the parent tube to leave the cut tube
		w = turn_left(cut_v)
		if dot(w, par_v) < 0: w *= -1
		b = cut_r / dot(w, par_v)

		centre = self.parent.get_end(self.parent_end)

		# The first corner is the inside one, the second one the outside one.
		# More by luck than judgment
		self.corners = (centre + a*cut_v + b*par_v,
				centre + a*cut_v - b*par_v)

	def reference_points(self):
		"""Three distances measured from the top of the tube, to the top of the
		mitre, the centre, and the bottom"""
		ret = []
		middle = (self.corners[0] + self.corners[1]) / 2
		for point in (middle,) + self.corners:
			ret.append(dot(point - self.parent.tube_top, -self.parent.vecn))
		ret.sort()
		return ret

class RearTriangle:
	def __init__(self, frame, config):
		self.frame = frame

		self._place_tubes()

		# Draw things from a top view, and portray z on the x axis
		self.projection = array([[0, 0, -1], [1, 0, 0], [0, 0, 0]])

		th = arctan(frame.left_cs.vecn[1] / frame.left_cs.vecn[0])
		self.rotation = array([
			[cos(th), -sin(th), 0],
			[sin(th), cos(th), 0],
			[0, 0, 1]]).transpose()
		self.projection = dot(self.projection, self.rotation)

		self._calc_mitres()

	def _ss_transform(self):
		"""Transform to put the left SS and ST into the z=0 plane"""
		f = self.frame
		v = f.left_ss.vecn
		w = cross(v, f.seat_tube.vecn); w /= norm(w)
		u = cross(v, w)

		return array([u, v, w])

	def _place_css(self):
		f = self.frame

		bb_radius = f.bb_tube.get_radius(BOTTOM)
		rear_offset = (f.rear_axle_spacing + f.dropout_thickness) / 2

		# Chainstay length is defined as the length on the centre-line as if
		# the CS ran from the centre of the BB shell to the centre of the rear
		# axle. So use that to solve for the end position.
		x = sqrt(f.cs_length**2 - f.bb_drop**2)
		end = array([-x, f.bb_drop, rear_offset])

		# Now find the start to get bb_ext right
		pos = 0
		x = array([-1, 0, 0])
		c = f.bb_shell/2 - f.bb_ext

		while True:
			start = array([0, 0, pos])
			v = end - start
			v /= norm(v)

			# Angle between CS and BB shell
			angle = arccos(dot(v, array([0, 0, 1])))

			# "Radius" of the CS at the angle it makes with the BB shell
			rad = f.left_cs.get_radius(BOTTOM) / sin(angle)

			error = c - (pos + bb_radius * cos(angle) + rad)
			if abs(error) < EPSILON:
				break
			pos += error

		# Now adjust the end position for the dropout.

		# The DO runs in the same direction as the CS in the xy plane but is
		# flat in the z.
		w = v * -1
		w[2] = 0
		w *= f.left_drop.cs_length / norm(w)
		end += w

		f.left_cs.place(end, start)

		# The right one is the same only reflected in the xy plane
		start, end = start.copy(), end.copy()
		start[2] *= -1
		end[2] *= -1

		f.right_cs.place(end, start)

	def _update_bb(self):
		"""Now that we know where the CS intersects the BB we can work out its
		Tube Top and Tube Bottom properly"""
		f = self.frame
		top, bottom = f.bb_tube.top, f.bb_tube.bottom

		f.bb_tube.top = f.right_cs.bottom
		f.bb_tube.bottom = f.left_cs.bottom

		f.bb_tube.tube_top = top
		f.bb_tube.tube_bottom = bottom

	def _place_dropouts(self):
		f = self.frame
		f.left_drop.place(f.left_cs)
		f.right_drop.place(f.right_cs)

	def _place_sss(self):
		f = self.frame
		t = f.top_tube.bottom
		top = array([t[0], t[1], 0])

		bottom = f.left_drop.end.copy()
		v = bottom - top
		bottom -= (v * f.left_drop.ss_length) / norm(v)
		f.left_ss.place(top, bottom)

		bottom = bottom.copy()
		bottom[2] *= -1
		f.right_ss.place(top, bottom)

	def _calc_mitres(self):
		f = self.frame
		cs = f.left_cs.transform(self.projection)
		bb = f.bb_tube.transform(self.projection)
		self.cs_bb = Mitre(cs, BOTTOM, bb, BOTTOM)
		self.cs_bb.make_template("cs_bb", "Chain Stay to BB Shell", False)

		self.cs_in_mitre_length = abs(dot(cs.top - self.cs_bb.corners[INSIDE],
			cs.vecn))
		self.cs_out_mitre_length = abs(dot(cs.top - self.cs_bb.corners[OUTSIDE],
			cs.vecn))

		t = self._ss_transform()
		ss = f.left_ss.transform(t)
		st = f.seat_tube.transform(t)
		self.ss_st = Mitre(ss, TOP, st, TOP)
		self.ss_st.make_template("ss_st", "Seat Stay to Seat Tube", False)

		self.ss_in_mitre_length = abs(dot(ss.bottom - self.ss_st.corners[INSIDE],
			ss.vecn))
		self.ss_out_mitre_length = abs(dot(ss.bottom - self.ss_st.corners[OUTSIDE],
			ss.vecn))

	def _place_tubes(self):
		self._place_css()
		self._update_bb()
		self._place_dropouts()
		self._place_sss()

	def _draw_tyre(self, r):
		f = self.frame
		centre = (f.left_drop.end + f.right_drop.end) / 2
		v = f.left_cs.top - f.right_cs.top
		v /= norm(v)

		left = centre - v * (f.tyre_width / 2)
		right = left + v * f.tyre_width
		rad = f.wheel_diameter / 2 + f.tyre_height

		ends = []
		w = array([1, 0, 0])
		for p in (left, right):
			start = dot(self.projection, p)
			end = dot(self.projection, p + w * rad)
			ends.append(end)
			r.polyline((start, end), "black")
		r.polyline(ends, "black")

	def _draw_dropouts(self, r):
		f = self.frame
		for drop in (f.left_drop, f.right_drop):
			p = drop.start
			p[2] += drop.thickness / 2

			z = array([0, 0, -drop.thickness])
			x = array([-drop.cs_length, 0, 0])

			points = [p, p + z, p + z + x, p + x, p]
			points = [dot(self.projection, p) for p in points]

			r.polyline(points, "black")

		# Draw the axle as well
		start = f.left_drop.end
		end = f.right_drop.end
		points = [dot(self.projection, p) for p in (start, end)]
		r.polyline(points, "black")

	def _draw_chainrings(self, r):
		f = self.frame
		n = float(len(f.chainring_radii))

		# The distance from the inside to the middle of the chainset
		if n > 1:
			m = f.chainring_spacing * (n / 2)
		else:
			# This should be half the thickness of the chainring. We'll
			# over-estimate a bit to be safe.
			m = 2.5

		offset = f.chainline - m
		for rad in f.chainring_radii:
			line_start = dot(self.projection, array([-rad, 0, offset]))
			line_end = dot(self.projection, array([rad, 0, offset]))
			r.polyline((line_start, line_end), "green")
			offset += f.chainring_spacing

	def display(self):
		f = self.frame

		print("""
Rear Triangle
=============

Tube Cuts
---------""")

		print("Chain Stay length from inside mitre to dropout: {:.2f}".format(
			self.cs_in_mitre_length))

		print("Chain Stay length from outside mitre to dropout: {:.2f}".format(
			self.cs_out_mitre_length))

		print("Seat Stay length centre-centre: {:.2f}".format(
			f.left_ss.length()))

		print("Seat Stay length from inside mitre to dropout: {:.2f}".format(
			self.ss_in_mitre_length))

		print("Seat Stay length from outside mitre to dropout: {:.2f}".format(
			self.ss_out_mitre_length))

		print("Angle between CS and BB: {:.2f}deg".format(
			rad2deg(self.cs_bb.angle)))

		print("Angle between CS and dropout: {:.2f}deg".format(
			rad2deg(pi - self.cs_bb.angle + pi/2)))

		print("Angle between SS and ST: {:.2f}deg".format(
			rad2deg(arccos(dot(f.left_ss.vecn, f.seat_tube.vecn)))))

	def render_scale_diagram(self, outfile, px_per_mm=2):
		f = self.frame

		left_cs = f.left_cs.transform(self.projection)
		right_cs = f.right_cs.transform(self.projection)
		bb = f.bb_tube.transform(self.projection)

		inf = left_cs.top
		sup = bb.top

		r = Render(inf, sup, px_per_mm)

		r.draw_tube(bb)
		r.draw_tube(left_cs)
		r.draw_tube(right_cs)

		self._draw_tyre(r)
		self._draw_dropouts(r)
		self._draw_chainrings(r)
		r.draw_mitre(self.cs_bb)

		r.save(outfile)

class Frame:
	def __init__(self, config):
		self._parse_section(config["Basic Dimensions"], (
			("st", "seat tube length"),
			("tt", "top tube length"),
			("bb_drop", "bb drop"),
			("bb_shell", "bb shell"),
			("cs_length", "chain stay length"),
			))

		self._parse_section(config["Angles"], (
			("seat_angle", "seat"),
			("head_angle", "head"),
			("top_angle", "top"),
			), convert=lambda s: deg2rad(float(s)))

		self._parse_section(config["Fork"], (
			("fork_len", "fork length"),
			("fork_offset", "fork offset"),
			("lower_stack", "lower stack"),
			))

		self._parse_section(config["Extensions"], (
			# By "edge" we mean these are the distances from the top of the TT
			# where it touches the HT and from the bottom of the DT where it
			# touches the HT.
			("upper_ht_ext_edge", "head tube upper"),
			("lower_ht_ext_edge", "head tube lower"),
			("st_stickout", "seat tube"),
			("bb_ext", "bb shell"),
			))

		diams = {k: float(v) for k, v in dict(
			config["Tube Diameters"]).items()}
		walls = {k: float(v) for k, v in dict(
			config["Tube Walls"]).items()}

		self.top_tube = Tube(diams["top tube"], walls["top tube"])
		self.head_tube = Tube(diams["head tube"], walls["head tube"])
		self.down_tube = Tube(diams["down tube"], walls["down tube"])
		self.bb_tube = Tube(diams["bb shell"], walls["bb shell"])

		self.left_cs = EllipticalTube((diams["chain stay major"],
			diams["chain stay minor"]), walls["chain stay"])
		self.right_cs = deepcopy(self.left_cs)

		self.left_ss = Tube(diams["seat stay"], walls["seat stay"])
		self.right_ss = deepcopy(self.left_cs)

		self.seat_tube = ExternallyButtedTube(
			(diams["seat tube top"], diams["seat tube bottom"]),
			(walls["seat tube top"], walls["seat tube bottom"]))

		self._parse_section(config["Wheels"], (
			("rear_axle_spacing", "rear axle spacing"),
			("wheel_diameter", "diameter"),
			("tyre_height", "tyre height"),
			("tyre_width", "tyre width"),
			))

		self._parse_section(config["Dropouts"], (
			("dropout_thickness", "thickness"),
			("dropout_cs_length", "cs_length"),
			("dropout_ss_length", "ss_length"),
			))

		self._parse_section(config["Chainrings"], (
			("chainline", "chainline"),
			("chainring_spacing", "chainring spacing"),
			))

		chainrings = [int(x) for x in config["Chainrings"]["teeth"].split(',')]
		self.chainring_radii = [(x * 25.4) / (4*pi) for x in chainrings]

		self.left_drop = Dropout(self.dropout_thickness,
				self.dropout_cs_length, self.dropout_ss_length)
		self.right_drop = Dropout(self.dropout_thickness,
				self.dropout_cs_length, self.dropout_ss_length)

		self._place_tubes()
		self._calc_mitres()
		self.rear_triangle = RearTriangle(self, config)

	def _parse_section(self, config_section, mapping, convert=float):
		for our_name, their_name in mapping:
			setattr(self, our_name, convert(config_section[their_name]))

	def _place_st(self):
		bottom = array([0, 0, 0])
		angle = self.seat_angle
		vec = array([-cos(angle), sin(angle), 0])
		top = bottom + vec * self.st
		tube_top = bottom + vec * (self.st + self.st_stickout)
		self.seat_tube.place(top, bottom, tube_top)

	def _place_tt(self):
		"""The "bottom" is the seat tube end, although we assume the TT is
		actually horizontal"""
		bottom = self.seat_tube.top
		top = bottom + self.tt * \
				array([cos(self.top_angle), sin(self.top_angle), 0])
		self.top_tube.place(top, bottom)

	def _place_ht_dt(self):
		"""HT and DT are placed together, as we need to iterate to get the
		desired stickout below the HT, as it depends on the angle with the
		DT"""
		top = self.top_tube.top
		angle = self.head_angle
		vec = array([cos(angle), -sin(angle), 0])

		# Fork rake seems to be measured perpendicular to the head-tube line.
		rake_vec = turn_left(vec)
		fork_vec = vec * self.fork_len + rake_vec * self.fork_offset

		# x is the distance from crown_ix to the bottom of the head-tube. Our
		# input value instead is the distance from the bottom of the DT where
		# it touches the HT. We need to iterate a bit to find x as it depends
		# on the angle of the DT which also depends on the position of
		# crown_ix. You could solve this without iterating I'm sure but the
		# math is a bit harder.
		x = 0
		while True:
			# top + (la + lower_stack + x) * vec + fork_vec takes us to the axle
			# line. The crown intersection is at top + la * vec.
			la = (self.bb_drop - top[1] - fork_vec[1]) / vec[1] \
					- self.lower_stack - x

			# Provisional position of intersection of DT and HT
			crown_ix = top + la * vec

			d = crown_ix - self.seat_tube.bottom
			d /= norm(d)

			# How far do we need to travel down DT from crown_ix to leave the HT?
			a = self.head_tube.radius / dot(turn_right(vec), -d)

			# And how far down HT to leave DT
			b = self.down_tube.radius / dot(turn_right(d), vec)

			# Where the lower edge of DT touches HT
			corner = crown_ix - a*d + b*vec

			error = self.lower_ht_ext_edge - x - dot(crown_ix-corner, vec)
			if abs(error) < EPSILON: break
			x += error

		# We'll set tube_top and tube_bottom for the HT after computing the
		# Mitres.
		self.head_tube.place(top, crown_ix)
		self.down_tube.place(crown_ix, self.seat_tube.bottom)

	def _place_tubes(self):
		self._place_st()
		self._place_tt()
		self._place_ht_dt()

		self.bb_tube.place(array([0, 0, -self.bb_shell/2]),
				array([0, 0, self.bb_shell/2]))

	def _calc_dt_bb_mitre_corner(self):
		# The top of the down-tube is a line defined by a + l*e
		dt = self.down_tube
		a = dt.tube_top + turn_left(dt.vecn) * dt.get_radius(BOTTOM)
		e = -dt.vecn

		# Now solve for where that intersects the circle defined by the outside
		# edge of the BB shell (which we assume is centred at the origin,
		# because it is).
		r = self.bb_tube.radius
		sols = solve_quadratic(dot(e, e), 2*dot(a, e), dot(a, a) - r*r)

		# The line intersects the circle in two places, we want the one nearest
		# the head-tube, which is the shorter of the two.
		l = min(sols)
		inside = a + l * e

		# The outside corner is just the same on the other side of the tube
		outside = inside + turn_right(dt.vecn) * dt.get_diameter(BOTTOM)
		return inside, outside

	def _calc_mitres(self):
		# Put this into the XY plane so we can mitre it with the ST and DT. The
		# Mitre is the same (the builder will just rotate the template 90
		# degrees around the tube)
		rotated_bb = self.bb_tube.transform(array([[0, 0, 1],
			[0, 1, 0],
			[1, 0, 0]]))

		for datum in (
				("dt_ht", self.down_tube, TOP, self.head_tube, BOTTOM),
				("tt_ht", self.top_tube, TOP, self.head_tube, TOP),
				("st_bb", self.seat_tube, BOTTOM, rotated_bb, BOTTOM),
				("tt_st", self.top_tube, BOTTOM, self.seat_tube, TOP),
				("dt_bb", self.down_tube, BOTTOM, rotated_bb, BOTTOM),
				("dt_st", self.down_tube, BOTTOM, self.seat_tube, BOTTOM),
				):
			name, args = datum[0], datum[1:]
			mitre = Mitre(*args)
			setattr(self, name, mitre)

		self.head_tube.tube_top = (self.tt_ht.corners[1] +
				self.head_tube.vecn * self.upper_ht_ext_edge +
				turn_right(self.head_tube.vecn) * self.head_tube.get_radius(TOP))

		self.head_tube.tube_bottom = (self.dt_ht.corners[1] -
				self.head_tube.vecn * self.lower_ht_ext_edge +
				turn_right(self.head_tube.vecn) * self.head_tube.get_radius(BOTTOM))

		# Solve for the top of the DT/BB mitre as a special case (it would be
		# one of the Mitre corners, but our BB mitres are twisted because we're
		# only really working in 2D).
		self.dt_bb_mitre_corner = self._calc_dt_bb_mitre_corner()

	def calc_trail_and_wb(self):
		"""Work out the trail and the wheelbase. Easiest to do these
		together"""
		ground = self.bb_drop - (self.wheel_diameter / 2 + self.tyre_height)
		y = self.head_tube.tube_bottom[1]
		l = (ground - y) / self.head_tube.vecn[1]

		# The point where the projection of the head tube meets the ground
		p = self.head_tube.tube_bottom + l * self.head_tube.vecn

		# The centre of the front wheel
		fwc = self._fork_path()[-1]
		trail = p[0] - fwc[0]
		wb = fwc[0] - self.left_drop.end[0]
		return trail, wb

	def render_mitre_templates(self):
		print("""
Mitres
======""")

		for attr, caption, print_ref in (
				("tt_ht", "Top Tube to Head Tube", True),
				("tt_st", "Top Tube to Seat Tube", True),
				("dt_st", "Down Tube to Seat Tube", True),
				("dt_ht", "Down Tube to Head Tube", True),
				("st_bb", "Seat Tube to BB Shell", False),
				("dt_bb", "Down Tube to BB Shell", False),
				):
			mitre = getattr(self, attr)
			mitre.make_template(attr, caption, False)
			if print_ref:
				print(caption, "Mitres at ({}) from tube top".format(
						", ".join(["{:.2f}".format(x)
							for x in mitre.reference_points()])))

	def display(self):
		"""Show key metrics"""
		print("""
Front Triangle
==============

Tube Cuts
---------""")

		print("Head Tube cut square {:.2f}".format(self.head_tube.length()))

		l = abs(dot(self.tt_st.corners[OUTSIDE]
				- self.tt_ht.corners[OUTSIDE], self.top_tube.vecn))
		print("Top tube length between outside mitres: {:.2f}".format(l))

		l = abs(dot(self.tt_st.corners[INSIDE]
				- self.tt_ht.corners[INSIDE], self.top_tube.vecn))
		print("Top tube length between inside mitres: {:.2f}".format(l))

		print("Seat Tube length from top to BB mitre: {:.2f}".format(
			self.seat_tube.length() - self.bb_tube.get_radius(BOTTOM)))

		a = abs(dot(self.dt_ht.corners[INSIDE] - self.dt_st.corners[INSIDE],
			self.down_tube.vecn))
		print("Down Tube length from inside " \
				"ST mitre to inside HT mitre: {:.2f}".format(a))

		b = abs(dot(self.dt_ht.corners[INSIDE] - self.dt_bb_mitre_corner[INSIDE],
			self.down_tube.vecn))
		print("Down Tube length from inside " \
				"BB mitre to inside ST mitre: {:.2f}".format(b))

		print("Offset between DT/BB mitre and DT/ST mitre on " \
				"DT inside centreline: {:.2f}".format(b - a))

		l = abs(dot(self.dt_ht.corners[OUTSIDE] - self.dt_bb_mitre_corner[OUTSIDE],
			self.down_tube.vecn))
		print("Down Tube length from outside BB mitre to " \
				"outside HT mitre: {:.2f}".format(l))

		print("""
Other Metrics
-------------""")

		reach_stack = self.head_tube.tube_top
		print("Stack: {:.2f}".format(reach_stack[1]))
		print("Reach: {:.2f}".format(reach_stack[0]))

		trail, wb = self.calc_trail_and_wb()
		print("Trail: {:.2f}".format(trail))
		print("Wheelbase: {:.2f}".format(wb))

		print("Angle between ST and DT: {:.2f}deg".format(
			rad2deg(self.dt_st.angle)))

		print("Angle between DT and HT: {:.2f}deg".format(
			rad2deg(pi - self.dt_ht.angle)))

		print("Angle between ST and TT: {:.2f}deg".format(
			rad2deg(pi - self.tt_st.angle)))

		print("Angle between TT and HT: {:.2f}deg".format(
			rad2deg(self.tt_ht.angle)))

		self.rear_triangle.display()

	def _lower_stack_path(self):
		v = self.head_tube.vecn
		w = turn_left(v)
		b = self.head_tube.tube_bottom
		r = self.head_tube.get_radius(BOTTOM)

		start = b - v * self.lower_stack + w * r

		long_side = -w * r * 2
		short_side = v * self.lower_stack

		return [start, start + long_side, start + long_side + short_side,
				start + short_side, start]

	def _fork_path(self):
		v = -self.head_tube.vecn
		w = turn_left(v)

		start = self.head_tube.tube_bottom + self.lower_stack * v
		middle = start + self.fork_len * v
		end = middle + w * self.fork_offset

		return [start, end]

	def render_scale_diagram(self, outfile, px_per_mm=2):
		fork_path = self._fork_path()

		inf = array([self.left_cs.tube_top[0], 0])
		sup = array([fork_path[-1][0], self.head_tube.tube_top[1]])

		r = Render(inf, sup, px_per_mm)

		for tube in (self.seat_tube, self.top_tube, self.head_tube,
				self.down_tube):
			r.draw_tube(tube)

		for mitre in (self.tt_ht, self.dt_ht, self.dt_st, self.tt_st):
			r.draw_mitre(mitre)

		# Draw the lower stack
		path = self._lower_stack_path()
		r.polyline(self._lower_stack_path(), "black")

		# The fork
		r.polyline(fork_path, "darkblue")

		# And the centreline
		start = array([self.right_drop.end[0], self.bb_drop])
		end = array([fork_path[-1][0], self.bb_drop])
		r.polyline((start, end), "black")

		# The chainstays
		r.draw_tube(self.right_cs)
		r.draw_tube(self.right_ss)

		# The CS dropout
		do_top = self.right_drop.end + \
				self.right_ss.vecn * self.right_drop.ss_length
		r.polyline((self.right_drop.start, self.right_drop.end,
			do_top), "purple")

		# The wheels
		r.draw_wheel(self.right_cs.top + self.right_drop.vec,
				self.wheel_diameter, self.tyre_height)
		r.draw_wheel(fork_path[1], self.wheel_diameter, self.tyre_height)

		# And finally the BB shell
		r.circle(array([0, 0]), self.bb_tube.radius, "brown")
		r.draw_point(self.dt_bb_mitre_corner[INSIDE], "red")
		r.draw_point(self.dt_bb_mitre_corner[OUTSIDE], "cyan")

		r.save(outfile)
		print("Rendered diagram to {}".format(outfile))

def main():
	ap = ArgumentParser()
	ap.add_argument("-c", "--config", type=str, required=True)
	ap.add_argument("-r", "--resolution", type=int, default=2)

	args = ap.parse_args()
	config = configparser.ConfigParser()
	config.read(args.config)
	frame = Frame(config)
	frame.render_mitre_templates()
	frame.render_scale_diagram("side_view.png", args.resolution)
	frame.display()

	frame.rear_triangle.render_scale_diagram(
			"chainstays.png", args.resolution)

if __name__ == "__main__":
	main()
