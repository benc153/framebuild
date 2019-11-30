What's this?
============

There are three programs you can run here-- `framebuild.py`, `copecalc.py` and
`riderfit.py`.

## How to run them

You need [Python](https://www.python.org/), specifically Python 3 (any version
that starts with a 3 is fine, e.g. 3.7.1 is good) which is freely available for
all widely used operating systems. You also need a rough idea how to run
programs from the "command line" or "command prompt". You will also need to be
able to open and print out `.png` files if you want to use the mitre templates. 

## copecalc.py

This generates mitre templates, similar to
[this](http://www.metalgeek.com/static/cope.pcgi) online tool.

The idea is you generate a template for the join you're trying to make, print
it out, cut it out with scissors, wrap it around the tube and trace the curved
edge onto the tube with a pen or something. Then cut/sand/file the tube to that
line.

This version can join round tubes to round tubes and also elliptical to round.
The elliptical tube can come out at any angle but it has to have either its
semi-major axis or its semi-minor axis aligned with the parent tube. This is
all you need if what you're making is fairly symmetrical. This can be useful
for chainstays if they're elliptical (or just generally squashed, which is
close enough to ellipitcal) at the front.

It can also do offsets-- if you want the cut tube not to aim to the centre of
the parent tube but off to one side. This is useful for seatstays.

It can also generate a template for a tapered parent tube, assuming a linear
taper. You specify the ratio of diameter change to length change.

You can also adjust the resolution at which it prints. It defaults to 100px per
inch which is the same as the one on metalgeek. I print from [the
GIMP](https://www.gimp.org/) which has an Image Settings tab where you can set
the pixels per inch.

## framebuild.py

This takes a configuration file like `default.ini`, in which you specify the
dimensions of all your tubes and your desired geometry.

It outputs a couple of diagrams, a bunch of mitre templates, all the
measurements you need to cut the tubes to the right lengths, and a few other
useful metrics.

## riderfit.py

This takes a rider configuration file like `default_rider.ini` and optionally
also a frame configation file. If you provide a frame configuration it just
means it can draw a diagram of the rider actually sitting on the bike.

To get a sense of what size bike you want for a particular rider first put her
measurements, and those of the handlebar and stem and the seat tube angle into
the configuration file. Then decide on how much bar drop and what sort of
shoulder angle you want based on how aggressive a position you want. The
program will output the stack and reach you want to aim for. The shoulder angle
is the angle between the rider's upper arm and his torso, and should be a
maximum of 80 degrees for the most aggressive position.

### Designing your bike

The following is all assuming a lugless frame-- i.e. welded or fillet brazed.
If you're making a lugged frame, I guess the lugs you have already define the
angles, and you need to solve for the lengths. This program does output the
angles, so maybe you could tune the lengths to get the angles you need.

#### Choose some numbers

1. Pick a standard wheel size.
2. Decide what stack and reach you want (possibly using `riderfit.py`) -- these
   are the two measurements that depend the most on the size and proportions of
   the rider and what sort of riding position you're looking for.
3. Decide on a seat angle-- usually 73 or 74 degrees.
4. Decide on a head angle based on your understanding of the black art of
   steering geometry and what fork you're going to be using.
5. Pick a chainstay length taking into consideration whether you want to fit
   mudguards and things.
6. Decide on a BB drop. Most bikes use 70mm. Reduce this if you want more
   ground clearance, maybe for a cross bike. Increase it for a lower centre of
   gravity, or if you want it to be easier to put your feet on the ground.
   Don't go too low (too much drop) for a fixed-gear bike in case of pedal
   strike around corners.
7. The `[Extensions]` section is based on safe margins so you aren't welding
   too close to the ends of the tubes, have room for a seat-clamp etc. Probably
   better not to reduce them by much but you can make them bigger. I've seen
   Surly touring bikes with an increased head tube upper extension presumably
   so they could get a higher stack while keeping the top tube horizontal.
8. Think about where you stand on aesthetic issues-- do you want a horizontal
   top tube? Maybe a top tube and seat tube the same length for a classic
   Italian "square" look? Do you want to avoid a really short or really long
   head-tube?

#### Fine-tune them

1. Copy default.ini to mybike.ini
2. Set up any numbers you have pretty much decided on in the Basic Dimensions
   and Angles sections.
2. Run `python3 framebuild.py -c mybike.ini`
3. Look at the "Other Metrics" section and see what stack and reach you ended
   up with.
4. Adjust basic dimensions and repeat until you get the stack and reach you
   want while maximizing your aesthetic considerations. To increase reach, grow
   the top tube length. To increase stack, increase the top tube angle or the
   seat tube length.
5. Check `side_view.png` to see what your bike's looking like and to check the
   wheels actually fit etc.
6. Check `chainstays.png` for tyre and chainring clearance. Note that it will
   look tighter in the diagram than it really is because the chainstays are
   drawn as a constant width tube. In reality they will have an oval section
   about where the tyre hits.
7. Repeat steps 2 to 6.

Note that `side_view.png` is drawn completely to scale, but is projected flat
onto the paper. The main triangle is in the plane of the paper anyway, so if
you measure pixels there, or even print the diagram out and measure it with a
ruler, it should match up perfectly. But the rear triangle measurements will
appear a bit shorter than the lengths you're going to actually cut things to
because of the projection. The lengths output by the program are the right
ones.

`chainstays.png` is drawn in a reference frame in which the chainstays are flat
on the paper, so that should have the right lengths.

### Tapered Head-Tubes

`framebuild` doesn't support tapered head-tubes but `copecalc` does. You can
work around this as follows:

1. Design your frame using the diameter that the HT has at the top.
2. Once you've cut the actual tapered HT to length, mark the positions of the
   DT mitre on it.
3. Measure the diameter of the HT at the top and bottom of the mitre, and work
   out the diameter change per unit length.
4. Run `copecalc` with the -t option to make your template.
5. Replace the HT diameter in your ini file with the diameter at the top of the
   mitre. Run `framebuild` again and this will give you the correct DT length
   from the inside of the ST mitre to the inside of the HT mitre.
6. Replace the HT diameter with the diameter at the top of the mitre, run
   `framebuild` again, and this will give you the correct DT length from the
   outside of the BB mitre to the outside of the HT mitre. You don't strictly
   need both these measurements-- the one in step 5 is enough-- but it's good
   to double-check.

When you change the HT length and rerun `framebuild` the angle of the DT will
actually change very slightly because it's solving to keep your lower extension
at the value you asked for. But the difference is not significant and the cut
lengths you get from changing that diameter will be close enough.

### Flexural Rigidity and Deflection Estimates

For each tube the program works out the [Flexural
Rigidity](https://en.wikipedia.org/wiki/Flexural_rigidity) and an estimate of
the deflection of each tube under a 100kg load at its centre. For these
calculations it isn't sophisticated enough to understand butted tubes properly
so it works with the thinner wall thickness in the middle.

It's not that you expect to get that load in particular in normal use but it
provides an interesting point of comparison. You probably want the deflection
estimate to be similar for all frames-- around 3 or 4mm seems to be about right
for a nice lively feeling steel frame.

For example, a long low slack MTB with a 750mm down tube deflects around 4.5mm
under this theoretical test load in spite of having a diameter of 35mm,
compared to 4.08mm for a road frame with a 620mm DT which only has a 28.6mm
diameter.

### Building it

Once you're happy with the design, you can use the mitre templates and the
distances between mitres to cut up the tubes. Note that the down-tube has two
mitres at the bottom end-- where it joins the BB shell and where it joins the
seat-tube. The program gives you the offset between the two mitres, so you just
trace both templates onto the tube with that offset between then.

The first thing to do with the tubes is mark a centre-line, for example by
putting the tube inside a piece of angle iron and using the edge as a ruler.
You will match up the centre-line with the guidelines on the mitre templates to
make sure that when you have mitres at both ends of a tube that they aren't
rotated with respect to each other.

The program also gives you reference points to mark on the tubes, measured from
the top, showing where the top and bottom edges of the joining tubes go. Mark
these on the tubes with a pen and then when you're fitting things up you can
put everything in the right place.

Note that double-butted tubes are supplied with a longer butt on one end. The
idea is you cut that end off to get the length you need. So you will want to
put one mitre right at the short end (which has a bit of red paint on it in the
case of Columbus tubing), and then measure from there to where the other mitre
needs to go.

The program outputs mitre templates for the seat-stays on the assumption that
they're offset from the centre of the seat-tube by their diameter. Otherwise
they would bump into each other. With the templates generated from the program
and typical seat-tube and seat-stay diameters the tops of the two seat-stays
will fit next to each other very well.

Chain-stays usually have a flattened bit in them where the tyre goes. If you're
using 700C wheels this probably means you want to cut the chainstays to length
from the front since the flattened bit will be in the right place for that
wheel size.

Before cutting anything for the rear triangle, decide how you're going to
attach the dropouts you have, measure them, and set `cs_length` and `ss_length`
in the `[Dropouts]` section. Run the program with the new data to get the final
lengths.

## Why not just use BikeCad?

[BikeCad](https://www.bikecad.ca/) looks pretty awesome. I started with the
free version, but either it didn't give me all the dimensions I actually needed
to build the bike, or I couldn't find them. I also wanted to double-check
things anyway before actually cutting into my shiny tubes.

I could have bought the pro version at this point... But then I could have just
bought a frame, or a whole bike, and where would be the fun in that?
