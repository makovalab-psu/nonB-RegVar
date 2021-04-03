#!/usr/bin/env python
"""
Sub-annotate g-quadruplex motifs.
"""

from sys    import argv,stdin,stderr,exit
from string import maketrans
from re     import compile as re_compile


programName    = "gee_kwad"
programVersion = "0.1.2"


def usage(s=None):
	message = """
Parse g-quadruplex motif annotations into sub-annotations.

usage: cat bed_file | %s [options]
  --allow:bulges       allow bulges in stems
  --disallow:bulges    don't allow bulges in stems
                       (this is the default)
  --allow:gloop        allow loops consisting of a single G; this allows
                       a long stem to be split into stem-loop-stem
                       (this is the default)
  --disallow:gloop     don't allow loops consisting of a single G
  --warn:longloops     warn the user about sequences with long loops
  --nowarn:longloops   don't warn the user about sequences with long loops
                       (this is the default)
  --warn:bulges        warn the user about sequences with bulged stems
  --nowarn:bulges      don't warn the user about sequences with bulged stems
                       (this is the default)
  --warn:tail          warn the user about sequences with tails
                       (this is the default)
  --nowarn:tail        don't warn the user about sequences with tails
  --parse=fourstems    parse with a preference for exactly four stems
  --copyinput          copy the input lines to the output, as comments with
                       a "#" prefix
  --head=<number>      limit the number of input lines
  --progress=<number>  periodically report how many lines we've read
  --version            show version number and quit

The <bed_file> contains lines that look like this:
  chr1 11008 11026 GGGCGGGGGTTGGGGGGG                       18 + 34.98
  chr1 11058 11078 GGGCTGGGGCGGGGGGAGGG                     20 + 33.43
  chr1 12959 12999 CCCCTTCACTCCCAGCTCAGAGCCCAGGCCAGGGGCCCCC 40 - 3.37
  chr1 14236 14273 CCCCTGCTGACTGCCCTTCTCTCCTCCCTCTCATCCC    37 - 3.26
  chr1 15096 15129 CCCTGTCCCCACCCCCATGACACTCCCCAGCCC        33 - 33.58
The 5th, 6th, and 7th columns are optional. If the 6th column, the strand, is
absent, we try parsing for either strand. The 5th and 7th columns are ignored
in any case.

Note that this program does NOT search genomes for g-quadruplexes. It is
assumed that some other program has done this. What this program does is
break those g-quadruplexes into stems and loops.""" \
% programName

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():
	global allowBulges,allowGLoops,allowBadLength
	global debug

    # parse the command line

	allowBulges     = False
	allowGLoops     = True
	allowBadLength  = False
	warnOnLongLoops = False
	warnOnBulges    = False
	warnOnTails     = True
	parseAs         = "latest version"
	copyInputLines  = False
	headLimit       = None
	reportProgress  = None
	debug           = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg in ["--allow:bulges","--allowbulges","--allow_bulges",
		            "--allow:bulge", "--allowbulge", "--allow_bulge"]):
			allowBulges = True
		elif (arg in ["--disallow:bulges","--disallowbulges","--disallow_bulges",
		              "--disallow:bulge", "--disallowbulge", "--disallow_bulge"]):
			allowBulges = False
		elif (arg in ["--allow:gloops","--allowgloops","--allow_gloops",
		              "--allow:gloop", "--allowgloop", "--allow_gloop",
		              "--allow:Gloops","--allowGloops","--allow_Gloops",
		              "--allow:Gloop", "--allowGloop", "--allow_Gloop"]):
			allowGLoops = True
		elif (arg in ["--disallow:gloops","--disallowgloops","--disallow_gloops",
		              "--disallow:gloop", "--disallowgloop", "--disallow_gloop",
		              "--disallow:Gloops","--disallowGloops","--disallow_Gloops",
		              "--disallow:Gloop", "--disallowGloop", "--disallow_Gloop"]):
			allowGLoops = False
		elif (arg in ["--allow:badlengths","--allow:bad_lengths","--allowbadlengths","--allow_bad_lengths",
		              "--allow:badlength", "--allow:bad_length", "--allowbadlength", "--allow_bad_length"]):
			allowBadLength = True
		elif (arg in ["--warn:longloops","--warn:long_loops",
		              "--warn:longloop", "--warn:long_loop"]):
			warnOnLongLoops = True
		elif (arg in ["--nowarn:longloops","--nowarn:long_loops",
		              "--nowarn:longloop", "--nowarn:long_loop"]):
			warnOnLongLoops = False
		elif (arg in ["--warn:bulges","--warn:bulge"]):
			warnOnBulges = True
		elif (arg in ["--nowarn:bulges","--nowarn:bulge"]):
			warnOnBulges = False
		elif (arg in ["--warn:tails","--warn:tail"]):
			warnOnTails = True
		elif (arg in ["--nowarn:tails","--nowarn:tail"]):
			warnOnTails = False
		elif (arg in ["--parse=fourstems","--parse=4stems"]):
			parseAs = "4 stems"
		elif (arg in ["--copyinput","--copylines"]):
			copyInputLines = True
		elif (arg.startswith("--head=")):
			headLimit = int_with_unit(argVal)
		elif (arg.startswith("--progress=")):
			reportProgress = int_with_unit(argVal)
		elif (arg in ["--version","--v","--V","-version","-v","-V"]):
			exit("%s, version %s" % (programName,programVersion))
		elif (arg == "--debug"):
			debug += ["debug"]
		elif (arg.startswith("--debug=")):
			debug += argVal.split(",")
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			usage("unrecognized option: %s" % arg)

	# process the putative g-quadruplex motifs

	gQuadParser = parse_as_g_quad
	if (parseAs == "4 stems"): gQuadParser = parse_as_g_quad_4_stems

	itemNum = 0
	for g4 in read_gquad_bed(stdin):
		itemNum += 1
		if (headLimit != None) and (itemNum > headLimit):
			print >>stderr, "limit of %s items reached" % (commatize(headLimit))
			break
		if (reportProgress != None) and (itemNum % reportProgress == 0):
			print >>stderr, "progress: item %s (%s %d %d)" \
			              % (commatize(itemNum),g4.chrom,g4.start,g4.end)

		# parse the motif

		strand = g4.strand
		if (strand == "+"):
			parts = gQuadParser(g4.motifSeq)
		elif (strand == "-"):
			motifRev = reverse_complement(g4.motifSeq)
			parts = gQuadParser(motifRev)
		else: # if (strand == None):
			parts = gQuadParser(g4.motifSeq)
			if (parts != None):
				strand = "+"
			else:
				motifRev = reverse_complement(g4.motifSeq)
				parts = gQuadParser(motifRev)
				if (parts != None):
					strand = "-"

		# report any warnings to the user and/or to the output

		if (warnOnBulges) and (parts != None) and (parts.hasBulge):
			print >>stderr, "WARNING: bulge in %s %d %d: %s" \
			              % (g4.chrom,g4.start,g4.end,g4.motifSeq)

		if (warnOnLongLoops) and (parts != None) and (parts.hasLongLoop):
			print >>stderr, "WARNING: long loop in %s %d %d: %s" \
			              % (g4.chrom,g4.start,g4.end,g4.motifSeq)

		if (warnOnTails) and (parts != None) and (parts.tail != None):
			print >>stderr, "WARNING: tail in %s %d %d: %s" \
			              % (g4.chrom,g4.start,g4.end,g4.motifSeq)

		if (copyInputLines):
			print "# %s" % g4.line

		stemLoopInconsistency = (parts != None) and (len(parts.stem) != len(parts.loop)+1)
		if (parts == None) or (stemLoopInconsistency):
			if (stemLoopInconsistency):
				message = "sub-annotation problem: %d stems and %d loops" % (len(parts.stem),len(parts.loop))
			else:
				message = "unable to sub-annotate"

			if (copyInputLines):
				print "# (%s)" % message
			else:
				print "# %s %s" % (message,g4.line)

			if (strand != None):
				print >>stderr, "WARNING: unable to sub-annotate %s %d %d %s: %s" \
				              % (g4.chrom,g4.start,g4.end,g4.strand,g4.motifSeq)
			else:
				print >>stderr, "WARNING: unable to sub-annotate %s %d %d: %s" \
				              % (g4.chrom,g4.start,g4.end,g4.motifSeq)
			continue

		# output the sub-annotations

		numStems = len(parts.stem)

		if (strand == "+"):
			stemStart = g4.start
			for stemIx in xrange(numStems):
				stemEnd = stemStart + len(parts.stem[stemIx])
				print "%s\t%d\t%d\t%s\t%d\t%s\tstem%d" \
				    % (g4.chrom,stemStart,stemEnd,
				       parts.stem[stemIx],len(parts.stem[stemIx]),
				       strand,1+stemIx)

				if (stemIx == numStems-1): break

				loopStart = stemEnd
				loopEnd   = loopStart + len(parts.loop[stemIx])
				print "%s\t%d\t%d\t%s\t%d\t%s\tloop%d" \
				    % (g4.chrom,loopStart,loopEnd,
				       parts.loop[stemIx],len(parts.loop[stemIx]),
				       strand,1+stemIx)

				stemStart = loopEnd

			if (parts.tail == None):
				assert (stemEnd == g4.end)
			else:
				tailStart = stemEnd
				tailEnd   = tailStart + len(parts.tail)
				assert (tailEnd == g4.end)

			if (parts.tail != None):
				print "%s\t%d\t%d\t%s\t%d\t%s\ttail" \
				    % (g4.chrom,tailStart,tailEnd,
				       parts.tail,len(parts.tail),
				       strand)

		else: # if (strand == "-"):
			stemStart = g4.end
			for stemIx in xrange(numStems):
				stemEnd = stemStart - len(parts.stem[stemIx])
				print "%s\t%d\t%d\t%s\t%d\t%s\tstem%d" \
				    % (g4.chrom,stemEnd,stemStart,
				       reverse_complement(parts.stem[stemIx]),len(parts.stem[stemIx]),
				       strand,1+stemIx)

				if (stemIx == numStems-1): break

				loopStart = stemEnd
				loopEnd   = loopStart - len(parts.loop[stemIx])
				print "%s\t%d\t%d\t%s\t%d\t%s\tloop%d" \
				    % (g4.chrom,loopEnd,loopStart,
				       reverse_complement(parts.loop[stemIx]),len(parts.loop[stemIx]),
				       strand,1+stemIx)

				stemStart = loopEnd

			if (parts.tail == None):
				assert (stemEnd == g4.start)
			else:
				tailStart = stemEnd
				tailEnd   = tailStart - len(parts.tail)
				assert (tailEnd == g4.start)

			if (parts.tail != None):
				print "%s\t%d\t%d\t%s\t%d\t%s\ttail" \
				    % (g4.chrom,tailEnd,tailStart,
				       reverse_complement(parts.tail),len(parts.tail),
				       strand)


# parse_as_g_quad--
#	Try to parse a sequence, in its entirety, as a g-qudruplex motif.  If
#	succesful, return an object describing the parts of the motif. Otherwise,
#	return None.
#
# Shoutout to
#	https://docs.python.org/3/howto/regex.html#greedy-versus-non-greedy
# for how to make regular expressions that give the shortest match instead of
# the longest

reNt         = "[ACGTNacgtn]"
reNonG       = "[ACTNactn]"
reStem       = "[Gg]{3,}"
reBulgedStem = "[Gg]("+reNonG+"{0,2}[Gg]){2,}"
reLoop       = reNt+"{0,6}"+reNonG  # this used to be reNt+"{1,7}"
reLongLoop   = reNt+"{0,}"+reNonG   # this used to be reNt+"{1,}"
reTail       = reNt+"*"


gQuad43Full \
    = re_compile("^"
               + "(?P<stem1>"+reStem+")(?P<loop1>"+reLoop+")"
               + "(?P<stem2>"+reStem+")(?P<loop2>"+reLoop+")"
               + "(?P<stem3>"+reStem+")(?P<loop3>"+reLoop+")"
               + "(?P<stem4>"+reStem+")(?P<tail>" +reTail+")"
               + "$")

gQuad43BulgesFull \
    = re_compile("^"
               + "(?P<stem1>"+reBulgedStem+")(?P<loop1>"+reLoop+")"
               + "(?P<stem2>"+reBulgedStem+")(?P<loop2>"+reLoop+")"
               + "(?P<stem3>"+reBulgedStem+")(?P<loop3>"+reLoop+")"
               + "(?P<stem4>"+reBulgedStem+")(?P<tail>" +reTail+")"
               + "$")

gQuad43LongLoopFull \
    = re_compile("^"
               + "(?P<stem1>"+reStem+")(?P<loop1>"+reLongLoop+")"
               + "(?P<stem2>"+reStem+")(?P<loop2>"+reLongLoop+")"
               + "(?P<stem3>"+reStem+")(?P<loop3>"+reLongLoop+")"
               + "(?P<stem4>"+reStem+")(?P<tail>" +reTail+")"
               + "$")

gQuad43BulgesLongLoopFull \
    = re_compile("^"
               + "(?P<stem1>"+reBulgedStem+")(?P<loop1>"+reLongLoop+")"
               + "(?P<stem2>"+reBulgedStem+")(?P<loop2>"+reLongLoop+")"
               + "(?P<stem3>"+reBulgedStem+")(?P<loop3>"+reLongLoop+")"
               + "(?P<stem4>"+reBulgedStem+")(?P<tail>" +reTail+")"
               + "$")

gQuad32Full \
    = re_compile("^"
               + "(?P<stem1>"+reStem+")(?P<loop1>"+reLoop+")"
               + "(?P<stem2>"+reStem+")(?P<loop2>"+reLoop+")"
               + "(?P<stem3>"+reStem+")(?P<tail>" +reTail+")"
               + "$")

gQuad32BulgesFull \
    = re_compile("^"
               + "(?P<stem1>"+reBulgedStem+")(?P<loop1>"+reLoop+")"
               + "(?P<stem2>"+reBulgedStem+")(?P<loop2>"+reLoop+")"
               + "(?P<stem3>"+reBulgedStem+")(?P<tail>" +reTail+")"
               + "$")

gQuad32LongLoopFull \
    = re_compile("^"
               + "(?P<stem1>"+reStem+")(?P<loop1>"+reLongLoop+")"
               + "(?P<stem2>"+reStem+")(?P<loop2>"+reLongLoop+")"
               + "(?P<stem3>"+reStem+")(?P<tail>" +reTail+")"
               + "$")

gQuad32BulgesLongLoopFull \
    = re_compile("^"
               + "(?P<stem1>"+reBulgedStem+")(?P<loop1>"+reLongLoop+")"
               + "(?P<stem2>"+reBulgedStem+")(?P<loop2>"+reLongLoop+")"
               + "(?P<stem3>"+reBulgedStem+")(?P<tail>" +reTail+")"
               + "$")

gQuad21Full \
    = re_compile("^"
               + "(?P<stem1>"+reStem+")(?P<loop1>"+reLoop+")"
               + "(?P<stem2>"+reStem+")(?P<tail>" +reTail+")"
               + "$")

gQuad21BulgesFull \
    = re_compile("^"
               + "(?P<stem1>"+reBulgedStem+")(?P<loop1>"+reLoop+")"
               + "(?P<stem2>"+reBulgedStem+")(?P<tail>" +reTail+")"
               + "$")

gQuad21LongLoopFull \
    = re_compile("^"
               + "(?P<stem1>"+reStem+")(?P<loop1>"+reLongLoop+")"
               + "(?P<stem2>"+reStem+")(?P<tail>" +reTail+")"
               + "$")

gQuad21BulgesLongLoopFull \
    = re_compile("^"
               + "(?P<stem1>"+reBulgedStem+")(?P<loop1>"+reLongLoop+")"
               + "(?P<stem2>"+reBulgedStem+")(?P<tail>" +reTail+")"
               + "$")

gQuad10Full \
    = re_compile("^"
               + "(?P<stem1>"+reStem+")(?P<tail>" +reTail+")"
               + "$")

gQuad10BulgesFull \
    = re_compile("^"
               + "(?P<stem1>"+reBulgedStem+")(?P<tail>" +reTail+")"
               + "$")

gQuad10LongLoopFull \
    = re_compile("^"
               + "(?P<stem1>"+reStem+")(?P<tail>" +reTail+")"
               + "$")

gQuad10BulgesLongLoopFull \
    = re_compile("^"
               + "(?P<stem1>"+reBulgedStem+")(?P<tail>" +reTail+")"
               + "$")

loopStemLongLoopFull \
    = re_compile("^"
               + "(?P<loop1>"+reLoop+")"
               + "(?P<stem1>"+reStem+")"
               + "(?P<loop2>"+reLongLoop+")"
               + "$")

loopStemLongLoopBulgesFull \
    = re_compile("^"
               + "(?P<loop1>"+reLoop+")"
               + "(?P<stem1>"+reBulgedStem+")"
               + "(?P<loop2>"+reLongLoop+")"
               + "$")

longLoopStemLongLoopFull \
    = re_compile("^"
               + "(?P<loop1>"+reLongLoop+")"
               + "(?P<stem1>"+reStem+")"
               + "(?P<loop2>"+reLongLoop+")"
               + "$")

longLoopStemLongLoopBulgesFull \
    = re_compile("^"
               + "(?P<loop1>"+reLongLoop+")"
               + "(?P<stem1>"+reBulgedStem+")"
               + "(?P<loop2>"+reLongLoop+")"
               + "$")


loopAndStemShortest \
    = re_compile("^(?P<loop>"+reLoop+"?)(?P<stem>"+reStem+"?)")

loopAndStemBulgesShortest \
    = re_compile("^(?P<loop>"+reLoop+"?)(?P<stem>"+reBulgedStem+"?)")

loopAndStemLongLoopShortest \
    = re_compile("^(?P<loop>"+reLongLoop+"?)(?P<stem>"+reStem+"?)")

loopAndStemBulgesLongLoopShortest \
    = re_compile("^(?P<loop>"+reLongLoop+"?)(?P<stem>"+reBulgedStem+"?)")


stemLongest \
    = re_compile("^(?P<stem>"+reStem+")")

stemBulgesLongest \
    = re_compile("^(?P<stem>"+reBulgedStem+")")


class GQuadParts: pass

def parse_as_g_quad(seq):

	# try to parse it with preference for exactly four stems

	parts = parse_as_g_quad_4_stems(seq)
	if (allowGLoops):
		if (parts == None):
			parts = parse_as_g_quad_3_stems(seq)
		if (parts == None):
			parts = parse_as_g_quad_2_stems(seq)
		if (parts == None):
			parts = parse_as_g_quad_1_stem(seq)

	if (parts == None):
		return None

	# try to re-parse each loop into a loop-stem-loop
	#
	# note that we don't use the typical "for ix in xrange(parts.loop)" loop
	# because we are potentially modifying the lists as we go
	#
	# also note that we can change an original loop into a series of
	# loop-stem-loop-...-stem-loop of any length; this sorta happens right to
	# left, since the regex finds the rightmost stem

	ix = 0
	while (ix < len(parts.loop)):
		loop = parts.loop[ix]

		m = loopStemLongLoopFull.match(loop)
		if (m == None):
			m = longLoopStemLongLoopFull.match(loop)

		if (allowBulges):
			if (m == None):
				m = loopStemLongLoopBulgesFull.match(loop)
				if (m != None): parts.hasBulge = True
			if (m == None):
				m = longLoopStemLongLoopBulgesFull.match(loop)
				if (m != None): parts.hasBulge = True

		if (m != None):
			loop1 = m.group("loop1")
			stem1 = m.group("stem1")
			loop2 = m.group("loop2")

			if ("loop-stem-loop" in debug):
				print >>stderr, "%s becomes %s/%s/%s" \
				              % (loop,loop1,stem1,loop2)

			parts.loop[ix] = loop1  # replaces loop
			parts.stem.insert(ix+1,stem1)
			parts.loop.insert(ix+1,loop2)
			# we *don't* increment ix in this case, because we may still have
			# a stem embedded in loop1

		else:
			ix += 1

	# if we don't have at least four stems, try to split stems into stem-G-stem,
	# using the middle G as a loop
	#
	# note that we don't use the typical "for ix in xrange(parts.stem)" stem
	# because we are potentially modifying the lists as we go


	if (allowGLoops) and (len(parts.stem) < 4):
		ix = 0
		while (ix < len(parts.stem)) and (len(parts.stem) < 4):
			stem = parts.stem[ix]
			if (len(stem) < 7):
				ix += 1
				continue

			isAllGs = True
			for nuc in stem:
				if (nuc != "G"):
					isAllGs = False
					break
			if (not isAllGs):
				ix += 1
				continue

			stem1 = stem[:3]
			loop1 = stem[3]
			stem2 = stem[4:]

			if ("stem-loop-stem" in debug):
				print >>stderr, "%s becomes %s/%s/%s" \
				              % (loop,stem1,loop1,stem2)

			parts.stem[ix] = stem1  # replaces loop
			parts.loop.insert(ix  ,loop1)
			parts.stem.insert(ix+1,stem2)
			ix += 1

	# if we still don't have at least four stems, give up

	if (len(parts.stem) < 4):
		return None

	# check whether the result has any long loops; it might have originally
	# had some, but they could have been shortened to loop-stem-loop

	parts.hasLongLoop = False
	for loop in parts.loop:
		if (len(loop) > 7):
			parts.hasLongLoop = True
			break

	return parts


# parse_as_g_quad_4_stems--
#	Gives preference to parsing the string into a four-stem object.

def parse_as_g_quad_4_stems(seq):
	if ("regex" in debug):
		print >>stderr, "seq = \"%s\"" % seq

	hasLongLoop = hasBulge = False
	m = gQuad43Full.match(seq)

	if (m == None):
		m = gQuad43LongLoopFull.match(seq)
		if (m != None): hasLongLoop = True

	if (allowBulges):
		if (m == None):
			m = gQuad43BulgesFull.match(seq)
			if (m != None): hasBulge = True
		if (m == None):
			m = gQuad43BulgesLongLoopFull.match(seq)
			if (m != None): hasLongLoop = hasBulge = True

	if (m == None):     # shouldn't happen, assuming the inputs are really
		return None     # .. g-quadruplex motifs

	parts = GQuadParts()
	parts.hasLongLoop = hasLongLoop
	parts.hasBulge    = hasBulge
	parts.stem  = []
	parts.loop  = []
	parts.stem += [m.group("stem1")]
	parts.loop += [m.group("loop1")]
	parts.stem += [m.group("stem2")]
	parts.loop += [m.group("loop2")]
	parts.stem += [m.group("stem3")]
	parts.loop += [m.group("loop3")]
	parts.stem += [m.group("stem4")]
	leftover   =   m.group("tail")

	if ("regex" in debug):
		print >>stderr, "  stem1: \"%s\"" % parts.stem1
		print >>stderr, "  loop1: \"%s\"" % parts.loop1
		print >>stderr, "  stem2: \"%s\"" % parts.stem2
		print >>stderr, "  loop2: \"%s\"" % parts.loop2
		print >>stderr, "  stem3: \"%s\"" % parts.stem3
		print >>stderr, "  loop3: \"%s\"" % parts.loop3
		print >>stderr, "  stem4: \"%s\"" % parts.stem4
		print >>stderr, "  tail:  \"%s\"" % leftover

	# if there's any left over, try to re-parse it into more loops and stems

	if (leftover != ""):
		leftover = reparse_leftover(parts,leftover)

	parts.tail = None if (leftover == "") else leftover

	return parts


# parse_as_g_quad_3_stems--
#	Gives preference to parsing the string into a three-stem object.

def parse_as_g_quad_3_stems(seq):
	if ("regex" in debug):
		print >>stderr, "seq = \"%s\"" % seq

	hasLongLoop = hasBulge = False
	m = gQuad32Full.match(seq)

	if (m == None):
		m = gQuad32LongLoopFull.match(seq)
		if (m != None): hasLongLoop = True

	if (allowBulges):
		if (m == None):
			m = gQuad32BulgesFull.match(seq)
			if (m != None): hasBulge = True
		if (m == None):
			m = gQuad32BulgesLongLoopFull.match(seq)
			if (m != None): hasLongLoop = hasBulge = True

	if (m == None):     # shouldn't happen, assuming the inputs are really
		return None     # .. g-quadruplex motifs

	parts = GQuadParts()
	parts.hasLongLoop = hasLongLoop
	parts.hasBulge    = hasBulge
	parts.stem  = []
	parts.loop  = []
	parts.stem += [m.group("stem1")]
	parts.loop += [m.group("loop1")]
	parts.stem += [m.group("stem2")]
	parts.loop += [m.group("loop2")]
	parts.stem += [m.group("stem3")]
	leftover   =   m.group("tail")

	if ("regex" in debug):
		print >>stderr, "  stem1: \"%s\"" % parts.stem1
		print >>stderr, "  loop1: \"%s\"" % parts.loop1
		print >>stderr, "  stem2: \"%s\"" % parts.stem2
		print >>stderr, "  loop2: \"%s\"" % parts.loop2
		print >>stderr, "  stem3: \"%s\"" % parts.stem3
		print >>stderr, "  tail:  \"%s\"" % leftover

	# if there's any left over, try to re-parse it into more loops and stems

	if (leftover != ""):
		leftover = reparse_leftover(parts,leftover)

	parts.tail = None if (leftover == "") else leftover

	return parts


# parse_as_g_quad_2_stems--
#	Gives preference to parsing the string into a two-stem object.

def parse_as_g_quad_2_stems(seq):
	if ("regex" in debug):
		print >>stderr, "seq = \"%s\"" % seq

	hasLongLoop = hasBulge = False
	m = gQuad21Full.match(seq)

	if (m == None):
		m = gQuad21LongLoopFull.match(seq)
		if (m != None): hasLongLoop = True

	if (allowBulges):
		if (m == None):
			m = gQuad21BulgesFull.match(seq)
			if (m != None): hasBulge = True
		if (m == None):
			m = gQuad21BulgesLongLoopFull.match(seq)
			if (m != None): hasLongLoop = hasBulge = True

	if (m == None):     # shouldn't happen, assuming the inputs are really
		return None     # .. g-quadruplex motifs

	parts = GQuadParts()
	parts.hasLongLoop = hasLongLoop
	parts.hasBulge    = hasBulge
	parts.stem  = []
	parts.loop  = []
	parts.stem += [m.group("stem1")]
	parts.loop += [m.group("loop1")]
	parts.stem += [m.group("stem2")]
	leftover   =   m.group("tail")

	if ("regex" in debug):
		print >>stderr, "  stem1: \"%s\"" % parts.stem1
		print >>stderr, "  loop1: \"%s\"" % parts.loop1
		print >>stderr, "  stem2: \"%s\"" % parts.stem2
		print >>stderr, "  tail:  \"%s\"" % leftover

	# if there's any left over, try to re-parse it into more loops and stems

	if (leftover != ""):
		leftover = reparse_leftover(parts,leftover)

	parts.tail = None if (leftover == "") else leftover

	return parts


# parse_as_g_quad_1_stem--
#	Gives preference to parsing the string into a one-stem object.

def parse_as_g_quad_1_stem(seq):
	if ("regex" in debug):
		print >>stderr, "seq = \"%s\"" % seq

	hasLongLoop = hasBulge = False
	m = gQuad10Full.match(seq)

	if (m == None):
		m = gQuad10LongLoopFull.match(seq)
		if (m != None): hasLongLoop = True

	if (allowBulges):
		if (m == None):
			m = gQuad10BulgesFull.match(seq)
			if (m != None): hasBulge = True
		if (m == None):
			m = gQuad10BulgesLongLoopFull.match(seq)
			if (m != None): hasLongLoop = hasBulge = True

	if (m == None):     # shouldn't happen, assuming the inputs are really
		return None     # .. g-quadruplex motifs

	parts = GQuadParts()
	parts.hasLongLoop = hasLongLoop
	parts.hasBulge    = hasBulge
	parts.stem  = []
	parts.loop  = []
	parts.stem += [m.group("stem1")]
	leftover   =   m.group("tail")

	if ("regex" in debug):
		print >>stderr, "  stem1: \"%s\"" % parts.stem1
		print >>stderr, "  tail:  \"%s\"" % leftover

	# if there's any left over, try to re-parse it into more loops and stems

	if (leftover != ""):
		leftover = reparse_leftover(parts,leftover)

	parts.tail = None if (leftover == "") else leftover

	return parts


# reparse_leftover--
#	Try to re-parse a leftover tail into more loops and stems
#
# Note that this may modify parts.

def reparse_leftover(parts,leftover):

	# as long as there's any left over, try to re-parse it into more loops and
	# stems
	#
	# nota bene: the regular expressions used here are designed to find the
	#            *shortest* match

	while (leftover != ""):
		m = loopAndStemShortest.match(leftover)

		if (m == None):
			m = loopAndStemLongLoopShortest.match(leftover)
			if (m != None): hasLongLoop = True

		if (allowBulges):
			if (m == None):
				m = loopAndStemBulgesShortest.match(leftover)
				if (m != None): hasBulge = True
			if (m == None):
				m = loopAndStemBulgesLongLoopShortest.match(leftover)
				if (m != None): hasLongLoop = hasBulge = True

		if (m == None):
			break

		loop = m.group("loop")
		stem = m.group("stem")
		leftover = leftover[len(loop)+len(stem):]
		parts.loop += [loop]
		parts.stem += [stem]

	# if there's still any left over, try to re-parse it, in combination with
	# the final stem, into a longer stem; this is to resolve the issue of the
	# shortest match use above not including everything it might in the final
	# stem

	if (leftover != ""):
		stemPlusLeftover = parts.stem[-1] + leftover

		m = stemLongest.match(stemPlusLeftover)

		if (m == None):
			m = stemLongLoopLongest.match(stemPlusLeftover)
			if (m != None): hasLongLoop = True

		if (m != None):
			stem = m.group("stem")
			leftover = stemPlusLeftover[len(stem):]
			parts.stem[-1] = stem

	return leftover


# read_gquad_bed--
#	Yield the next g-quadruplex from a bed file

class GQuad: pass

def read_gquad_bed(f,fName=None):
	if (fName == None): fName = "input"

	lineNumber = 0
	for line in f:
		lineNumber += 1
		line = line.strip()
		if (line == "") or (line.startswith("#")):
			continue

		fields = line.split()
		assert (len(fields) in [4,5,6,7]), \
		      "wrong number of fields at line %s in %s (got %d expected 4, 5, 6, or 7):\n%s" \
		    % (lineNumber,fName,len(fields),line)

		try:
			chrom    =     fields[0]
			start    = int(fields[1])
			end      = int(fields[2])
			motifSeq =     fields[3]
			strand   =     fields[5] if (len(fields) >= 6) else None
		except ValueError:
			assert (False), \
			      "bad line, interval is not integers (line %s in %s):\n%s" \
			    % (lineNumber,fName,line)

		assert (start < end), \
		      "bad line, empty interval (line %s in %s):\n%s" \
		    % (lineNumber,fName,line)

		if (allowBadLength):
			if (len(motifSeq) != end-start): end = start + len(motifSeq)
		else:
			assert (len(motifSeq) == end-start), \
			      "bad line, sequence length doesn't match interval length (line %s in %s):\n%s" \
			    % (lineNumber,fName,line)

		assert (strand in [None,"+","-"]), \
		      "bad line, strand is not + nor - (line %s in %s):\n%s" \
		    % (lineNumber,fName,line)

		g4 = GQuad()
		g4.line     = " ".join(line.split())
		g4.chrom    = chrom
		g4.start    = start
		g4.end      = end
		g4.strand   = strand
		g4.motifSeq = motifSeq
		yield g4


# reverse_complement--

complementMap = maketrans("ACGTSWRYMKBDHVNacgtswrymkbdhvn",
                          "TGCASWYRKMVHDBNtgcaswyrkmvhdbn")

def reverse_complement(nukes):
	return nukes[::-1].translate(complementMap)


# int_with_unit--
#	Parse a string as an integer, allowing unit suffixes

def int_with_unit(s):
	if (s.endswith("K")):
		multiplier = 1000
		s = s[:-1]
	elif (s.endswith("M")):
		multiplier = 1000 * 1000
		s = s[:-1]
	elif (s.endswith("G")):
		multiplier = 1000 * 1000 * 1000
		s = s[:-1]
	else:
		multiplier = 1

	try:               return          int(s)   * multiplier
	except ValueError: return int(ceil(float(s) * multiplier))


# commatize--
#	Convert a numeric string into one with commas.

def commatize(s):
	if (type(s) != str): s = str(s)
	(prefix,val,suffix) = ("",s,"")
	if (val.startswith("-")): (prefix,val) = ("-",val[1:])
	if ("." in val):
		(val,suffix) = val.split(".",1)
		suffix = "." + suffix

	try:    int(val)
	except: return s

	digits = len(val)
	if (digits > 3):
		leader = digits % 3
		chunks = []
		if (leader != 0):
			chunks += [val[:leader]]
		chunks += [val[ix:ix+3] for ix in xrange(leader,digits,3)]
		val = ",".join(chunks)

	return prefix + val + suffix


if __name__ == "__main__": main()
