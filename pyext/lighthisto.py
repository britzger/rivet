#!/usr/bin/env python

import os
import logging

## Work around lacking functionality in old Python versions:
## Make "sorted" a builtin function on Python < 2.4
if not 'sorted' in dir(__builtins__):
    def sorted(iterable, cmp=None, key=None, reverse=None):
        rtn = iterable
        rtn.sort(cmp)
        return rtn


from htmlentitydefs import codepoint2name
unichr2entity = dict((unichr(code), u'&%s;' % name)
                         for code,name in codepoint2name.iteritems()
                         if code != 38) # exclude "&"
def htmlescape(text, d=unichr2entity):
    if u"&" in text:
        text = text.replace(u"&", u"&amp;")
    for key, value in d.iteritems():
        if key in text:
            text = text.replace(key, value)
    return text

# Histo and Bin classes were copied from aida2flat
class Histo:
    aidaindent = "  "
    def __init__(self):
        self._bins = []
        self.path = None
        self.name = None
        self.title = None
        self.xlabel = None
        self.ylabel = None
        self._sorted = False

    def __cmp__(self, other):
        """Sort by $path/$name string"""
        return self.fullPath() > other.fullPath()

    def __str__(self):
        out = "Histogram '%s' with %d bins\n" % (self.fullPath(), self.numBins())
        out += "Title: %s\n" % self.title
        out += "XLabel: %s\n" % self.xlabel
        out += "YLabel: %s\n" % self.ylabel
        out += "\n".join([str(b) for b in self.getBins()])
        return out

    def fullPath(self):
        return os.path.join(self.path, self.name)
    fullpath = property(fullPath)

    def header(self):
        out = "# BEGIN PLOT\n"
        out += "LogY=1\n"
        out += "Title=%s\n" % self.title
        out += "XLabel=%s\n" % self.xlabel
        out += "YLabel=%s\n" % self.ylabel
        out += "# END PLOT\n"
        return out

    def asFlat(self, gnuplot=False):
        out = "# BEGIN HISTOGRAM %s\n" % self.fullPath()
        out += "AidaPath=%s\n" % self.fullPath()
        out += "Title=%s\n" % self.title
        out += "XLabel=%s\n" % self.xlabel
        out += "YLabel=%s\n" % self.ylabel
        if self.fullPath().startswith('/REF'):
            out += "PolyMarker=*\n"
            out += "ErrorBars=1\n"
        out += "## Area: %s\n" % self.area()
        out += "## Num bins: %d\n" % self.numBins()
        if gnuplot:
            out += "## xval  \tyval    \txlow    \txhigh    \tylow     \tyhigh\n"
        else:
            out += "## xlow  \txhigh   \tyval    \tyerrminus\tyerrplus\n"
        out += "\n".join([b.asFlat(gnuplot) for b in self.getBins()])
        out += "\n# END HISTOGRAM"
        return out

    def asAIDA(self):
        ind = self.aidaindent
        r = ind + '<dataPointSet name="%s" dimension="2"\n' % (
                self.name)
        if self.title is not None:
            r += ind + '    path="%s" title="%s">\n' % (
                    self.path, htmlescape(self.title))
        else:
            r += ind + '    path="%s" title="">\n' % (
                    os.path.dirname(self.path))
        if self.xlabel is not None:
            r += ind + '  <dimenstion dim="0" title="%s"/>\n' % (
                    htmlescape(self.xlabel))
        if self.ylabel is not None:
            r += ind + '  <dimenstion dim="1" title="%s"/>\n' % (
                    htmlescape(self.ylabel))
        r += ind + "  <annotation>\n"
        if self.title is not None:
            r += ind + '    <item key="Title" value="%s" sticky="true"/>\n' % (
                    htmlescape(self.title))
        else:
            r += ind + '    <item key="Title" value="" sticky="true"/>\n'
        r += ind + '    <item key="AidaPath" value="%s" sticky="true"/>\n' % (
                self.fullPath())
        # TODO: FullPath annotation item?
        # r += ind + '    <item key="FullPath" value
        r += ind + "  </annotation>\n"
        for b in self:
            r += b.asAIDA()
        r += ind + "</dataPointSet>\n"
        return r

    def numBins(self):
        return len(self._bins)

    def getBins(self):
        if not self._sorted:
            self._bins.sort()
            self._sorted = True
        return self._bins

    def setBins(self, bins):
        self._bins = bins
        self._sorted = False
        return self

    def addBin(self, bin):
        self._bins.append(bin)
        self._sorted = False
        return self

    def getBin(self, index):
        if not self._sorted:
            self._bins.sort()
            self._sorted = True
        return self.getBins()[index]

    bins = property(getBins, setBins)

    def area(self):
        return sum([bin.area() for bin in self.bins])
    getArea = area

    def __iter__(self):
        return iter(self.getBins())

    def __len__(self):
        return len(self._bins)

    def __getitem__(self, index):
        return self.getBin(index)

    def chop(self, *xranges):
        """Return a chopped histogram.

        The kept bins are defined by (xstart, xstop) pairs. The first xstart
        and last xstop can be None meaning that all is included from the
        first or up to the last bin respectively.
        Example:
            >>> hist.chop((2.5, 5.5), (7.5, None))
        """
        if len(xranges) == 0:
            raise ValueError("At least one (xstart, xstop) range is needed!")
        # check that xranges is
        laststop = xranges[0][1]
        for xr in xranges[1:]:
            if laststop >= xr[0]:
                raise ValueError("(xstart, xstop) ranges must be in numerical order!")
            laststop = xr[1]

        new = Histo()
        new.path = self.path
        new.name = self.name
        new.title = self.title
        new.xlabel = self.xlabel
        new.ylabel = self.ylabel

        irange = 0
        curran = xranges[irange]
        for b in self:
            lowok = False
            highok = False
            br = b.getXRange()
            # print curran, "->", br
            # update the current range used if we exceed the current upper
            # limit
            while (curran[1] is not None and
                    irange < len(xranges) - 1 and
                    br[0] > curran[1]):
                irange += 1
                curran = xranges[irange]

            if ((curran[0] is None or curran[0] <= br[0] or
                        br[0] <= curran[0] <= br[1]) and

                (curran[1] is None or curran[1] >= br[1] or
                        br[0] <= curran[1] <= br[1])):
                new.addBin(b)
            else:
                logging.info("Chopping bin %s:%f" % (self.fullPath(),
                             b.getBinCenter()))
        return new

    ## decorators are only available since Python 2.4
    def fromDPS(cls, dps):
        """Build a histogram from a xml dataPointSet."""
        new = cls()
        new.name = dps.get("name")
        new.title = dps.get("title")
        new.path = dps.get("path")
        axes = dps.findall("dimension")
        if (len(axes)==2):
            for a in axes:
                if (a.get("dim")=="0"):
                    new.xlabel = a.get("title")
                elif (a.get("dim")=="1"):
                    new.ylabel = a.get("title")
        points = dps.findall("dataPoint")
        numbins = len(points)
        for binnum, point in enumerate(points):
            bin = Bin()
            for d, m in enumerate(point.findall("measurement")):
                val  = float(m.get("value"))
                down = float(m.get("errorMinus"))
                up = float(m.get("errorPlus"))
                if d == 0:
                    low  = val - down
                    high = val + up
                    bin.setXRange(low, high)
                elif d == 1:
                    bin.yval = val
                    bin.yerrplus = up
                    bin.yerrminus = down
            new.addBin(bin)
        return new
    fromDPS = classmethod(fromDPS)

    def fromFlat(cls, stringbuf):
        """Build a histogram from a string buffer containing flat-format."""
        desc = {}
        new = cls()
        for line in stringbuf:
            line = line.rstrip()
            if "=" in line:
                linearray = line.split("=", 0)
                desc[linearray[0]] = linearray[1]
            else:
                linearray = line.split()
                if len(linearray) == 4:
                    new.addBin(Bin(float(linearray[0]), float(linearray[1]),
                                   float(linearray[2]),
                                   float(linearray[3]), float(linearray[3])))
                elif len(linearray) == 5:
                    new.addBin(Bin(float(linearray[0]), float(linearray[1]),
                                   float(linearray[2]),
                                   float(linearray[3]), float(linearray[4])))
                else:
                    logging.error("Unknown line format in '%s'" % (line))
        new.path, new.name = os.path.split(desc["AidaPath"])
        if desc.has_key("Title"):
            new.title = desc["Title"]
        if desc.has_key("XLabel"):
            new.title = desc["XLabel"]
        if desc.has_key("YLabel"):
            new.title = desc["YLabel"]
        return new
    fromFlat = classmethod(fromFlat)


class Bin:
    """A simple container for a binned value with an error."""
    aidaindent = "    "
    def __init__(self, xlow=None, xhigh=None, yval=0, yerrplus=0, yerrminus=0, focus=None):
        self.xlow = xlow
        self.xhigh= xhigh
        self.yval = yval
        self.yerrplus = yerrplus
        self.yerrminus = yerrminus
        self.focus= focus

    def __str__(self):
        out = "%e to %e: %e +%e-%e" % (self.xlow, self.xhigh,
                self.yval, self.yerrplus, self.yerrminus)
        return out

    def asFlat(self, gnuplot):
        if gnuplot:
            out = "%e\t%e\t%e\t%e\t%e\t%e" % (self.getBinCenter(), self.yval,
                                              self.xlow, self.xhigh, 
                                              self.yval-self.yerrminus, self.yval+self.yerrplus)
        else:
            out = "%e\t%e\t%e\t%e\t%e" % (self.xlow, self.xhigh, self.yval, self.yerrminus, self.yerrplus)
        return out

    def asAIDA(self):
        "Return this bin as AIDA formatted string."
        ind = self.aidaindent
        return (ind + "<dataPoint>\n"
            + ind
            + '  <measurement errorPlus="%e" value="%e" errorMinus="%e"/>\n' % (
                self.getXErrPlus(), self.getBinCenter(), self.getXErrMinus())
            + ind
            + '  <measurement errorPlus="%e" value="%e" errorMinus="%e"/>\n' % (
                self.yerrplus, self.yval, self.yerrminus)
            + ind + "</dataPoint>\n")

    def __cmp__(self, other):
        """Sort by mean x value (yeah, I know...)"""
        return (self.xlow + self.xhigh) > (other.xlow + other.xhigh)

    def getXRange(self):
        return (self.xlow, self.xhigh)

    def setXRange(self, xlow, xhigh):
        self.xlow = xlow
        self.xhigh = xhigh
        return self

    def getBinCenter(self):
        """Geometric middle of the bin range."""
        return self.xlow + .5*(self.xhigh - self.xlow)

    def getXErrPlus(self):
        return self.xhigh - self.getBinCenter()

    def getXErrMinus(self):
        return self.getBinCenter() - self.xlow

    def getFocus(self):
        """Mean x-value of the bin."""
        if self.focus is not None:
            return (self.xlow + self.xhigh)/2.0
        else:
            return self.focus

    def area(self):
        return self.yval * (self.xhigh - self.xlow)
    getArea = area

    def getYErr(self):
        """Get mean of +ve and -ve y-errors."""
        return (self.yerrplus + self.yerrminus)/2.0

    def setYErr(self, yerr):
        """Set both +ve and -ve y-errors simultaneously."""
        self.yerrplus = yerr
        self.yerrminus = yerr
        return self

    yerr = property(getYErr, setYErr)
