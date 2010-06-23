# File: coordinate.py
# Coordinate module used to transform coordinates and perform calculations.
#
# About:
# Author: Jonathan Sawyer
#
# Copyright (C) 2009, 2010, Geographic Information Network of Alaska, University of Alaska Fairbanks
# All rights reserved. You may not use any portion of this code without prior written
# permission fromn the Geographic Information Network of Alaska or the University of Alaska
# Fairbanks.

#
# Module: coordinate
#
# Description:
# Coordinate library that contains these classes:
# {<Coordinate>} - Parses and represents geographic coordinates (latitude and
#                  longitude).
# {<Point>} - Contains methods to calculate range, bearing, and create points
#             from those.
#
# Decimal degrees:
#    123.12341234[NSEW]
#    [NSEW+-]123.12341234
#
# Degrees-minutes-seconds:
#     [NSEW+-]11d 22m 33.33333s
#     11d 22m 33.33333s[NSEW]
#     11:22:33.33333[NSEW]
#     [NSEW+-]11:22:33.33333
#
# Examples:
#     -64.12347874
#     44D 10M 32.123S S
#     N 89.1234
#

import re
import math

#
# Constants: Radius
# EARTH_RADIUS_MI - in miles (According to IUGG)
# EARTH_RADIUS_FT - derived from <EARTH_RADIUS_MI>
# EARTH_RADIUS_KM - in kilometers (According to IUGG)
# EARTH_RADIUS_M  - derived from <EARTH_RADIUS_KM>
# EARTH_RADIUS_NMI - in nautical miles (According to IUGG)
# EARTH_RADIUS - default radius
# DEG2RAD - factor to convert degrees to radians (PI/180)
# RAD2DEG - factor to convert radians to degrees (180/PI)
#
EARTH_RADIUS_MI = 3958.761
EARTH_RADIUS_FT = EARTH_RADIUS_MI * 5280.0
EARTH_RADIUS_KM = 6371.009
EARTH_RADIUS_M = EARTH_RADIUS_KM * 1000.0
EARTH_RADIUS_NMI = 3440.069
EARTH_RADIUS = EARTH_RADIUS_KM
DEG2RAD =  0.01745329252
RAD2DEG = 57.29577951308

#
# Class: Coordinate
#
# Parse an input into a known coordinate.
#
class Coordinate:
    coord = None     # string
    
    coord_dd = None  # real number
    coord_dms = None # 3-tuple
    direction = None # string in [N,S,E,W]
    
    groups = None    # tuple
    
    # NOTE: Order Matters. The first two produce a 4-tuple, the other four
    # produce a 6-tuple
    regs = [
        # 123.12341234[NSEW]
        re.compile(r'''^\s*()([0-9]{0,3}(\.[0-9]{1,15})?)\s*([nNsSeEwW]?)\s*$'''),
        
        # [NSEW+-]123.12341234
        re.compile(r'''^\s*([nNsSeEwW+-]?)\s*([0-9]{0,3}(\.[0-9]{1,15})?)()\s*$'''),
        
        # [NSEW+-]11d22m33.333333s
        re.compile(r'''\s*([nNsSeEwW+-]?)\s*([0-9]{0,3})\s*[dD]\s*([0-9]{0,2})\s*[mM]\s*([0-9]{0,2}(\.[0-9]{1,15})?)\s*[sS]()\s*$'''),
        
        # 11d22m33.333333s[NSEW]
        re.compile(r'''\s*()([0-9]{0,3})\s*[dD]\s*([0-9]{0,2})\s*[mM]\s*([0-9]{0,2}(\.[0-9]{1,15})?)\s*[sS]\s*([nNsSeEwW]?)\s*$'''),
        
        # [NSEW+-]11:22:33.333333
        re.compile(r'''\s*([nNsSeEwW+-]?)\s*([0-9]{0,3})\s*:\s*([0-9]{0,2})\s*:\s*([0-9]{0,2}(\.[0-9]{1,15})?)()\s*$'''),
        
        # 11:22:33.333333[NSEW]
        re.compile(r'''\s*()([0-9]{0,3})\s*:\s*([0-9]{0,2})\s*:\s*([0-9]{0,2}(\.[0-9]{1,15})?)\s*([nNsSeEwW]?)\s*$'''),
        
        ]
    
    
    # 
    # Method: Constructor
    # 
    # Parameters:
    # coord - {string} The unformatted coordinate to parse
    # axis  - {string} The axis in which the unformatted coordinate is in.
    #                  Either Longitude or Latitude.
    #
    def __init__(self, coord, axis=''):
        try:
            c = float(coord)
            self.coord = "%.12f" % c
        except:
            self.coord = str(coord)
        self._parse(axis)
    
    
    #
    # Method: _parse
    #
    # (Private) Parse the objects own data.
    #
    # Parameters:
    # axis - {string} The axis in which the unformatted coordinate is in.
    #                 Either Longitude or Latitude
    #
    def _parse(self, axis=''):
        # Match the unformatted coordinate by the given regular expressions
        # list.
        position = 1
        for r in self.regs:
            m = r.match(self.coord)
            if m:
                self.groups = m.groups()
                break
            else:
                position += 1
                continue
        
        # If there was no match, return.
        if not self.groups:
            return
        
        # Parse direction and coordinate.
        if position <= 2:
            if self.groups[0]:
                self._parse_direction(self.groups[0], axis)
            elif self.groups[3]:
                self._parse_direction(self.groups[3], axis)
            else:
                self._parse_direction("+", axis)
            self._parse_dd(self.groups[1])
        elif position <= 6:
            if self.groups[0]:
                self._parse_direction(self.groups[0], axis)
            elif self.groups[5]:
                self._parse_direction(self.groups[5], axis)
            else:
                self._parse_direction("+", axis)
            self._parse_dms((self.groups[1], self.groups[2], self.groups[3]))
    
    
    #
    # Method: _parse_direction
    #
    # (Private) Parse the direction given by `direction` and `axis`.
    #
    # Parameters:
    # direction - {string} The direction, given by a string of either +, -, N,
    #                      S, E, or W.
    # axis - {string} The axis in which the unformatted coordinate is in.
    #                 Either Longitude or Latitude.
    #
    def _parse_direction(self, direction, axis):
        direction = direction.lower()
        axis = axis.lower()
        
        if direction in ['n', 's', 'e', 'w']:
            self.direction = direction.upper()
            return
        
        if axis in ['latitude', 'lat']:
            if direction == '+' or direction == '':
                self.direction = 'N'
            else:
                self.direction = 'S'
        elif axis in ['longitude', 'lon', 'long']:
            if direction == '+' or direction == '':
                self.direction = 'E'
            else:
                self.direction = 'W'
        else:
            if direction == '+' or direction == '':
                self.direction = 'N'
            else:
                self.direction = 'S'
    
    
    #
    # Method: _parse_dms
    #
    # (Private) Parse the degrees-minutes-seconds of a 3-tuple. Returns nothing.
    #
    # Parameters:
    # dms - {3-tuple} A 3-tuple representing the degrees, minutes, and seconds
    #                 of a coordinate.
    #
    def _parse_dms(self, dms):
        if len(dms) != 3 or not isinstance(dms, tuple):
            return
        
        self.coord_dms = (
            int(dms[0]),
            int(dms[1]),
            float(dms[2])
            )
        
        self.coord_dd = self._dms2dd(self.coord_dms)
    
    
    #
    # Method: _parse_dd
    #
    # (Private) Parse the decimal degrees of a float. Returns nothing.
    #
    # Parameters:
    # dd - {float} The decimal degree representation of the coordinate.
    #
    def _parse_dd(self, dd):
        self.coord_dd = float(dd)
        
        self.coord_dms = self._dd2dms(self.coord_dd)
    
    
    #
    # Method: _dd2dms
    #
    # (Private) Convert decimal degrees into degrees-minutes-seconds.
    #
    # Parameters:
    # coord - {float} The decimal degree representation of the coordinate.
    #
    # Returns:
    # {3-tuple} A 3-tuple representing the degrees, minutes, and seconds of a
    # coordinate.
    #
    def _dd2dms(self, coord):
        import math
        
        if not isinstance(coord, float):
            return
        
        d = math.floor(coord)
        ms = coord - d
        m = 60 * ms
        s = 60 * (m - math.floor(m))
        
        return (int(d), int(m), round(float(s), 6))
    
    
    #
    # Method: _dms2dd
    #
    # (Private) Convert degrees-minutes-seconds into decimal degrees.
    #
    # Parameters:
    # coord - {3-tuple} A 3-tuple representing the degrees, minutes, and
    #                   seconds of a coordinate.
    #
    # Returns:
    # {float} The decimal degree representation of the coordinate.
    #
    def _dms2dd(self, coord):
        if len(coord) != 3 or not isinstance(coord, tuple):
            return 0.0
        
        if not isinstance(coord[0], int) or \
           not isinstance(coord[1], int) or \
           not isinstance(coord[2], float):
            return 0.0
        
        dd = coord[0] + (coord[1]/60.0) + (coord[2]/(60.0*60.0))
        
        return round(float(dd), 6)
    
    
    #
    # Method: getDms
    #
    # Returns:
    # {string} The degrees-minutes-seconds representation of the coordinate.
    #
    def getDms(self):
        if self.coord_dms and self.direction:
            return "%dd %dm %.6fs %s" % (
                self.coord_dms[0],
                self.coord_dms[1],
                self.coord_dms[2],
                self.direction
                )
        else:
            return None
    
    
    #
    # Method: getDdString
    #
    # Returns:
    # {string} The decimal degree representation of the coordinate with
    #          direction. Does not use the negative sign but instead either
    #          'W' or 'S'.
    #
    def getDdString(self):
        if self.coord_dd or (self.coord_dd != None and int(self.coord_dd) == 0) and self.direction:
            return "%.12f %s" % (
                self.coord_dd,
                self.direction
                )
        else:
            return None
    
    
    #
    # Method: getDdNum
    #
    # Returns:
    # {string} The decimal degree representation of the coordinate.
    #
    def getDdNum(self):
        if self.direction in ['W', 'S']:
            return "-%.12f" % self.coord_dd
        else:
            if self.coord_dd or (self.coord_dd != None and int(self.coord_dd) == 0):
                return "%.12f" % self.coord_dd
            else:
                return None
    
    
    #
    # Method: getDd
    #
    # Returns:
    # {string} The decimal degree representation of the coordinate.
    #
    def getDd(self):
        return self.getDdNum()
    
    
    #
    # Method: getCoord
    #
    # Returns:
    # {string} The coordinate as entered in by the user.
    #
    def getCoord(self):
        return "%s" % self.coord
    
    
    #
    # Method: isValid
    #
    # Returns:
    # {boolean} True if the coordinate parsing succeeded, False otherwise.
    #
    def isValid(self):
        try:
            return (self.coord_dd or int(self.coord_dd) == 0) and self.coord_dms and self.direction
        except:
            return False
    
    
    #
    # Method: debug
    #
    # Description:
    # Prints debug information.
    #
    def debug(self):
        print "Coord:    %.12f" % self.getCoord()
        print "DD:       %.12f" % self.getDd()
        print "DDNum:    %.12f" % self.getDdNum()
        print "DDString: %.12f" % self.getDdString()
        print "DMS:      %.12f" % self.getDms()


#
# Class: Point
#
# Requires:
# <Coordinate> - Coordinate parsing class used in the constructor of this
#                class.
#
# Reference:
# Williams, Ed, 2000, "Aviation Formulary V1.43" web page
# http://williams.best.vwh.net/avform.htm
#
class Point:
    x = None
    y = None
    
    
    #
    # Method: Constructor
    #
    # Parameters:
    # x - {string} The x coordinate of the point.
    # y - {string} The y coordinate of the point.
    #
    def __init__(self, x, y):
        self.x = float(Coordinate(x, "Lon").getDd())
        self.y = float(Coordinate(y, "Lat").getDd())
    
    
    #
    # Method: __str__
    #
    # Description:
    # Overloaded print function
    #
    # Returns:
    # {string} The string representation of the point.
    #
    def __str__(self):
        return "(%.12f, %12f)" % (self.x, self.y)
    
    
    #
    # Method: geoDistanceTo
    #
    # Parameters:
    # point - {<Point>}
    #
    # Returns:
    # {float} Great Circle distance to Point. Coordinates must be in decimal
    # degrees.
    #
    def geoDistanceTo(self, point, units='km'):
        global EARTH_RADIUS_KM, EARTH_RADIUS_MI, EARTH_RADIUS_NMI, DEG2RAD
        
        x = [0, 0]
        y = [0, 0]
        radius = 0
        
        # Calculates the radius of the earth used
        if units and units.lower() == 'km':
            radius = EARTH_RADIUS_KM
        elif units and units.lower() == 'm':
            radius = EARTH_RADIUS_M
        elif units and units.lower() == 'mi':
            radius = EARTH_RADIUS_MI
        elif units and units.lower() == 'ft':
            radius = EARTH_RADIUS_FT
        elif units and units.lower() == 'nmi':
            radius = EARTH_RADIUS_NMI
        else:
            radius = EARTH_RADIUS_KM
        
        x[0] = self.x * DEG2RAD
        x[1] = point.x * DEG2RAD
        y[0] = self.y * DEG2RAD
        y[1] = point.y * DEG2RAD
        
        a = math.pow( math.sin(( y[1]-y[0] ) / 2.0 ), 2)
        b = math.pow( math.sin(( x[1]-x[0] ) / 2.0 ), 2)
        c = math.pow(( a + math.cos( y[1] ) * math.cos( y[0] ) * b ), 0.5)
        
        return 2 * math.asin( c ) * radius
    
    
    #
    # Method: geoBearingTo
    #
    # Parameters:
    # point - {<Point>}
    #
    # Returns:
    # {float} - The bearing clockwise from North in degrees.
    #
    def geoBearingTo(self, point):
        global EARTH_RADIUS_KM, EARTH_RADIUS_MI, EARTH_RADIUS_NMI
        global DEG2RAD, RAD2DEG
        
        x = [0, 0]
        y = [0, 0]
        bearing = None
        adjust = None
        
        x[0] = self.x * DEG2RAD
        x[1] = point.x * DEG2RAD
        y[0] = self.y * DEG2RAD
        y[1] = point.y * DEG2RAD
        
        a = math.cos(y[1]) * math.sin(x[1] - x[0])
        b = math.cos(y[0]) * math.sin(y[1]) - math.sin(y[0]) * math.cos(y[1]) * math.cos(x[1] - x[0])
        
        if a == 0 and b == 0:
            return 0.0
        
        if b == 0:
            if a < 0:
                return 270.0
            else:
                return 90.0
        
        if b < 0:
            adjust = math.pi
        else:
            if a < 0:
                adjust = 2 * math.pi
            else:
                adjust = 0
        
        return (math.atan(a/b) + adjust) * RAD2DEG
    
    
    #
    # Method: geoWaypoint
    #
    # Parameters:
    # distance - {float}
    # bearing  - {float}
    # units    - {string}
    #
    # Returns:
    # {<Point>} - The generated point given by distance and bearing.
    #
    def geoWaypoint(self, distance, bearing, units='km'):
        global EARTH_RADIUS_KM, EARTH_RADIUS_MI, EARTH_RADIUS_NMI
        global DEG2RAD, RAD2DEG
        
        wp = Point(0, 0)
        radius = 0
        
        # Calculates the radius of the earth used
        if units and units.lower() == 'km':
            radius = EARTH_RADIUS_KM
        elif units and units.lower() == 'm':
            radius = EARTH_RADIUS_M
        elif units and units.lower() == 'mi':
            radius = EARTH_RADIUS_MI
        elif units and units.lower() == 'ft':
            radius = EARTH_RADIUS_FT
        elif units and units.lower() == 'nmi':
            radius = EARTH_RADIUS_NMI
        else:
            radius = EARTH_RADIUS_KM
        
        x = self.x * DEG2RAD
        y = self.y * DEG2RAD
        radBearing = bearing * DEG2RAD
        
        # Convert arc distance to radians
        c = distance / radius
        
        wp.y = math.asin( math.sin(y) * math.cos(c) + math.cos(y) * math.sin(c) * math.cos(radBearing)) * RAD2DEG
        
        a = math.sin(c) * math.sin(radBearing)
        b = math.cos(y) * math.cos(c) - math.sin(y) * math.sin(c) * math.cos(radBearing)
        
        if b == 0:
            wp.x = self.x
        else:
            wp.x = self.x + math.atan(a/b) * RAD2DEG
        
        return wp
    
    #
    # Method: rotate
    #
    # Parameters:
    # point   - {<Point>} The point in which to rotate about.
    # degrees - {float} The rotation angle in degrees.
    #
    # Returns:
    # {<Point>} - The generated point given by point and degrees.
    #
    def rotate(self, point, degrees):
        if degrees == 0.0 or degrees == 360.0:
            return self
        
        bearing = point.geoBearingTo(self)
        distance = point.geoDistanceTo(self)
        
        return point.geoWaypoint(distance, bearing+degrees)
