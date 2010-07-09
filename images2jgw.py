#!/usr/bin/env python

# images2jgw.py

import os, sys, commands, threading, re, math
from optparse import OptionParser

from coordinate import Point

DEBUG = False
prev_point = None
tif = False
rotate = False
idir = None
odir = './'
progress = False

__VERSION__ = '0.0.1'

def debug(txt):
    global DEBUG
    if DEBUG:
        print "%s" % txt

class FovBox(object):
    cp = {'lon': 0, 'lat': 0}
    ul = {'lon': 0, 'lat': 0}
    ll = {'lon': 0, 'lat': 0}
    ur = {'lon': 0, 'lat': 0}
    lr = {'lon': 0, 'lat': 0}
    dw = {'value': 0, 'units': ''}
    dh = {'value': 0, 'units': ''}
    az = {'value': 0, 'units': ''}
    
    def __str__(self):
        txt  = "CP %s\n" % (self.cp,)
        txt += "UL %s\n" % (self.ul,)
        txt += "LL %s\n" % (self.ll,)
        txt += "UR %s\n" % (self.ur,)
        txt += "LR %s\n" % (self.lr,)
        txt += "DW %s\n" % (self.dw,)
        txt += "DH %s\n" % (self.dh,)
        return txt


class JheadParser:
    """Data class for storing Jhead output. This will be put into SwathViewer
    XML Curser on Target data"""
    re_Filename      = r'^\s*File name\s*:\s*(.*)\n'
    re_DateTime      = r'^\s*Date/Time\s*:\s*(.*)\n'
    re_GPSLatitude   = r'^\s*GPS Latitude\s*:\s*(.*)\n'
    re_GPSLongitude  = r'^\s*GPS Longitude\s*:\s*(.*)\n'
    re_GPSAltitude   = r'^\s*GPS Altitude\s*:\s*(.*)\n'
    re_CameraMake    = r'^\s*Camera make\s*:\s*(.*)\n'
    re_CameraModel   = r'^\s*Camera model\s*:\s*(.*)\n'
    re_ExposureTime  = r'^\s*Exposure time\s*:\s*(.*)\n'
    re_GPSTimeStamp  = r'^\s*GPSTimeStamp\s*=\s*(.*)\n'
    re_GPSDateStamp  = r'^\s*GPSDateStamp\s*="(.*)"\n'
    re_Dimensions    = r'^\s*Resolution\s*:\s*([0-9]+)\sx\s([0-9]+)\n'
    
    Filename = None
    DateTime = None
    CameraMake = None
    CameraModel = None
    ExposureTime = None
    GPSLatitude = None
    GPSLongitude = None
    GPSAltitude = None 
    GPSTimeStamp = None
    GPSDateStamp = None
    GPSVersion = '2.2.0.0'
    Width = None
    Height = None
    
    Filename_parsed = False
    DateTime_parsed = False
    GPSLatitude_parsed = False
    GPSLongitude_parsed = False
    GPSAltitude_parsed = False
    CameraMake_parsed = False
    CameraModel_parsed = False
    ExposureTime_parsed = False
    GPSTimeStamp_parsed = False
    GPSDateStamp_parsed = False
    Dimensions_parsed = False
    
    VertexOffset = 0.001717906999999
    
    rotate = False
    
    #
    # Public methods
    #
    
    def __init__(self, rotate=False):
        if rotate == True:
            self.rotate = True
    
    def setControlPoints(p0, p1, p2):
        if instanceof(p0, Point) and instanceof(p1, Point) and instanceof(p2, Point):
            self.p0 = p0
            self.p1 = p1
            self.p2 = p2
            return True
        else
            return False
    
    def getRotationDegree(t=.5):
        
    
    def parse(self, line):
        if not self.Filename_parsed and self.__parse_Filename(line):
            self.Filename_parsed = True
            return True
        if not self.CameraMake_parsed and self.__parse_CameraMake(line):
            self.CameraMake_parsed = True
            return True
        if not self.CameraModel_parsed and self.__parse_CameraModel(line):
            self.CameraModel_parsed = True
            return True
        if not self.ExposureTime_parsed and self.__parse_ExposureTime(line):
            self.ExposureTime_parsed = True
            return True
        if not self.GPSLatitude_parsed and self.__parse_GPSLatitude(line):
            self.GPSLatitude_parsed = True
            return True
        if not self.GPSLongitude_parsed and self.__parse_GPSLongitude(line):
            self.GPSLongitude_parsed = True
            return True
        if not self.GPSAltitude_parsed and self.__parse_GPSAltitude(line):
            self.GPSAltitude_parsed = True
            return True
        if not self.GPSTimeStamp_parsed and self.__parse_GPSTimeStamp(line):
            self.GPSTimeStamp_parsed = True
            return True
        if not self.GPSDateStamp_parsed and self.__parse_GPSDateStamp(line):
            self.GPSDateStamp_parsed = True
            return True
        if not self.Dimensions_parsed and self.__parse_Dimensions(line):
            self.Dimensions_parsed = True
            return True
        return False
    
    def write_xml(self):
        self.__combine_DateTimeStamp()
        if not self.GPSLatitude or not self.GPSLongitude:
            return False
        
        debug('Writing FOV Box...')
        try:
            xmlfile = open('fovbox_info.txt', 'a+')
        except:
            print 'Could not open text file for writing. Exiting.'
            raise SystemExit
        
        if float(self.GPSAltitude) < 0:
            return False
        
        txt = "%s\n" % (self.Filename)
        txt += commands.getoutput('fovbox.py --lon %.8f --lat %.8f --alt %s --fovx 25.55 --fovy 17.13' %\
                                  (self.GPSLongitude,
                                   self.GPSLatitude,
                                   self.GPSAltitude))
        txt += "\n\n"
        
        xmlfile.write(txt)
        xmlfile.close()
        return True
    
    def getFov(self):
        
        global prev_point
        
        if self.rotate:
            if self.GPSLongitude and self.GPSLatitude:
                if prev_point:
                    new_point = Point(self.GPSLongitude, self.GPSLatitude)
                    azi = prev_point.geoBearingTo(new_point)
                    
                    prev_point = new_point
                else:
                    prev_point = Point(self.GPSLongitude, self.GPSLatitude)
                    azi = 0
        else:
            prev_point = None
            azi = 0
        
        self.__combine_DateTimeStamp()
        if not self.GPSLatitude or not self.GPSLongitude:
            return False
        
        if float(self.GPSAltitude) < 0:
            return None
        
        try:
            txt = commands.getoutput('fovbox.py --lon %.8f --lat %.8f --alt %s --fovx 25.55 --fovy 17.13 --azi %.2f' %\
                                      (self.GPSLongitude,
                                       self.GPSLatitude,
                                       self.GPSAltitude,
                                       azi))
            txt = txt.split('\n')
            
            fovBox = FovBox()
            
            # CP
            cp = txt[0].split()
            fovBox.cp['lon'] = cp[1]
            fovBox.cp['lat'] = cp[2]
            
            # UL
            ul = txt[1].split()
            fovBox.ul['lon'] = ul[1]
            fovBox.ul['lat'] = ul[2]
            
            # LL
            ll = txt[2].split()
            fovBox.ll['lon'] = ll[1]
            fovBox.ll['lat'] = ll[2]
            
            # UR
            ur = txt[3].split()
            fovBox.ur['lon'] = ur[1]
            fovBox.ur['lat'] = ur[2]
            
            # LR
            lr = txt[4].split()
            fovBox.lr['lon'] = lr[1]
            fovBox.lr['lat'] = lr[2]
            
            # DW
            dw = txt[5].split()
            fovBox.dw['value'] = dw[1]
            fovBox.dw['units'] = dw[2]
            
            # DH
            dh = txt[6].split()
            fovBox.dh['value'] = dh[1]
            fovBox.dh['units'] = dh[2]
            
            # AZ
            fovBox.az['value'] = azi
            fovBox.az['units'] = 'degrees'
            
            return fovBox
        except:
            return None
    
    
    #
    # Private methods
    #
    def __combine_DateTimeStamp(self):
        self.DateTime = str(self.GPSDateStamp) + 'T' + str(self.GPSTimeStamp) + '.0Z'
        return True    
    
    def __parse_Filename(self, line):
        #print "Parsing filename...", line
        match = re.match(self.re_Filename, line)
        #print self.re_Filename, line, match
        if match:
            debug('File name: ' + match.group(1))
            self.Filename = os.path.basename(match.group(1))
            return True
        return False
    
    def __parse_DateTime(self, line):
        match = re.match(self.re_DateTime, line)
        if match:
            # match.group(1) looks like "2009:04:21 09:09:59"
            # we need it to look like "2009-04-21T09:09:59.76Z"
            txt = match.group(1)
            arr = txt.split()
            arr[0] = '-'.join(arr[0].split(':'))
            txt = 'T'.join(arr)+'.0Z'
            debug('Date/Time: ' + txt)
            self.DateTime = txt
            return True
        return False
    
    def __parse_GPSTimeStamp(self, line):
        match = re.match(self.re_GPSTimeStamp, line)
        if match and not re.match(r'\?', match.group(1)):
            arr = match.group(1).split(', ')[0:3]
            hour = int((arr[0])[0:-2]) #optimization, assumes arr[0] is in form "hh/1"
            minute = int((arr[1])[0:-2]) #optimization, assumes arr[1] is in form "mm/1"
            second = arr[2].split('/')
            second = int(second[0]) / int(second[1])
            self.GPSTimeStamp = str(hour) + ':' + str(minute) + ':' + str(second)
            debug("GPSTimeStamp: %s" % self.GPSTimeStamp)
            return True
        return False
    
    def __parse_GPSDateStamp(self, line):
        match = re.match(self.re_GPSDateStamp, line)
        if match and not re.match(r'\?', match.group(1)):
           self.GPSDateStamp = match.group(1)
           debug("GPSDateStamp: %s" % self.GPSDateStamp)
           return True
        return False
    
    def __parse_GPSLatitude(self, line):
        match = re.match(self.re_GPSLatitude, line)
        if match and not re.match(r'\?', match.group(1)):
            dec = self.__parse_DMS_to_decimal(match.group(1))
            debug('GPSLatitude: %(#).15f' % {'#': dec})
            self.GPSLatitude = dec
            return True
        return False
    
    def __parse_GPSLongitude(self, line):
        match = re.match(self.re_GPSLongitude, line)
        if match and not re.match(r'\?', match.group(1)):
            dec = self.__parse_DMS_to_decimal(match.group(1))
            debug('GPSLongitude: %(#).15f' % {'#': dec})
            self.GPSLongitude = dec
            return True
        return False
    
    def __parse_GPSAltitude(self, line):
        """line is in the form of `208.0m`. We need this into float format."""
        match = re.match(self.re_GPSAltitude, line)
        if match and not re.match(r'\?', match.group(1)):
            txt = match.group(1)
            height = float(txt[0:-1])
            debug('GPSAltitude: %s m' % height)
            self.GPSAltitude = height
            return True
        return False
    
    def __parse_CameraMake(self, line):
        match = re.match(self.re_CameraMake, line)
        if match and not re.match(r'\?', match.group(1)):
            debug('Camera make: %s' % match.group(1))
            self.CameraMake = match.group(1)
            return True
        return False
    
    def __parse_CameraModel(self, line):
        match = re.match(self.re_CameraModel, line)
        if match and not re.match(r'\?', match.group(1)):
            debug('Camera model: %s' % match.group(1))
            self.CameraModel = match.group(1)
            return True
        return False
    
    def __parse_ExposureTime(self, line):
        """line is in the form of `0.0001 s  (1/8000)`, we need just the
        exposure time in seconds"""
        match = re.match(self.re_ExposureTime, line)
        if match and not re.match(r'\?', match.group(1)):
            txt = match.group(1)
            arr = txt.split()
            sec = float(arr[0])
            debug('Exposure time: %s s' % sec)
            self.ExposureTime = sec
            return True
        return False
    
    def __parse_DMS_to_decimal(self, line):
        """line is in the form of `N 64d 51.5500m  0s`, we need this in decimal
        notation, i.e., `64.8473`"""
        txt = line.split()
        sign = 0
        
        pos = txt[0]
        if pos == 'N' or pos == 'E':
            sign = 1
        if pos == 'S' or pos == 'W':
            sign = -1
        if sign == 0:
            return None
        
        # convert '147d' to the integer 147
        deg = int((txt[1])[0:-1])
        # convert '51.5500m' to the float 51.5500
        min = float((txt[2])[0:-1])
        # etc
        sec = float((txt[3])[0:-1])
        # add minutes and seconds
        minsec = min+sec
        # divide by 60*60 to get the ratio
        minsec_dec = minsec / 60.0
        # add the result to the degrees
        dec_dms = deg + minsec_dec
        # return the decimal number with the proper sign as noted above
        # if sign is 0, the result is 0, and thus the string is ill formed.
        return dec_dms*sign
    
    def __parse_Dimensions(self, line):
        match = re.match(self.re_Dimensions, line)
        if match:
            self.Width = match.group(1)
            self.Height = match.group(2)
            debug('Width: %s, Height: %s' % (self.Width, self.Height))
            return True
        else:
            return False


def parse_args():
    global __VERSION__, tif, rotate, idir, odir, progress, DEBUG
    
    # option parser (see http://docs.python.org/library/optparse.html)
    usage = "%prog [options] FILE.jpg [FILE2.jpg ...]"
    parser = OptionParser(usage=usage, version="%prog "+__VERSION__)
    parser.add_option("-t", "--tif", action="store_true", dest="tif", help="generate GeoTiffs (must have GDAL in your path)")
    parser.add_option("-r", "--rotate", action="store_true", dest="rotate", help="account for rotation from the sequence of images")
    parser.add_option("-i", "--idir", dest="idir", help="input directory to scan for jpeg images", metavar="IDIR")
    parser.add_option("-o", "--odir", dest="odir", help="output directory to place .jgw files (default: './')", default="./", metavar="ODIR")
    parser.add_option("-p", "--progress-indicator", dest="p", help="show the progress indicator", action="store_true")
    parser.add_option("-v", "--verbose", dest="v", help="show lots of verbose output; disables option -p", action="store_true")
    
    (options, args) = parser.parse_args()
    
    # Parse arguments
    if options.tif == True:
        tif = True
    
    if options.rotate == True:
        rotate = True
    
    if options.idir != None:
        if not os.path.isdir(options.idir):
            parser.error("Input directory `%s' does not exist." % (options.idir,))
        idir = options.idir
        if not idir.endswith('/'):
            idir += '/'
    
    if options.odir != None:
        if not os.path.isdir(options.odir):
            parser.error("Output directory `%s' does not exist." % (options.odir,))
        odir = options.odir
        if not odir.endswith('/'):
            odir += '/'
    
    if options.p == True:
        progress = True
    
    if options.v == True:
        DEBUG = True
        progress = False
    
    images = []
    
    if idir:
        dirList = os.listdir(idir)
        for fname in dirList:
            if os.path.isfile(idir+fname) and (
                fname.lower().endswith('.jpg') or
                fname.lower().endswith('.jpeg')
            ):
                images.append(idir+fname)
    
    for i in args:
        if os.path.isfile(i) and (
            i.lower().endswith('.jpg') or
            i.lower().endswith('.jpeg')
        ):
            images.append(i)
    if not images:
        parser.error("Please include some images to parse: FILE.jpg [FILE2.jpg ...] or give me an input directory to grab images from.")
    return images


def prog(txt):
    global progress
    if progress:
        print txt,
        sys.stdout.flush()


def main():
    global tif, rotate, odir
    
    images = parse_args()
    errors_output = ""
    images_line = (" ").join(images)
    prog("Parsing %d file(s)...\n" % (len(images),))
    prog("Gathering jhead output...")
    jhead = commands.getoutput('jhead -v -autorot %s' % (images_line,))
    prog("Done.\n")
    #print jhead
    exif = jhead.split("\n\n")
    
    count = 0
    percent = 0
    maxnum = len(exif)
    prog(" 0%...")
    
    for i in exif:
        count += 1
        perc = int((float(count)/float(maxnum)) * 10)
        if perc > percent:
            while percent < perc:
                percent += 1;
                txt = "%d" % (percent*10,)
                txt += "%..."
                prog(txt)
        jheadParser = JheadParser(rotate=rotate)
        #jhead = commands.getoutput('jhead -v %s' % (i,))
        #print jhead
        ex = i.split('\n')
        #print jhead
        for line in ex:
            jheadParser.parse(line+'\n')
        fov = jheadParser.getFov()
        debug(fov)
        if isinstance(fov, FovBox):
            #print jheadParser.Filename
            #print jheadParser.Width
            #print jheadParser.Height
            dw = math.sqrt(
                    math.pow(float(fov.ur['lon']) - float(fov.ul['lon']), 2) +
                    math.pow(float(fov.ur['lat']) - float(fov.ul['lat']), 2)
                 ) / float(jheadParser.Width)
            dh = math.sqrt(
                    math.pow(float(fov.ul['lat']) - float(fov.ll['lat']), 2) +
                    math.pow(float(fov.ul['lon']) - float(fov.ll['lon']), 2)
                 ) / float(jheadParser.Height)
            
            if jheadParser.Filename.lower().endswith('.jpg'):
                fn = odir+jheadParser.Filename[:-4] + '.jgw'
                tf = odir+jheadParser.Filename[:-4] + '.tif'
                tfw = odir+jheadParser.Filename[:-4] + '.tfw'
            elif jheadParser.Filename.lower().endswith('.jpeg'):
                fn = odir+jheadParser.Filename[:-5] + '.jgw'
                tf = odir+jheadParser.Filename[:-5] + '.tif'
                tfw = odir+jheadParser.Filename[:-5] + '.tfw'
            else:
                fn = odir+jheadParser.Filename + '.jgw'
                tf = odir+jheadParser.Filename + '.tif'
                tfw = odir+jheadParser.Filename + '.tfw'
            
            if tif:
                GCP1 = "%d %d %s %s" % (int(jheadParser.Width)/2, int(jheadParser.Height)/2, fov.cp['lon'], fov.cp['lat'])
                GCP2 = "0 0 %s %s" % (fov.ul['lon'], fov.ul['lat'])
                GCP3 = "0 %d %s %s" % (int(jheadParser.Height), fov.ll['lon'], fov.ll['lat'])
                GCP4 = "%d 0 %s %s" % (int(jheadParser.Width), fov.ur['lon'], fov.ur['lat'])
                GCP5 = "%d %d %s %s" % (int(jheadParser.Width), int(jheadParser.Height), fov.lr['lon'], fov.lr['lat'])
                
                #                                     images[count-1]     GCP1    GCP2    GCP3    GCP4    GCP5                                             odir
                debug(commands.getoutput("gdal_translate -of GTiff %s -gcp %s -gcp %s -gcp %s -gcp %s -gcp %s -a_srs \"+proj=latlong +datum=WGS84 +no_defs\" %stemp.tif" % (images[count-1], GCP1, GCP2, GCP3, GCP4, GCP5, odir)))
                debug(commands.getoutput("gdalwarp -t_srs \"+proj=laea +lat_0=60.00 +lon_0=-180.000 +x_0=0 +y_0=0 +datum=WGS84 +units=m  +no_defs\" -dstnodata 0 -dstalpha -co \"TILED=YES\" -co \"ALPHA=YES\" -multi %stemp.tif %s" % (odir, tf)))
                debug(commands.getoutput("rm %stemp.tif" % (odir,)))
                debug(commands.getoutput("listgeo -tfw %s" % (tf,)))
            
            try:
                f = open(fn, 'w')
                
                # See http://www.omg.unb.ca/~jonnyb/processing/geotiff_tifw_format.html
                txt = "%.12f\r\n%.12f\r\n%.12f\r\n%.12f\r\n%s\r\n%s\r\n" % (
                         dw *  math.cos(math.radians(fov.az['value'])),
                         dw * -math.sin(math.radians(fov.az['value'])),
                         dh * -math.sin(math.radians(fov.az['value'])),
                        -dh *  math.cos(math.radians(fov.az['value'])),
                        fov.ul['lon'],
                        fov.ul['lat']
                )
                
                debug("#%s" % (fn))
                debug(txt)
                debug(" ")
                f.write(txt)
                f.close()
            except:
                pass
        else:
            errors_output += "Could not parse file `%s'... (check attributes such as correct altitude, longitude, or latitude)\r\n" % (os.path.basename(images[count-1]),)
    
    if errors_output:
        try:
            f = open(odir+"images2jgw_errors.txt", "w")
            f.write("%d file(s) had parsing errors:\r\n" % (len(errors_output.split("\n"))-1,))
            f.write(errors_output)
            f.close()
            print "\nThere were some errors during processing. See `%simages2jgw_errors.txt' for details." % (odir,)
        except:
            print "\nThere were some errors during processing. I couldn't write the errors to a file, so the output is below:"
            sys.stderr.write("%d file(s) had parsing errors:\r\n" % (len(errors_output.split("\n"))-1,))
            sys.stderr.write(errors_output)

if __name__ == '__main__':
    main()
