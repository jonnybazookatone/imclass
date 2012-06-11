# Author: J Elliott
# Date: 07.11.2011
#-------------------

# imports
from pyraf import iraf
import matplotlib.pyplot as plt
import sys, os, signal
import subprocess
import pyfits
import scipy
import numpy
import glob

import python.lib.sextractor as sextractor

HOTPANTSPATH="/home/sw/bin"

## Usagephot

Usage = """FITS image class, to be used with pypants.py and pyremap.py
for example:
	from image.py import imFits
"""

# Image Class
# To be used with pypants.py and pyremap.py

class imObject(object):
	
	def __init__(self):

		self._id = ""
		self._pixelx = ""
		self._pixely = ""
		self._mag = ""
		self._ra = ""
		self._dec = ""
		self._parentimage = ""
		self._appMag = ""
		self._appMagErr = ""
		self._Name = ""
		self._Box = ""
                self._NPIX = ""
                self._MEAN = ""
                self._STDDEV = ""
                self._MIN = ""
                self._MAX = ""
		self._MEDIAN = ""
		self._midMJD = ""
		self._skyFlux = ""
		self._flux = ""
	
		# Apertures for photometry
                self._fap = 30.0
                self._fdan = 50.0
                self._fan = 55.0


	def setAps(self, fap, fdan, fan, scale=1.0):
		self._fap = fap * scale
		self._fdan = fdan * scale
		self._fan = fan * scale
	def getAps(self):
		return self._fap, self._fdan, self._fan

	def setBox(self, Box):
		self._Box = Box
	def setMEAN(self, MEAN):
		self._MEAN = MEAN
	def setNPIX(self, NPIX):
		self._NPIX = NPIX
	def setSTDDEV(self, STDDEV):
		self._STDDEV = STDDEV
	def setMIN(self, MIN):
		self._MIN = MIN
	def setMAX(self, MAX):
		self._MAX = MAX
	def setMEDIAN(self, MEDIAN):
		self._MEDIAN = MEDIAN            
    
	def copyObjectProperties(self, newObject):
		self._Name = newObject._Name
		self._ra = float(newObject._ra)
		self._dec = float(newObject._dec)

	def cleanOutputFiles(self, output):
		# Cleans .coo files because IRAF does not like overwriting
                oldfiles = glob.glob(output)
                if oldfiles:
                        print "Found old .coo files"
                        print "Deleting..."
                        for old in oldfiles:
                                os.remove(old)
                                print "Removed: %s" % old
                return True
	
	def printInfo(self):
		print ""
		print "###############################"
		print "GENERAL INFO"
		print "###############################"
		print "Name: %s" % (self._Name)
		print ""
		print "Object ID: %s" % (self._id)
		print "Parent Image: %s" % (self._parentimage)
		print "WCS co-ordinates: %s %s" % (self._ra, self._dec)
		print "Pixel co-ordinates: %s %s" % (self._pixelx, self._pixely)
		print "Magnitude: %s" % (self._mag)
		print "Apperture Magnitude: %s +/- %s" % (self._appMag, self._appMagErr)
		print "Time: %s" % (self._midMJD)
		print "fap: %s" % (self._fap)
		print "fdan: %s" % (self._fdan)
		print "fan: %s" % (self._fan)
		print ""
		print "###############################"
		print "BOX STATISTICS"
		print "###############################"
		print "Box co-ordinates: %s" % self._Box
		print ""
		print "# of pixels: %s" % self._NPIX
		print "Mean: %s" % self._MEAN
		print "Median: %s" % self._MEDIAN
		print "STDDEV: %s" % self._STDDEV
		print "Min: %s" % self._MIN
		print "Max: %s" % self._MAX


	def distanceApart(self, imObject2):

		r = scipy.sqrt( (self._pixelx - imObject2._pixelx)**2 + (self._pixely - imObject2._pixely)**2 )
		return r

	def findLogicalPosition(self, write=False):

		# Convert the GRB position to the pixel x,y position
		# Make the input file for the IRAF task
		if self._ra == "" or self._dec == "":
			print "You did not set an ra or dec"
			return 1
		coordfilenamein = "%s" % (self._parentimage.replace(".fits", "_skyctran_in.coo"))
		coordfilenameout = "%s" % (self._parentimage.replace(".fits", "_skyctran_out.coo")) 

		try:
			self.cleanOutputFiles(coordfilenameout)
		except:
			print "Error:"
			print sys.exc_info()[0]

		try:
			coordout = open(coordfilenamein, "w")
			coordout.write("%10.7f %10.7f" % (self._ra, self._dec))
			coordout.close()
		except:
			print "File opening failed:"
			print "Did you reference ._parentimage???"
			print sys.exc_info()[0]

		# This subroutine lies in imcoords on IRAF
		iraf.image()
		iraf.imcoords()

		# Set the parameters of skyctran
		# 
		# input
		# output
		# insystem: input coordinate system
		# outsystem: output coordinate system
		# lngcolu: RA column no. 
		# latcolu: DEC column no.
		# ilnguni: RA input coord type
		# ilnlatuni: DEC input coord type
		# olnguni: RA output coord type
		# olatuni: DEC output coord type
		# ilngfor: RA input float format
		# ilatfor: DEC input float format
		# olngfor: RA output float format
		# olatfor: DEC output float format
		iraf.skyctran.setParam('input', coordfilenamein)
		iraf.skyctran.setParam('output', coordfilenameout)
		iraf.skyctran.setParam('insystem', '%s world' % self._parentimage)
		iraf.skyctran.setParam('outsystem', '%s logical' % self._parentimage)
		iraf.skyctran.setParam('lngcolu', '1')
		iraf.skyctran.setParam('latcolu', '2')
		iraf.skyctran.setParam('ilnguni', 'degrees')
		iraf.skyctran.setParam('ilatuni', 'degrees')
		iraf.skyctran.setParam('olnguni', 'logical')
		iraf.skyctran.setParam('olatuni', 'logical')
		iraf.skyctran.setParam('ilngfor', '%10.7f')
		iraf.skyctran.setParam('ilatfor', '%10.7f')
		iraf.skyctran.setParam('olngfor', '%10.3f')
		iraf.skyctran.setParam('olatfor', '%10.3f')
		iraf.skyctran.setParam('verbose', 'yes')

		# savepars
		skyctran = "skyctranpars.par"
#		iraf.skyctran.saveParList(skyctran)

		# run skyctran
		iraf.skyctran(mode='h', Stdout=1)
			
		# Set properties of this object from the output
		try:
			coordfile = open(coordfilenameout, "r")

			# parse the input file
			coords = coordfile.readlines()
			coordfile.close()
			for coord in coords:
				if coord[0] != "#" and coord[0] != "\n":
					coordinates = [i for i in coord.replace("\n","").split(" ") if i != ""]

			# Set the object parameters
			self._pixelx = float(coordinates[2])
			self._pixely = float(coordinates[3])

		except:
			print "File opening failed"
			print sys.exc_info()[0] 


		# region file for object of interest
		if write:
			try:
				regionout = open("%s" % (self._parentimage.replace(".fits", "_skyctran_object.reg")), "w")
                                regionout.write("image; circle(%f,%f,4) # color = red\n" % (self._pixelx, self._pixely))
                                regionout.close()
			except:
				print "Unexpected error:"
				print sys.exc_info()[0]

		self.printInfo()
		
	def findWorldPosition(self, write=False):

                # Convert the GRB pixel position to WCS
                # Make the input file for the IRAF task
                if self._pixelx == "" or self._pixely == "":
                        print "You did not set an pixel x or pixel y"
                        return 1
                coordfilenamein = "%s" % (self._parentimage.replace(".fits", "_skyctran_pixel_in.coo"))
                coordfilenameout = "%s" % (self._parentimage.replace(".fits", "_skyctran_pixel_out.coo"))

		try:
			self.cleanOutputFiles(coordfilenameout)
		except:
			print "Error cleaning files"
			print sys.exc_info()[0]

                try:
                        coordout = open(coordfilenamein, "w")
                        coordout.write("%10.7f %10.7f" % (self._pixelx, self._pixely))
                        coordout.close()
                except:
                        print "File opening failed:"
                        print sys.exc_info()[0]

                # This subroutine lies in imcoords on IRAF
                iraf.image()
                iraf.imcoords()

                # Set the parameters of skyctran
                # 
                # input
                # output
                # insystem: input coordinate system
                # outsystem: output coordinate system
                # lngcolu: RA column no. 
                # latcolu: DEC column no.
                # ilnguni: RA input coord type
                # ilnlatuni: DEC input coord type
                # olnguni: RA output coord type
                # olatuni: DEC output coord type
                # ilngfor: RA input float format
                # ilatfor: DEC input float format
                # olngfor: RA output float format
                # olatfor: DEC output float format
                iraf.skyctran.setParam('input', coordfilenamein)
                iraf.skyctran.setParam('output', coordfilenameout)
                iraf.skyctran.setParam('insystem', '%s logical' % self._parentimage)
                iraf.skyctran.setParam('outsystem', '%s world' % self._parentimage)
                iraf.skyctran.setParam('lngcolu', '1')
                iraf.skyctran.setParam('latcolu', '2')
                iraf.skyctran.setParam('ilnguni', 'logical') 
                iraf.skyctran.setParam('ilatuni', 'logical') 
                iraf.skyctran.setParam('olnguni', 'degrees')
                iraf.skyctran.setParam('olatuni', 'degrees')
                iraf.skyctran.setParam('ilngfor', '%10.7f')
                iraf.skyctran.setParam('ilatfor', '%10.7f')
                iraf.skyctran.setParam('olngfor', '%10.3f')
                iraf.skyctran.setParam('olatfor', '%10.3f')
                iraf.skyctran.setParam('verbose', 'yes')

                # savepars
                skyctran = "skyctranpars.par"
                #iraf.skyctran.saveParList(skyctran)

                # run skyctran
                tmp = iraf.skyctran(mode='h', Stdout=1) 

                # Set properties of this object from the output
                try:
                        coordfile = open(coordfilenameout, "r")

                        # parse the input file
                        coords = coordfile.readlines()
                        coordfile.close()
                        for coord in coords:
                                if coord[0] != "#" and coord[0] != "\n":
                                        coordinates = [i for i in coord.replace("\n","").split(" ") if i != ""]

                        # Set the object parameters
                        self._ra = float(coordinates[2])
                        self._dec = float(coordinates[3])

                except:
                        print "File opening failed"
                        print sys.exc_info()[0]


                # region file for object of interest
                if write:
                        try:
                                regionout = open("%s" % (self._parentimage.replace(".fits", "_skyctran_object.reg")), "w")
				if self._MEDFWHM == "":
					circ = 4.0
				else:
					circ = self._MEDFWHM*1.5
                                regionout.write("image; circle(%f,%f,%f) # color = red\n" % (self._pixelx, self._pixely, circ))
                                regionout.close()
                        except:
                                print "Unexpected error:"
                                print sys.exc_info()[0]


        def calculateError(self):

		# Airmass coefficients
		airmassCoefficients = {"g": 0.16, \
					"r": 0.09,\
					"i": 0.04,\
					"z": 0.04,\
					"J": 0.12,\
					"H": 0.06,\
					"K": 0.07\
					}

                # We assume it is running on a squared image
                print "######################################"
                print "INFO: YOU MUST FIRST SQUARE YOUR IMAGE"
                print "######################################"

                # Check the fluxes exist
		if self._skyFlux == "" or self._flux == "" or self._Band == "":
			print "Please calculate fluxes or specify the band"
			return 0

		# Calculation
                sqrd_error = self._skyflux + self._flux
		fluxerror = scipy.sqrt(sqrd_error)

		AIRMASS = getHeader("AIRMASS")
		EXPTIME = getHeader("EXPTIME")
		COEFF = airmassCoefficients[self._Band]

		fluxmoderror = zeropoint[self._ZP] + (1-AIRMASS)*COEFF + 2.5*scipy.log10(EXPTIME) - 2.5*scipy.log10(fluxerror)

		return fluxmoderror


class imFits(object):
      
   
	def __init__(self):
		self._Name = ""
		self._Band = ""
		self._logger = {}
		self._header = {}
		self._ObjectCounter = 0
		self._ObjectList = []
		self._MEDFWHM = ""
		self._skySTDEV = ""
		self._ZPDICT = {
				"g": 24.99, \
				"r": 25.21, \
				}

	def printInfo(self):
		print "Image: %s" % (self._Name)
		print "Band: %s" % (self._Band)
		print "# of objects: %s" % (self._ObjectCounter)
		print "Median FWHM: %s" % (self._MEDFWHM)
		print "Sky background STDDEV: %s" % (self._skySTDEV)

		# Auto functions
		# self.loadHeader()

	def giveObjectID(self):
		# return old value, but increment the counter
		tmp = self._ObjectCounter
		self._ObjectCounter = tmp + 1
		
		return tmp

	def loadHeader(self):
		hdulist = pyfits.open(self._Name)
		header = hdulist[0].header
		self._header = header
		hdulist.close()
		return self._header

	def getHeader(self, keyword='FITRAD', extra=False):
		try:
			if extra:
				headeritem = self._header.get(keyword, extra)
				return headeritem
			else:
				headeritem = self._header.get(keyword)
				return headeritem
		except:
			print "WARNING"
			sys.exc_info()[0]
			return "EMPTY"

	def getMidMJD(self, seconds=False, zero=False):
		self.loadHeader()
		MJD = float(self.getHeader("MJD-MID"))

		if seconds and zero:
			MJD = (MJD - zero)*60*60*24
		elif seconds:
			MJD = (MJD) * 60*60*24
		elif zero:
			MJD = MJD - zero

		return MJD

	def squareMyself(self):

		# Square the image
		newObjectName = self._Name.replace(".fits", "_sqrd.fits")

		# Check for output
		self.cleanOutputFiles(newObjectName)

		# Run iraf
		iraf.imarith(operand1=self._Name,
				operand2=self._Name,
				op="*",
				result=newObjectName
			)

		return newObjectName 

	def returnInfo(self, verbose=False):
		if verbose:
			print "Image name: _Name"
			
		return(self._Name)

	def makeObject(self, inObjectFile="mkobject.coo"):
		
		"""IRAF routine: mkobjects. This places a star in the fits image that is the parent object of this class, at the given ra and dec. It takes as input a file of "ra dec magnitude" for a star type object.

		1. Get info from the header files if you want noise
		2. Get FWHM for size of the star for the given image
		3. Generate the star

		radius = FWHM / scale, from statistics, scale = 2*sqrt(2*ln(2))
		"""
	
		# 1.
		self.loadHeader()
		INFODICT = {
			"RON": self.getHeader("RON"), \
			"GAIN": self.getHeader("GAIN"), \
			"EXPTIME": self.getHeader("EXPTIME"), \
			"MAGZP": self._ZPDICT[self.getHeader("FILTER")], \
			"POISSON": "no", \
			"STAR": "gaussian", \
			}	
		# 2.
		scale = 2*scipy.sqrt(2*scipy.log(2))
		self.getMyMedianFWHM()
		INFODICT["RADIUS"] = self._MEDFWHM / scale

		# 3. 
		outputname = self._Name.replace(".fits", "_mko.fits")
		iraf.noao(_doprint=0)
		iraf.artdata(_doprint=0)
		
	        iraf.mkobjects(input=self._Name, \
                        output=outputname, \
                        # When creating new images
                        #title="",\
                        #ncols=self.xdim, \
                        #nlines=self.ydim, \
                        #header=self.hdrfile, \
                        #background=0.0, \
                        # When creating objects
                        objects=inObjectFile, \
                        xoffset=0.0, \
                        yoffset=0.0, \
                        star=INFODICT["STAR"], \
                        radius=INFODICT["RADIUS"], # this is fwhm/scale (pixels), where scale = 2sqrt(2ln2)\
                        #beta=2.5,  # for star=moffat profile\
                        ar=1.0, \
                        pa=0.0, \
                        distance=1.0, \
                        exptime=INFODICT["EXPTIME"], \
                        magzero=INFODICT["MAGZP"], \
                        # Noise parameters
                        gain=INFODICT["GAIN"], \
                        rdnoise=INFODICT["RON"], \
                        poisson=INFODICT["POISSON"], \
                        seed=1, \
                        comments=1)

		return outputname 

	def reMap(self, remapFits, outname="input_remapped.fits", verbose=False, wcsregister=False):
	  
		# Create an object class for the output image for fluency of continuing code
		outFits = imFits()
		outFits._Name = outname
		outFits._logger["ReMapping"] = []
	  
		# remapFits = template image
		# this object = source image
		# outFits = output fits image
		
		if wcsregister:
			iraf.wregister(input=self._Name, output=outFits._Name, wcs="world", reference=remapFits._Name, Stdout=1)

		else:
			wcscmd = ["wcsremap","-template", self._Name, "-source", remapFits._Name, "-outIm", outFits._Name]
			wcsremap = subprocess.Popen(wcscmd, stdout=subprocess.PIPE)
			outFits._logger["ReMapping"].append(wcsremap.communicate()[0])
		
		if verbose:
			for log in outFits._logger["ReMapping"]:
				print log
		
		return outFits

	def trimMyself(self, outname="remcut.fits", region="[400:500,400:500]", verbose=False):
		
                self._logger["trimMyself"] = []

		if self._Name == "None":
			self._logger["trimMyself"].append("No filename set, please check.")

		imCopy = imFits()
		imCopy._Name = outname
		iraf.imcopy("%s%s" % (self._Name,region), imCopy._Name)
		
		if verbose:
			self._logger["trimMyself"].append(iraf.imstat(self._Name))	

			for log in self._logger["trimMyself"]:
				print log

		return imCopy

	def subtractTemplate(self, templateFits, outname, _tu="None", _tuk="None", _tl="None", _tg="None", _tr="None", _iu="None", _iuk="None", _il="None", _ig="None", _ir="None", _nsx="None", _nsy="None", _ng="None", band="r", verbose=False, _norm="None"):
		
		if self._Band == "None":
			print "You must set the band"
			return 1

		# Logger
                self._logger["hotpants"] = []

		# Output image object
		outFits = imFits()
		outFits._Name = outname
		outFits._Band = band
		outFits._Noise = outname.replace('.fits', '_noise.fits')

		# Band dictionary for gain and readout noise
		self.loadHeader()
		readgain_dict = {\
					"GAIN": self.getHeader("GAIN"), \
					"RON": self.getHeader("RON"),\
				}

		# Will fill the correct details for each band
		if _tg == "None":
			_tg = readgain_dict["GAIN"]
		if _tr == "None":
			_tr = readgain_dict["RON"]
		if _ig == "None":
			_ig = readgain_dict["GAIN"]
		if _ir == "None":
			_ir = readgain_dict["RON"]

		# Default values from Wiki page if no input values
		if _tu == "None":
			_tu = 25000
		if _tuk == "None":
			_tuk = 25000
		if _tl == "None":
			_tl = -80
		if _iu == "None":
			_iu = 15000
		if _iuk == "None":
			_iuk = 15000
		if _il == "None":
			_il = -80
		if _nsx == "None":
			_nsx = 20
		if _nsy == "None":
			_nsy = 20
		if _ng == "None":
			_ng = [3,6,0.6,4,1.5,2,3]
		if _norm == "None":
			_norm = "i"

		# [-inim fitsfile] : comparison image to be differenced
		# [-tmplim fitsfile]: template image
		# [-outim fitsfile] : output difference image
		# [-tu tuthresh] : upper valid data count, template (25000)
		# [-tuk tucthresh] : upper valid data count for kernel, template (tuthresh)
		# [-tl tlthresh] : lower valid data count, template (0)
		# [-tg tgain] : gain in template (1)
		# [-tr trdnoise] : e- readnoise in template (0)
		# [-iu iuthresh] : upper valid data count, image (25000)
		# [-iuk iucthresh] : upper valid data count for kernel, image (iuthresh)
		# [-il ilthresh] : lower valid data count, image (0)
		# [-ig igain] : gain in image (1)
		# [-ir irdnoise] : e- readnoise in image (0)
		# [-nsx xstamp] : number of each region's stamps in x dimension (10)
		# [-nsy xstamp] : number of each region's stamps in y dimension (10)
		# [-ng ngauss degree0 sigma0 .. degreeN sigmaN]
		
		#    * : ngauss = number of gaussians which compose kernel (3) : degree = degree of polynomial associated with gaussian #
		#          o (6 4 2) 
		#      : sigma = width of gaussian #
		#          o (0.70 1.50 3.00) 
		#      : N = 0 .. ngauss - 1 

		# Call hotpants
		cmd = ["%s/hotpants" % HOTPANTSPATH, "-v", "0", "-inim", self._Name, "-tmplim", templateFits._Name, "-outim", outFits._Name, "-oni", outFits._Noise, "-tu", "%s" % _tu, "-tuk", "%s" %  _tuk, "-tl", "%s" % _tl, "-tg", "%s" % _tg, "-tr", "%s" % _tr, "-iu", "%s" % _iu, "-iuk", "%s" % _iuk, "-il", "%s" % _il, "-ig", "%s" % _ig, "-ir", "%s" % _ir, "-nsx", "%s" % _nsx, "-nsy", "%s" % _nsy, "-n", "%s" % _norm]

		cmd.append("-ng")

		for inputc in _ng:
			cmd.append("%s" % inputc)

		#cmd.append(">")
		#cmd.append("/dev/null")

		print cmd
		tmp = ""
		for i in cmd:
			tmp = "%s %s" % (tmp, i)
		print tmp
		      
		hotpants = subprocess.Popen(cmd).wait() #, stdout=subprocess.PIPE)
		self._logger["hotpants"].append("%s" % hotpants.communicate()[0])

		try:
			os.kill(hotpants.pid, signal.SIGKILL)
		except:
			print "Process killed normally"

		for log in self._logger["hotpants"]:
			print log

		return outFits 

	def cleanOutputFiles(self, output):
		# Cleans .coo files because IRAF does not like overwriting
                oldfiles = glob.glob(output)
                if oldfiles:
                        print "Found old .coo files"
                        print "Deleting..."
                        for old in oldfiles:
                                os.remove(old)
                                print "Removed: %s" % old
		return True


	def findStars(self):

		# To look for the stars that exist in the image specified.
		# It will:
		#	Utilise the IRAF package daofind
		#	be used for this image subtraction package
		#	Gives ouput of co-ords, "#imagename_daofind.coo"
		output = self._Name.replace(".fits", "_daofind.coo")
		try:
			self.cleanOutputFiles(output)
		except:
			print "Error cleaning output files"
			print sys.exc_info()[0]

		# Load as in IRAF and neglect output:
		# NOAO
		# DIGIPHOT
		# APPHOT
		# DAOFIND
		iraf.noao(_doprint=0)
		iraf.digiphot(_doprint=0)
		iraf.apphot(_doprint=0)

		# Set the parameter list
		# image: input image
		# output: output co-ordinates
		# verify: don't know
		# threshold: sigma above which stars are determined (fixed in this script)
		# datapar: point to a file with the user input
		#
		#	ccdread: CCD readout noise (image header keyword)
		#	gain: CCD gain (image header keyword)
		#	readnoi: CCD readout noise in e
		#	epadu: gain in e- per count
		#
		# findpar: point to a file with the user input
		#
		#	threshold: in units of sigma
		#
		# centerpars:
		#
		# Save the parameter list

		# Pull information from the header of the FITS image
		#self.loadHeader()
		#datapars = iraf.datapars.getParList()
		#findpars = iraf.findpars.getParList()
		#ceneterpars = iraf.centerpars.getParList()
		#for par in datapars:
		#print par
		## Print to user
		#print "Information on images:"
		#print "Gain: %s" % self.getHeader('RON')
		#print "RON: %s" % self.getHeader('RON')
	
		# datapar:
		iraf.datapars.setParam('ccdread', 'RON')
		iraf.datapars.setParam('gain', 'GAIN')
		iraf.datapars.setParam('exposur','EXPTIME')
		iraf.datapars.setParam('airmass','AIRMASS')
		iraf.datapars.setParam('filter','SUBSET')
		iraf.datapars.setParam('obstime','HIERARCH ESO OBS START')
		iraf.datapars.setParam('datamin', '-100')
		iraf.datapars.setParam('datamax', '60000')
#		iraf.datapars.saveParList(filename="%s/uparm/datapars.par" % os.getcwd())

		# findpar:
		#threshold = 4.0
		#iraf.findpars.setParam('threshold', str(threshold))
		#iraf.findpars.saveParList(filename="uparm/findpars.par")

		#print glob.glob(output)
			

		# daofind:
		iraf.daofind(image=self._Name,
				output=output,
				verify=0,
				verbose=1,
				sigma=self._skySTDEV,
				nsigma=1,
				scale=1,
				fwhmpsf=self._MEDFWHM,
				sharplo='0.0',
				sharphi='1.0',
				roundlo='-1.0',
				roundhi='1.0',
				thresho='4.0',
				wcsout="logical",
				Stdout=1,
				mode='h'
			)


		# 1. Open the .coo file
		# 2. Pasre the input

		# 1. take input file
		filename = "%s" % (self._Name.replace(".fits", "_daofind.coo"))
		coords = open(filename)
		coord = coords.readlines()
		coords.close()

		# 2. parse out the lines beginning with "#"
		filename = "%s" % (self._Name.replace(".fits", "_daofind.reg"))
		coordregfile = open(filename, "w")
		coobject = []
		for co in coord:
			if co[0] != "#":

		# 2b. parse each element
			
				newCoordObject = imObject()
				coitem = [i for i in co.replace("\n","").split(" ") if i!= ""]
				newCoordObject._id = self.giveObjectID()
				newCoordObject._pixelx = float(coitem[0])
				newCoordObject._pixely = float(coitem[1])
				newCoordObject._mag = float(coitem[2])
				#newCoordObject.printInfo()

				coordregfile.write("image; circle(%f,%f,4) # color = white\n" % (newCoordObject._pixelx, newCoordObject._pixely))

				coobject.append(newCoordObject)

		coordregfile.close()
		self._ObjectList = coobject
				
	def findNearestNeighbours(self, CoordObject, rcirc=100.0, write=False):

		# temporary type check
		#print type(CoordObject._pixelx), type(CoordObject._pixely)

		# 1. Look at object list
		#	look for other objects closest to CoordObject
		dist = 100000
		NearestNeighbours = []
		for coord in self._ObjectList:

			# Not itself
			dist = CoordObject.distanceApart(coord)
			if CoordObject._id != coord._id and dist > rcirc and dist < rcirc+100:
				#print "Comparing: %s and %s" % (CoordObject._id, coord._id)
				#print "Distance: %s" % dist
				NearestNeighbours.append([coord,dist])

	
		if write:
			NNOutput = open(self._Name.replace(".fits", "_nn.reg"), "w")
			for i in NearestNeighbours:
				NNOutput.write("image; circle(%s,%s,10) # color = blue\n" % (i[0]._pixelx, i[0]._pixely))
			NNOutput.close()
				
	
		return NearestNeighbours

	def findNearestNeighbour(self, CoordObject, rcirc=100.0, write=False):

		# NearestNeighbours is a list of [imObject, dist]

		NearestNeighbours = self.findNearestNeighbours(CoordObject=CoordObject, rcirc=rcirc, write=write)

		# Find the smallest distance

		distance = 10000
		smallest = []
		for Neighbour in NearestNeighbours:
			if distance > Neighbour[1]:
			
				smallest = Neighbour
				distance = Neighbour[1]

		return smallest
			
	def runApperturePhotometryOnObject(self, photObject, write=False):

		# Aperture sizes
                #fap = 1.0
                #fdan = 2.0
                #fan = 3.0

		# To take photometry of a given object using DAOPHOT/APPHOT.PHOT

		# Check we have an object list
		coordfilename = self._Name.replace(".fits", "objectID_%s.coo" % photObject._id)
		coordfilename2 = self._Name.replace(".fits", "_apphot_%s.reg" % photObject._Name.replace(" ",""))

		try:
			regionOut = open(coordfilename2, "w")
			regionOut.write("image; circle(%s, %s, %s) # color = red\n"% (photObject._pixelx, photObject._pixely, photObject._fap))
			regionOut.write("image; circle(%s, %s, %s) # color = blue\n" % (photObject._pixelx, photObject._pixely, photObject._fdan))
			regionOut.write("image; circle(%s, %s, %s) # color = green\n" % (photObject._pixelx, photObject._pixely, photObject._fan))
			regionOut.close()
		except:
			print "ERROR opening region file to write"
			print sys.exc_info()[0]

		try:
			objectList = open(coordfilename,"r")
			objectList.close()
					
		except:
			objectList = open(coordfilename,"w")
			objectList.write("%s %s" % (photObject._pixelx, photObject._pixely))
			objectList.close()

                # Load as in IRAF and neglect output:
                # NOAO
                # DIGIPHOT
                # APPHOT
                # DAOFIND
                iraf.noao(_doprint=0)
                iraf.digiphot(_doprint=0)
                iraf.apphot(_doprint=0)		

		# DAOPHOT/APPHOT - PHOT Package
		# Parameters:
		#
		# image
		# skyfile???
		# coords - input file
		# output - output file
		# datapars
		# centerpars
		# fitskypars
		# photpars

		cwd = os.getcwd()
		# datapars:
		datapars = "%s/uparm/datapars.par" % cwd
		#datapars = "uparm/datapars.par"
		#iraf.datapars.saveParList(filename=datapars)
		# centerpars:
		
		centerpars = "%s/uparm/centerpars.par" % cwd
		#centerpars = "uparm/centerpars.par"
#		iraf.centerpars.saveParList(filename=centerpars)
		# fitskypars:
		fitskypars = "%s/uparm/fitskypars.par" % cwd
		#fitskypars = "uparm/fitskypars.par"
#		iraf.fitskypars.saveParList(filename=fitskypars)
		# photpars:
		photpars = "%s/uparm/photpars.par" % cwd
		#photpars = "uparm/photpars.par"
#		iraf.photpar.saveParList(filename=photpars)


		#iraf.datapar(sigma=1.5, exposure='IMGEXP', gain='GAIN', ccdread='RON')

		# standard
		#iraf.phot.setParam('image', self._Name)
		#iraf.phot.setParam("output", self._Name.replace(".fits", "_phot.mag"))
		#iraf.phot.setParam("datapar", datapars)
		#iraf.phot.setParam("centerp", centerpars)
		#iraf.phot.setParam("fitskyp", fitskypars)
		#iraf.phot.setParam("photpar", photpars)
		#iraf.phot.setParam("interac", "no")
		#iraf.phot.setParam("verify", "yes")

		# co-ordinate file
		iraf.phot.setParam("coords", coordfilename)

		# run phot
		apphotpar = "%s/uparm/phot.par" % cwd
#		iraf.phot.saveParList(filename=apphotpar)

		output = self._Name.replace(".fits", "_objectID%s_phot.mag" % photObject._id)
		try:
			self.cleanOutputFiles(output)
		except:
			print "Error cleaning failed"
			print sys.exc_info()[0]

		#apsizes= (.4*fap*self._MEDFWHM,.5*fap*self._MEDFWHM,.6*fap*self._MEDFWHM,.8*fap*self._MEDFWHM,
	        #        	1.*fap*self._MEDFWHM,1.2*fap*self._MEDFWHM, 1.5*fap*self._MEDFWHM,2.*fap*self._MEDFWHM,
	        #		2.5*fap*self._MEDFWHM,3.*fap*self._MEDFWHM
		#		)
		#irafapsizes = '%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f' % apsizes
		#irafapsizes = '%.2f' % (5*apsizes[4])

		self.loadHeader()
		if self._Band in ["J", "H", "K"]:
			ifilter = self.getHeader("FILTER")
		else:
			ifilter = self.getHeader("SUBSET")

                print "RON: %s" % self.getHeader("RON")
                print "GAIN: %s" % self.getHeader("GAIN")
                print "EXPTIME: %s" % self.getHeader("EXPTIME")
                print "AIRMASS: %s" % self.getHeader("AIRMASS")
                print "DATE: %s" % self.getHeader("DATE-OBS")

		iraf.phot(image=self._Name,
				coords=coordfilename, 
				output=output,
				verify=0,
				verbose=1,
				interac="no",
				scale=1,
				fwhmpsf=self._MEDFWHM,
				sigma=self._skySTDEV,
				wcsin="logical",
				wcsout="logical",
				datamin=-100,
				datamax="INDEF",
				zmag=22.703,	# to do: make function to calculate
				annulus=photObject._fan,
				dannulus=photObject._fdan,
				calgorithm="none",
				aperture=photObject._fap,
				cbox=1.5*self._MEDFWHM,
				maxshift=15,
				mode="ql",
				Stdout=1,
		 		readnoi=self.getHeader("RON"),
				epadu=(2/3)*self.getHeader("NIMGS")*self.getHeader("GAIN"),
				itime=self.getHeader("EXPTIME"),
				xairmass=self.getHeader("AIRMASS"),
				ifilter=ifilter,
				otime=self.getHeader("DATE-OBS"),
				salgori="mode",
				#skyvalu=0,
				smaxite=1
				)

		# Parse the output
		photout = open(output, "r")
		photoObjectProperties = photout.readlines() 
		photout.close()
		
		propList = []
		for Properties in photoObjectProperties:
			if Properties[0] != "#":
				Property = [i for i in Properties.replace("\n","").replace("\\","").split(" ") if i != ""]
				propList.append(Property)

		photObject._skyFlux = propList[3][0]
		photObject._flux = propList[4][1]
		photObject._appMag = propList[4][4]
		photObject._appMagErr = propList[4][5]
		photObject._midMJD = self.getMidMJD()#seconds=True,zero=55822.89371528)

	def getMyMedianFWHM(self):

		# Create a SExtractor instance
		sex = sextractor.SExtractor()

	        # Modify the SExtractor configuration
        	sex.config['DETECT_MINAREA'] = 8
	        sex.config['DETECT_THRESH'] = 8
		#MAXAREA is not known in version on faramir
	        #sex.config['DETECT_MAXAREA'] = 100

	        # Add a parameter to the parameter list
	        sex.config['PARAMETERS_LIST'].append('FWHM_IMAGE')

	        # Lauch SExtractor on a FITS file
	        sex.run(self._Name)

	        # Print FHWM information
	        catalog_name = sex.config['CATALOG_NAME']
	        catalog_f = sextractor.open(catalog_name)
	        catalog = catalog_f.readlines()
	        FWHMList = numpy.array([])

		catalogooname = "%s" % self._Name.replace(".fits", "_sex.reg")
	        catalogoo = open(catalogooname, "w")

	        for star in catalog:
	                catalogoo.write("image; circle(%f,%f,%f) # color = green\n" % (star['X_IMAGE'], star['Y_IMAGE'], 10))
	                FWHMList = numpy.append(FWHMList,float(star['FWHM_IMAGE']))

	        catalogoo.close()

	        # Take the smallest 20% of the FWHM
	        FWHMList = numpy.sort(FWHMList)
	        FWHMListLength = len(FWHMList)
	        FWHMTwenty = int(0.2*FWHMListLength)
	        FWHMTwentyList = numpy.array([])

	        for FWHM in range(FWHMTwenty):

	                FWHMTwentyList = numpy.append(FWHMTwentyList, FWHMList[FWHM])

	        medianfwhm = numpy.median(FWHMTwentyList)
		print "Median FWHM: %f" % medianfwhm
	
	        # Plot
	        fig = plt.figure(0)
	        ax = fig.add_subplot(111)
	        ax.plot(range(len(FWHMList)),FWHMList)
	        ax.set_xlabel('star')
	        ax.set_ylabel('FHWM/pixels')
	        savename = "%s" % (self._Name.replace(".fits", "_sex.png"))
	        plt.savefig(savename, format="png")

	        sex.clean(config=True, catalog=True, check=True)

		self._MEDFWHM = float(medianfwhm)

	def getBackgroundSTDEV(self):

		# load the header
		self.loadHeader()

		# Use iraf to get the sky background STDEV
	        skySTDEV = iraf.imstat(images=self._Name,
					fields='stddev',
					nclip = 25,
					format=0,
					Stdout=1
				)

		self._skySTDEV = float(skySTDEV[0])

		# Skynoise in JHK is interpolation smoothed due to blowing up the images
		# In the Drizzle algorithm the ratio of true to output noise is R = 1/(1-r/3) where r is pixfrac/scale

		if self._Band in 'JHK':
			self._skySTDEV = self._skySTDEV * self.getHeader("INTERPSM", 1.662)
		else:
			self._skySTDEV = self._skySTDEV * self.getHeader("INTERPSM", 1.07)

		print "SKY STDDEV: %f" % (self._skySTDEV)

		return self._skySTDEV

	def getObjectStatistics(self, photObject, multiple=2.5):

		#photObject.printInfo()

		# Will get statistics around a box of width = w + multiple/2., height = h + multiple/2.
		_boxX1 = int(photObject._pixelx - multiple*(self._MEDFWHM/2.))
		_boxX2 = int(photObject._pixelx + multiple*(self._MEDFWHM/2.))

		_boxY1 = int(photObject._pixely - multiple*(self._MEDFWHM/2.))
		_boxY2 = int(photObject._pixely + multiple*(self._MEDFWHM/2.))

		irafinput = "%s[%d:%d,%d:%d]" % (self._Name, _boxX1, _boxX2, _boxY1, _boxY2)

		tmp = iraf.imstat(irafinput, 
					fields='npix,mean,stddev,min,max,midpt',
					format=0,
					Stdout=1
				)

		tmp = [i for i in tmp[0].split(" ") if i != ""]
		photObject.setBox([_boxX1,_boxX2,_boxY1,_boxY2])
		photObject.setNPIX(float(tmp[0]))
		photObject.setMEAN(float(tmp[1]))
		photObject.setSTDDEV(float(tmp[2]))
		photObject.setMIN(float(tmp[3]))
		photObject.setMAX(float(tmp[4]))
		photObject.setMEDIAN(float(tmp[5]))
		
if __name__ == "__main__":
  
	print Usage
