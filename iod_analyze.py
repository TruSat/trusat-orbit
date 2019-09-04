#!/usr/bin/env python3
import sys
import os # Only needed for paths while our repos are not public

if sys.version_info[0] != 3 or sys.version_info[1] < 7:
	print("This script requires Python version 3.7")
	sys.exit(1)

import os, re, math, csv
from skyfield.api import Topos
from skyfield.api import load

import inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
iod_path = os.path.join(parentdir, "sathunt-iod")
sys.path.insert(1,iod_path) 
from iod import iod_parse_line

import logging
log = logging.getLogger(__name__)

class CosparSite:
	# Class variable

	def __init__(self,sitefilename="data/sites.txt"):
		self.sitefilename = sitefilename
		self.site_re = re.compile('^(\d{4})\s(\w{2})\s+([+-]*\d+\.\d+)\s+([+-]*\d+\.\d+)\s+([+-]*\d+)\.*\s+(\w+.*$)')

		# Read COSPAR stations (hack)
		self.site = {}
		
		self.sitefile = open(self.sitefilename, "r")
		for self.line in self.sitefile:
			self.match = self.site_re.search(self.line)
			if self.match:
				self.site[int(self.match[1])] = {
					'id' 	    : self.match[2],
					'latitude'  : float(self.match[3]),
					'longitude' : float(self.match[4]),
					'elevation' : int(self.match[5]),
					'observer'  : self.match[6].rstrip()
				}
		self.sitefile.close()

	"""Function to output Topos tuple for requested COSPAR station ID"""
	def topos(self,cospar=9999):
		try:
			return ( [self.site[cospar]['latitude'],self.site[cospar]['longitude'], self.site[cospar]['elevation']] )
		except:
			log.error("no site data for {:d}".format(cospar))
			sys.exit()

# MAIN
def main():

	#planets = load('de421.bsp')
	#earth = planets['earth']

	Site = CosparSite("data/sites.txt")
	duvall = Topos('47.77166666666667 N', '121.90416666666667 W')



	ts = load.timescale()

	#t = ts.now()

	#stations_url = 'http://celestrak.com/NORAD/elements/stations.txt'
	satellites_file = '../tle/bulk.tle'
	satellites = load.tle(satellites_file)

	#satellite = satellites['ARKYD 6A']
	satellite = satellites['TDRS 11']


	#astrometric = observer_location.at(t).observe(mars)
	#alt, az, d = astrometric.apparent().altaz()

	#print(" AZ: ",az)
	#print("ALT: ",alt)

	for line in sys.stdin:
		IOD_line = iod_parse_line(line)
		if (IOD_line):
			try:
				satellite = satellites[IOD_line.ObjectNumber]

				tm = IOD_line.DateTime
				tmseconds = tm.second + (tm.microsecond/1000000)

				t = ts.utc(tm.year,tm.month,tm.day,tm.hour,tm.minute,tmseconds)

				(lat, lon, elevation_meters) = Site.topos(IOD_line.Station)
				observer_location = Topos(latitude_degrees = lat, longitude_degrees = lon, elevation_m = elevation_meters)

				#print(IOD_line.line)
				#print(satellite.name)

				days = t - satellite.epoch
				#print('{:.3f} days away from epoch'.format(days))

				difference = satellite - observer_location
				#print(difference)

				topocentric = difference.at(t)
		#		print(topocentric.speed().km_per_s)
		#		print(topocentric.position.km)
	
				# (dx,dy,dz) = topocentric.position.km
				# print(dx,dy,dz)

				print("Station {} Observing Satellite {} : {:.3f} days away from epoch".format(IOD_line.Station,satellite.model.satnum,days))
				if(IOD_line.RA is not None):
					ra, dec, distance = topocentric.radec()  # ICRF ("J2000")
					print("RA:  {:7.4f} Obs: {:.4f} (delta: {:.3f})".format(ra._degrees, IOD_line.RA, (IOD_line.RA-ra._degrees) ))
					print("DEC: {:7.4f} Obs: {:.4f} (delta: {:.3f})".format(dec._degrees, IOD_line.DEC, (IOD_line.DEC-dec._degrees) ))

					# An attempt at computing residuals, based on sattools residuals.c
					# rx0 = ra._degrees
					# ry0 = dec._degrees

					# drx = IOD_line.RA  - rx0
					# dry = IOD_line.DEC - ry0

					# dt = -(rx0*drx+ry0*dry)/(drx*drx*dry*dry)
					# dr = math.sqrt( ((dry*rx0-drx*ry0)**2.0)/(drx*drx+dry*dry) )

					# if ((-rx0*drx-ry0*dry)<0.0):
					#  	dr *=-1
					# print("Residuals: ",dr,dt)
				else:
					alt, az, distance = topocentric.altaz()
					print(az.degrees, IOD_line.AZ)
					print(alt.degrees, IOD_line.EL)
				print()
			except:
				print("### Object not in catalog - new observation for {} ###".format(IOD_line.ObjectNumber))

if (__name__ == '__main__'):
  main()