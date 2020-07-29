import requests
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
from scipy.interpolate import UnivariateSpline
from astropy.time import Time

def getBbandinfo(jsonboi, allinfo):
	try:
		jsonupack = json.loads(jsonboi)
		for f in jsonupack:
			photometry = jsonupack[f]['photometry']
			if len(photometry) > 10:
				allinfo.append(jsonupack[f])
		return allinfo
	except:
		return allinfo

def findWavelengthfromvel(v, lam_r, c=2.998E5):
    lilboi = math.sqrt(((1+(v/c))/(1-(v/c)))*lam_r*lam_r)
    return lilboi


def trimdata(data, avgflux, xlow, xhi):
	waves = []
	flux = []
	for x in data:
		a = float(x[0])
		if a > xlow and a < xhi:
			try:
				flux.append(math.log10(float(x[1])))
				waves.append(float(x[0]))
			except Exception as e:
				pass
	thisavg = np.mean(flux)
	offset = avgflux-thisavg-0.8
	flux = [f+offset for f in flux]
	return waves, flux

class sndata:
	def __init__(self, name, sntype, absmag, absmag_e, sivel, sivel_e):
		self.name = name
		self.sntype = sntype
		self.absmag = absmag
		self.absmag_e = absmag_e
		self.sivel = sivel
		self.sivel_e = sivel_e

def main5():
	font = {'family': 'serif',
	        'color':  'black',
	        'weight': 'normal',
	        'size': 10}

	zhengdata = []
	zhengfile = 'zhengdata.txt'
	zheng_o = open(zhengfile).readlines()[1:]

	zdata = []
	for row in zheng_o:
		rs = row.split()

		if len(rs) == 20:
			zdata.append(sndata(rs[0], rs[1], float(rs[14]), float(rs[16].split(',')[0]), float(rs[7]), float(rs[9].split(',')[0])))
		elif len(rs) == 22:
			zdata.append(sndata(rs[0], rs[1], float(rs[15+1]), float(rs[17+1].split(',')[0]), float(rs[7+1]), float(rs[9+1].split(',')[0])))
		else:
			zdata.append(sndata(rs[0], rs[1], float(rs[15]), float(rs[17].split(',')[0]), float(rs[7]), float(rs[9].split(',')[0])))

	#Loop over each of them,
	#query photometry
	#plot their B band magnitude (I know that it exists) Measure dm15
	#plot NIR light curves relative to B band.
	#in second plot, show their early optical spectra

	photbands = ["B", "I", "i", "Y"]

	for sninfo in zdata:

		sn = sninfo.name
		r = requests.get('https://api.astrocats.space/'+sn)
		j = json.loads(r.content)

		try:
			figsize = plt.figaspect(1.0/0.8)

			fig, axes = plt.subplots(1, 1, figsize=figsize)

			specaxC = axes

			maxdate = j[sn]['maxdate'][0]['value'].split('/')
			t = Time([maxdate[0]+'-'+maxdate[1]+'-'+maxdate[2]])
			maxmjd = t.mjd[0]

			photoffset = 0
			newmaxmjd = 0
			for b in photbands:
				photrequest = requests.get('https://api.astrocats.space/'+sn+'/photometry/time+magnitude+band?band='+b+'')
				photjson = json.loads(photrequest.content)
				photometry = photjson[sn]['photometry']

				thresh = 25
				if len(photometry) > 7:
					dates = [float(x[0]) for x in photometry]
					mags = [float(x[1]) for x in photometry]
					if b == 'B':
						cutdates = [i for f,i in enumerate([d for d in dates if d < maxmjd + thresh])]
						cutmags = [mags[f] for f,i in enumerate([d for d in dates if d < maxmjd + thresh])]
						xinterp = np.linspace(min(cutdates),maxmjd+thresh, 1000)
						#try:
						spl = UnivariateSpline(cutdates, cutmags, k=5)
						fit = spl(xinterp)

						maxmag = min(fit)
						index = list(fit).index(maxmag)
						newmaxmjd = xinterp[index]

						day15offset = [abs(newmaxmjd+15-x) for x in xinterp]
						day15index = list(day15offset).index(min(day15offset))
						day15mag = fit[day15index]

						mb15 = abs(maxmag-day15mag)
					else:
						#get nir data
						dmag = maxmag - min(mags)
						photoffset = dmag - 1.0
						maxmag = min(mags)

						tdates = [x - newmaxmjd for x in dates]
						tmags = [x - photoffset for x in mags]

						pass

			avgfluxC = 0
			specrequest = requests.get('https://api.astrocats.space/'+sn+'/spectra/time+data')
			specjson = json.loads(specrequest.content)
			spectra = specjson[sn]['spectra']
			if len(spectra) > 0:
				if any(float(s[0])-newmaxmjd <=-6 for s in spectra):
					ss = [s for s in spectra if float(s[0]) -newmaxmjd <= 0]
					for s in ss:
						acqdate = s[0]
						if float(acqdate)-newmaxmjd <= 0:
							specdata = s[1]
							waves, flux = trimdata(specdata, avgfluxC, 6000, 7000)
							waves = [x/10000.0 for x in waves]
							avgfluxC = np.mean(flux)
							dt = str(round(float(acqdate)-newmaxmjd,2))
							specaxC.plot(waves, flux, '-k')
							specaxC.annotate(dt, xy=(0.6500, avgfluxC+0.15), fontsize=9)

			ylabel = '$log_{10}(Flux)$+constant'
			#specax.axvline(x=findWavelengthfromvel(-12000, 0.6588), linewidth=6, color='red', alpha=0.33)
			specaxC.yaxis.set_major_formatter(plt.NullFormatter())
			specaxC.yaxis.set_ticks_position('none')
			for tick in specaxC.xaxis.get_major_ticks():
				tick.label.set_fontsize(8)
			
			plt.subplots_adjust(wspace=0, hspace=0)
			plt.locator_params(numticks=12)
			fig.text(0.5, 0.93, sn, ha='center', fontdict=font)
			plt.show()

			hasC = input('Is There Carbon?: [y/n/m]')
			zhengdata.append([sninfo.name, sninfo.sntype, sninfo.absmag, sninfo.absmag_e, sninfo.sivel, sninfo.sivel_e, round(mb15, 2), hasC])
		except:
			print('Object not in sne.space')
			hasC = 'nd'
			zhengdata.append([sninfo.name, sninfo.sntype, sninfo.absmag, sninfo.absmag_e, sninfo.sivel, sninfo.sivel_e, round(mb15, 2), hasC])

	of = open('zheng_search.txt', 'w')
	for row in zhengdata:
		of.write('{}, {}, {}, {}, {}\n'.format(str(row[0]), str(row[1]), str(row[2]), str(row[3]), str(row[4]), str(row[5]), str(row[6]), str(row[7])))
	of.close()
	print(zhengdata)

main5()
