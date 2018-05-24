# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
distributor
===========

A simple Flask_ application to serve spectral data.

.. _Flask: http://flask.pocoo.org/
"""
import os
import glob
import re
import json

from flask import Flask, jsonify
from werkzeug.contrib.cache import SimpleCache

import numpy as np
from astropy.io import fits


app = Flask(__name__)
cache = None
healpixels = None


@app.before_first_request
def load_app_data():
    """Get the list of spectra-64 files and initialize the cache.
    """
    global cache, healpixels
    cache = SimpleCache()
    desi_spectro_redux = os.environ['DESI_SPECTRO_REDUX']
    specprod = os.environ['SPECPROD']
    spectra = glob.glob(os.path.join(os.environ['DESI_SPECTRO_REDUX'],
                                     os.environ['SPECPROD'],
                                     'spectra-64',
                                     '[0-9][0-9]',
                                     '*',
                                     'spectra-64-*.fits'))
    spectrare = re.compile(r'spectra-64-([0-9]+)\.fits')
    healpixels = dict([(int(spectrare.match(os.path.basename(s)).groups()[0]), s) for s in spectra])


@app.route('/')
def index():
    """Return all healpixels.
    """
    response = jsonify(sorted(list(healpixels.keys())))
    response.headers.add('Access-Control-Allow-Origin', '*')
    return response


@app.route('/<int:pixel>')
def spectra(pixel):
    """Return the number of spectra in a healpixel.
    """
    response = cache.get(str(pixel))
    if response is None:
        with fits.open(healpixels[pixel]) as hdulist:
            response = jsonify(hdulist['B_FLUX'].header['NAXIS2'])
            response.headers.add('Access-Control-Allow-Origin', '*')
        cache.set(str(pixel), response, timeout=3600)
    return response


@app.route('/<int:pixel>/<int:spectrum>')
def flux(pixel, spectrum):
    """Return an individual spectrum.
    """
    cacheID = "{0:d}-{1:d}".format(pixel, spectrum)
    response = cache.get(cacheID)
    if response is None:
        data = get_spectrum(pixel, spectrum)
        response = jsonify(data)
        response.headers.add('Access-Control-Allow-Origin', '*')
        cache.set(cacheID, response, timeout=3600)
    return response


def get_spectrum(pixel, spectrum):
    """Read spectrum data from files.

    Parameters
    ----------
    pixel : :class:`int`
        The healpixel.
    spectrum : :class:`int`
        The spectrum number.

    Returns
    -------
    :class:`dict`
        A dictionary suitable for conversion to JSON.
    """
    data = dict(pixel=pixel, spectrum=spectrum)
    with fits.open(healpixels[pixel]) as hdulist:
        # hdulist.info()
        for spectrograph in 'brz':
            data[spectrograph] = dict()
            data[spectrograph]['wavelength'] = hdulist['{0}_WAVELENGTH'.format(spectrograph.upper())].data.tolist()
            data[spectrograph]['flux'] = hdulist['{0}_FLUX'.format(spectrograph.upper())].data[data['spectrum']].tolist()
            data[spectrograph]['ivar'] = hdulist['{0}_IVAR'.format(spectrograph.upper())].data[data['spectrum']].tolist()
        data['desi_target'] = str(hdulist['FIBERMAP'].data[data['spectrum']]['DESI_TARGET'])
        data['bgs_target'] = str(hdulist['FIBERMAP'].data[data['spectrum']]['BGS_TARGET'])
        data['mws_target'] = str(hdulist['FIBERMAP'].data[data['spectrum']]['MWS_TARGET'])
        data['night'] = int(hdulist['FIBERMAP'].data[data['spectrum']]['NIGHT'])
        data['expid'] = int(hdulist['FIBERMAP'].data[data['spectrum']]['EXPID'])
        data['tileid'] = int(hdulist['FIBERMAP'].data[data['spectrum']]['TILEID'])
    with fits.open(os.path.join(os.path.dirname(healpixels[pixel]),
                                os.path.basename(healpixels[pixel]).replace('spectra', 'zbest'))) as hdulist:
        fibermap = hdulist['FIBERMAP'].data[data['spectrum']]
        targetid = int(fibermap['TARGETID'])
        data['targetid'] = str(targetid)
        data['ra'] = float(fibermap['RA_TARGET'])
        data['dec'] = float(fibermap['DEC_TARGET'])
        zcatalog = hdulist['ZBEST'].data
        w = zcatalog['TARGETID'] == targetid
        if w.any():
            row = zcatalog[w.nonzero()[0]]
            data['redshift'] = float(row['Z'])
            data['zerr'] = float(row['ZERR'])
            data['zwarn'] = float(row['ZWARN'])
            data['type'] = str(row['SPECTYPE'][0])
            data['subtype'] = str(zcatalog[w.nonzero()[0]]['SUBTYPE'][0])
            # data['desi_target'] = str(row['DESI_TARGET'])
            # data['bgs_target'] = str(row['BGS_TARGET'])
            # data['mws_target'] = str(row['MWS_TARGET'])
        # else:
        #     print('Target {0:d} not found in ZCATALOG.'.format(spectrum['targetid']))
    return data


if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=56789)
