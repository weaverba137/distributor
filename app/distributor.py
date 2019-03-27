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


class CompressedFloatJSONEncoder(json.JSONEncoder):
    _sig = 4  # Number of significant figures to keep.
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            if obj.dtype.kind == 'f':
                return round_to_sig_fig(obj, self._sig)
            else:
                return obj.tolist()
        return super(CompressedFloatJSONEncoder, self).default(obj)
    def encode(self, obj):
        e = super(CompressedFloatJSONEncoder, self).encode(obj)
        return e.replace('"@', '').replace('@"', '')


app = Flask(__name__)
app.json_encoder = CompressedFloatJSONEncoder
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


def get_spectrum(pixel, spectrum, tolist=False):
    """Read spectrum data from files.

    Parameters
    ----------
    pixel : :class:`int`
        The healpixel.
    spectrum : :class:`int`
        The spectrum number.
    tolist : :class:`bool`, optional
        If ``True``, coerce flux and wavelength values to Python :class:`list`.

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
            w = hdulist['{0}_WAVELENGTH'.format(spectrograph.upper())].data
            data[spectrograph]['w0'] = float(w[0])
            data[spectrograph]['dw'] = float(np.diff(w).mean())
            flux = hdulist['{0}_FLUX'.format(spectrograph.upper())].data[data['spectrum']]
            ivar = hdulist['{0}_IVAR'.format(spectrograph.upper())].data[data['spectrum']]
            if tolist:
                flux = flux.tolist()
                ivar = ivar.tolist()
            data[spectrograph]['flux'] = flux
            data[spectrograph]['ivar'] = ivar
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
        data['ra'] = float(fibermap['TARGET_RA'])
        data['dec'] = float(fibermap['TARGET_DEC'])
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


# The following constant was computed in maxima 5.35.1 using 64 bigfloat digits of precision
__logBase10of2 = 3.010299956639811952137388947244930267681898814621085413104274611e-1


def round_to_sig_fig(x, sig, format=True):
    """Round and format the values in `x` to the number of significant figures in `sig`.

    Parameters
    ----------
    x : :class:`numpy.ndarray`
        Values to round.
    sig : :class:`int`
        Number of significant figures.
    format : :class:`bool`, optional
        If ``True``, format the rounded values to string.

    Returns
    -------
    :class:`numpy.ndarray` or :class:`list`
        A rounded array of the same type as `x` or a list of strings.

    Notes
    -----
    Adapted from code found in `SELPythonLibs <https://github.com/odysseus9672/SELPythonLibs/blob/master/SigFigRounding.py>`_ .
    """
    if not np.isreal(x).all():
        raise TypeError("round_to_sig_fig: x must be real.")

    xsgn = np.sign(x)
    absx = xsgn * x
    mantissa, binaryExponent = np.frexp( absx )
    decimalExponent = __logBase10of2 * binaryExponent
    omag = np.floor(decimalExponent)
    mantissa *= 10.0**(decimalExponent - omag)
    w = mantissa < 1.0
    if w.any():
        mantissa[w] *= 10.0
        omag[w] -= 1.0
    mantissa = xsgn * np.around(mantissa, decimals=sig-1)
    fmt = "@{:." + str(sig-1) + "g}@"
    if format:
        return [fmt.format(m*10.0**e) for m, e in zip(mantissa.tolist(), omag.tolist())]
    return mantissa * 10.0**omag




if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=56789)
