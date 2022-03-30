"""
Thaddaeus mentioned only a small fraction (~10%) of the Kerr+21 Table 1 stars
had 2 minute data.  Let's verify: is it really that small?
"""
import os
from os.path import join
import pandas as pd, numpy as np, matplotlib.pyplot as plt, lightkurve as lk
from sunnyhills.paths import DATADIR
from astropy import units as u
from astropy.coordinates import SkyCoord
from astrobase.services.identifiers import gaiadr2_to_tic
from glob import glob

df = pd.read_csv(join(DATADIR, 'current', 'current_key.csv'))

cachedir = './temp_cache'
if not os.path.exists(cachedir):
    os.mkdir(cachedir)

# "cadences" is a string like "120,1800" listing out all available
# cadences, in units of seconds.

df['GDR2_ID'] = df['GDR2_ID'].astype(str)

for ix, r in df.iterrows():

    tic_id = gaiadr2_to_tic(r['GDR2_ID'])

    print(f"{ix}/{len(df)}...")
    cachepath = join(cachedir, tic_id+'.csv')
    if os.path.exists(cachepath):
        print(f"Found {cachepath}, continue.")
        continue

    search_str = 'TIC '+tic_id
    lcset = lk.search_lightcurve(search_str)

    if len(lcset) > 0:
        cadences = ','.join(np.unique(lcset.exptime.value).astype(int).astype(str))
    else:
        cadences = 'No TESS data found'

    r['cadences'] = cadences
    r['tic_id_corrected'] = tic_id
    pd.DataFrame(r).T.to_csv(cachepath, index=False, sep="|")
    print(f"Wrote {cachepath}")

cachepaths = glob(join(cachedir, '*.csv'))
_df = pd.concat((pd.read_csv(f, sep="|") for f in cachepaths))

c = SkyCoord(ra=np.array(_df.GDR2_RA)*u.deg, dec=np.array(_df.GDR2_DEC)*u.deg,
             frame='icrs')
_df['ecl_lat'] = c.barycentrictrueecliptic.lat.value
_df['ecl_lon'] = c.barycentrictrueecliptic.lon.value
_df['gal_lat'] = c.galactic.b.value
_df['gal_lon'] = c.galactic.l.value

outpath = join(DATADIR, 'current', 'current_key_with_cadences.csv')
_df.to_csv(outpath, sep="|", index=False)

has_2min = _df.cadences.astype(str).str.contains('120')
has_any = _df.cadences != 'No TESS data found'

print(f'N={len(_df)} objects')
print(f'N={len(_df[has_2min])} has 2 minute data')
print(f'N={len(_df[has_any])} has any data')
print(f'N={len(_df[~has_any])} has no data')

# quick-look plot to see where there is no data
xytuples = [
    ('GDR2_RA','GDR2_DEC', 'linear', 'linear'),
    ('ecl_lon','ecl_lat', 'linear', 'linear'),
    ('gal_lon','gal_lat', 'linear', 'linear'),
    ('Distance', 'g_mag', 'linear', 'linear'),
    ('bp-rp', 'g_mag', 'linear', 'linear'),
]

for xy in xytuples:

    xkey, ykey = xy[0], xy[1]
    xscale, yscale = xy[2], xy[3]

    plt.close('all')
    fig, ax = plt.subplots()
    ax.scatter(_df[xkey], _df[ykey],
               s=2, c='darkgray', zorder=1,
               label=f'all N={len(_df)} objects')
    ax.scatter(_df[has_any][xkey], _df[has_any][ykey],
               s=2, c='black', zorder=2,
               label=f'N={len(_df[has_any])} objects with data')
    ax.scatter(_df[has_2min][xkey], _df[has_2min][ykey],
               s=1, c='red', zorder=3,
               label=f'N={len(_df[has_any])} objects with 2min data')
    ax.update({
        'xlabel':xkey,
        'ylabel':ykey,
        'xscale':xscale,
        'yscale':yscale
    })
    fig.savefig(f'figures/verify_2min_{xkey}_{ykey}.png', dpi=400, bbox_inches='tight')

# quick analysis based on the results of those plots...
sel = (_df.ecl_lat > -5) & (_df.ecl_lat < 20) & (_df.ecl_lon > 260) & (_df.ecl_lon < 280)
sdf = _df[sel]

# nothing here should exist, based on my sky-maps, or Roland Vanderspek's (https://tess.mit.edu/tess-year-4-observations/, ditto years 1-4)
print(sdf[sdf.cadences.astype(str).str.contains('120')])

import IPython; IPython.embed()
