import json, urllib.request, numpy as np, rasterio
from rasterio.windows import from_bounds
from rasterio.warp import transform_bounds
BBOX=(39.05,21.25,39.40,21.85)  # Jeddah urban + lower wadis
url=("https://stac.eodc.eu/api/v1/search?collections=GFM&bbox=%s,%s,%s,%s"
     "&datetime=2022-11-15T00:00:00Z/2022-12-05T00:00:00Z&limit=50")%BBOX
items=json.load(urllib.request.urlopen(url,timeout=60))['features']
for it in sorted(items,key=lambda i:i['id']):
    res={}
    for k in ['ensemble_flood_extent','exclusion_mask']:
        with rasterio.open(it['assets'][k]['href']) as ds:
            b=transform_bounds('EPSG:4326',ds.crs,*BBOX)
            w=from_bounds(*b,ds.transform).intersection(rasterio.windows.Window(0,0,ds.width,ds.height))
            a=ds.read(1,window=w)
            res[k]=(a, ds.nodata)
    fe,nd=res['ensemble_flood_extent']; ex,_=res['exclusion_mask']
    vals,cnt=np.unique(fe,return_counts=True)
    print(it['id'][15:30], it['id'][-12:], 'win',fe.shape,'flood vals',dict(zip(vals.tolist(),cnt.tolist())),'excluded px',int((ex==1).sum()))
