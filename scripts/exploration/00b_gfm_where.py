import json, urllib.request, numpy as np, rasterio
from rasterio.windows import from_bounds
from rasterio.warp import transform_bounds
from pyproj import Transformer
BBOX=(39.05,21.25,39.40,21.85)
url=("https://stac.eodc.eu/api/v1/search?collections=GFM&bbox=%s,%s,%s,%s"
     "&datetime=2022-11-29T15:29:00Z/2022-11-29T15:30:00Z&limit=50")%BBOX
items=json.load(urllib.request.urlopen(url,timeout=60))['features']
pts=[]
for it in items:
    with rasterio.open(it['assets']['ensemble_flood_extent']['href']) as ds:
        b=transform_bounds('EPSG:4326',ds.crs,*BBOX)
        w=from_bounds(*b,ds.transform).intersection(rasterio.windows.Window(0,0,ds.width,ds.height))
        a=ds.read(1,window=w); T=ds.window_transform(w)
        r,c=np.nonzero(a==1)
        xs,ys=rasterio.transform.xy(T,r,c)
        lon,lat=Transformer.from_crs(ds.crs,'EPSG:4326',always_xy=True).transform(xs,ys)
        pts.append(np.c_[lon,lat])
        tot=(a!=255).sum()
        with rasterio.open(it['assets']['exclusion_mask']['href']) as dx:
            ex=dx.read(1,window=w)
        print(it['id'][-12:],'observed px',int(tot),'excluded share of observed %.1f%%'%(100*((ex==1)&(a!=255)).sum()/max(tot,1)))
P=np.vstack(pts); print('total flood px',len(P),'~area km2 %.2f'%(len(P)*400/1e6))
H,xe,ye=np.histogram2d(P[:,0],P[:,1],bins=[np.arange(39.05,39.45,0.05),np.arange(21.25,21.90,0.05)])
print('lon\lat '+' '.join('%5.2f'%y for y in ye[:-1]))
for i,x in enumerate(xe[:-1]): print('%6.2f  '%x+' '.join('%5d'%v for v in H[i]))
