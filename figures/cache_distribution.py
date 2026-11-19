#!/usr/bin/env python3
"""Cache the Araucaria distribution from the CONAF Catastro so F1 is reproducible without the
3.1 GB source. Two definitions are extracted: Araucaria as the DOMINANT forest type (tipo_fores,
the definition behind the published 201,264 ha figure) and Araucaria present as ANY component
species, which is 65% larger. Sampling coverage is computed under both."""
import shapefile, numpy as np, pickle, pyproj, time
FILES=["gis/conaf/cbn_09reg_2024.shp","gis/conaf/cbn_08reg_2024.shp"]
tf=pyproj.Transformer.from_crs("EPSG:32718","EPSG:4326",always_xy=True)
rings_dom, rings_extra, cent_dom, cent_any, land = [],[],[],[],[]
n_dom=n_any=0; ha_dom=ha_any=0.0
t0=time.time()
for f in FILES:
    r=shapefile.Reader(f); flds=[x[0] for x in r.fields[1:]]
    i_tf=flds.index("tipo_fores"); i_ha=flds.index("superf_ha")
    esp=[flds.index(c) for c in flds if c.startswith("especi") and c.endswith("_ci")]
    for sr in r.iterShapeRecords():
        rec=sr.record; bb=sr.shape.bbox
        land.append(bb)
        a_tf="araucaria" in str(rec[i_tf]).lower()
        a_any=a_tf or any("araucaria" in str(rec[j]).lower() for j in esp)
        if not a_any: continue
        h=float(rec[i_ha] or 0)
        cx,cy=tf.transform((bb[0]+bb[2])/2,(bb[1]+bb[3])/2)
        n_any+=1; ha_any+=h; cent_any.append((cx,cy,h))
        if a_tf: n_dom+=1; ha_dom+=h; cent_dom.append((cx,cy,h))
        pts=np.asarray(sr.shape.points,dtype=float); parts=list(sr.shape.parts)+[len(pts)]
        for a,b in zip(parts[:-1],parts[1:]):
            seg=pts[a:b]
            if len(seg)<3: continue
            lon,lat=tf.transform(seg[:,0],seg[:,1])
            ring=np.column_stack([lon,lat]).astype(np.float32)
            (rings_dom if a_tf else rings_extra).append(ring)
L=np.array(land)
lo1,la1=tf.transform(L[:,0],L[:,1]); lo2,la2=tf.transform(L[:,2],L[:,3])
gx=np.arange(-74.4,-70.6,0.02); gy=np.arange(-40.1,-36.9,0.02)
grid=np.zeros((len(gy),len(gx)),dtype=bool)
ix1=np.clip(np.searchsorted(gx,lo1)-1,0,len(gx)-1); ix2=np.clip(np.searchsorted(gx,lo2),0,len(gx)-1)
iy1=np.clip(np.searchsorted(gy,la1)-1,0,len(gy)-1); iy2=np.clip(np.searchsorted(gy,la2),0,len(gy)-1)
for a,b,c,d in zip(iy1,iy2,ix1,ix2): grid[a:b+1,c:d+1]=True

S=np.load("figdata/stand_coords.npz",allow_pickle=True); SC=S["coords"]
def km(lon1,lat1,lon2,lat2):
    R=6371.0; p1,p2=np.radians(lat1),np.radians(lat2)
    h=np.sin((p2-p1)/2)**2+np.cos(p1)*np.cos(p2)*np.sin(np.radians(lon2-lon1)/2)**2
    return 2*R*np.arcsin(np.sqrt(h))
cov={}
for lab,C in (("dominant",np.array(cent_dom)),("any component",np.array(cent_any))):
    dmin=np.min([km(C[:,0],C[:,1],s[0],s[1]) for s in SC],axis=0)
    tot=C[:,2].sum()
    cov[lab]=dict(pct50=100*C[dmin<=50,2].sum()/tot, maxd=float(dmin.max()),
                  pct25=100*C[dmin<=25,2].sum()/tot, ha=float(tot), n=len(C))
    print(f"{lab:>15}: {len(C):,} polygons  {tot:,.0f} ha | within 50 km {cov[lab]['pct50']:.1f}% | "
          f"within 25 km {cov[lab]['pct25']:.1f}% | farthest {cov[lab]['maxd']:.1f} km")
pickle.dump(dict(rings_dom=rings_dom,rings_extra=rings_extra,grid=grid,gx=gx,gy=gy,
                 n_dom=n_dom,ha_dom=ha_dom,n_any=n_any,ha_any=ha_any,cov=cov),
            open("figdata/distribution_cache.pkl","wb"))
print(f"\ncached in {time.time()-t0:.0f}s -> figdata/distribution_cache.pkl")
