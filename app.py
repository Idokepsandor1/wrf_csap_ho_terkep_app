import streamlit as st
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shapereader
from wrf import getvar, to_np, latlon_coords, get_cartopy
import geopandas as gpd
from cartopy.mpl.patch import geos_to_path
import matplotlib.patches as mpatches
from matplotlib import font_manager
import imageio.v2 as imageio
import requests
import tempfile
import os
from datetime import datetime, timedelta

# --- KONFIGURÁCIÓ ---
st.set_page_config(page_title="WRF Térkép Generátor", layout="wide")

# Fájlok (Mivel a GitHubon vannak, a kód mellett lesznek)
DESIGN_PNG_PATH = "logo_cim_idokep.png"
FONT_TITLE_PATH = "Oswald-Bold.ttf"
FONT_SUBTITLE_PATH = "Inter_18pt-SemiBold.ttf"
ADMIN_GEOJSON_PATH = "hungary.geojson"
PNG_RELIEF_PATH = "height_map_ceu.png"

# Modellek
MODEL_URLS = {
    "ICON_PEST": "https://pest.idokep.hu/run/",
    "ICON_BUDA": "https://buda.idokep.hu/run/",
    "GFS_JUMBO": "https://jumbo.idokep.hu/run/",
    "GFS_BALLON": "https://ballon.idokep.hu/run/"
}

# Színek
PRECIP_HEX = ['#0233a3', '#0b55ff', '#00a0ff', '#0b55ff', '#00a0ff', '#089908', '#2be007', '#ffff00', '#ffc800', '#ff7700', '#ff7700', '#ff0000', '#c80000', '#7d0101', '#7d0101', '#fc0390', '#ff00ff', '#fc0390', '#f75cb4', '#f75cb4']
SNOW_HEX = ['#bccef7', '#7ca3fc', '#5285fa', '#2564f7', '#0048f0', '#0332a1', '#1d009e', '#5b00f7', '#9603ff', '#cf03fc', '#fc03ca', '#fc0388', '#fc0345', '#fc030b', '#fc5e03', '#d44a02', '#a83602', '#7f2602', '#5a1b02', '#3c1101']

# Cache-elt erőforrások (hogy gyorsabb legyen)
@st.cache_data
def get_geometries():
    budapest = None
    hungary = None
    
    if os.path.exists(ADMIN_GEOJSON_PATH):
        gdf = gpd.read_file(ADMIN_GEOJSON_PATH).to_crs("EPSG:4326")
        budapest = gdf[gdf['name'] == 'Budapest'].geometry.iloc[0]
        
    try:
        reader = shapereader.Reader(shapereader.natural_earth('10m', 'cultural', 'admin_0_countries'))
        hungary = next((c.geometry for c in reader.records() if c.attributes['NAME'] == 'Hungary'), None)
    except:
        pass
        
    return budapest, hungary

# --- FÜGGVÉNYEK ---
def get_urls(start, end, base_url):
    urls = []
    curr = start
    while curr <= end:
        url = f"{base_url}wrfout_d02_{curr.strftime('%Y-%m-%d_%H_%M_%S')}"
        urls.append(url)
        curr += timedelta(minutes=15)
    return urls

def download_nc(url):
    try:
        r = requests.head(url, timeout=2)
        if r.status_code != 200: return None
        r = requests.get(url, stream=True, timeout=30)
        tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".nc")
        with tmp as f:
            for chunk in r.iter_content(8192): f.write(chunk)
        return tmp.name
    except:
        return None

def get_hungarian_day_name(dt):
    return ["hétfő", "kedd", "szerda", "csütörtök", "péntek", "szombat", "vasárnap"][dt.weekday()]

def get_hungarian_month_name(m):
    return ["január", "február", "március", "április", "május", "június", "július", "augusztus", "szeptember", "október", "november", "december"][m-1]

def get_time_of_day_suffix(h):
    if 0<=h<1: return "éjfélig"
    elif 1<=h<6: return "hajnalig"
    elif 6<=h<12: return "reggelig"
    elif 12<=h<16: return "délutánig"
    else: return "estig"

def add_media_design(fig, start_dt, end_dt, p_type):
    if not os.path.exists(DESIGN_PNG_PATH): return
    
    try:
        ax = fig.add_axes([0, 0, 1, 1], zorder=50)
        ax.axis('off')
        img = plt.imread(DESIGN_PNG_PATH)
        ax.imshow(img, aspect='auto', extent=[0, 1, 0, 1], origin='upper')

        day = get_hungarian_day_name(end_dt).upper()
        suffix = get_time_of_day_suffix(end_dt.hour).upper()
        type_str = "HÓRÉTEG" if p_type == 'Hó' else "CSAPADÉK"
        full_title = f"{day} {suffix} VÁRHATÓ {type_str}"
        
        now = datetime.now()
        now = now.replace(minute=(now.minute//5)*5, second=0)
        month = get_hungarian_month_name(now.month)
        subtitle = f"Kiadva: {now.year}. {month} {now.day}. {now.strftime('%H:%M')}"

        if os.path.exists(FONT_TITLE_PATH):
            fp_t = font_manager.FontProperties(fname=FONT_TITLE_PATH, size=28)
            fp_s = font_manager.FontProperties(fname=FONT_SUBTITLE_PATH, size=24)
            ax.text(0.045, 0.922, full_title, transform=ax.transAxes, color='#FFFFFF', fontproperties=fp_t, ha='left', va='center')
            ax.text(0.045, 0.865, subtitle, transform=ax.transAxes, color='#1E8FDB', fontproperties=fp_s, ha='left', va='center')
        else:
            ax.text(0.045, 0.92, full_title, transform=ax.transAxes, color='#FFFFFF', size=22)
    except: pass

# --- FŐOLDAL ---
st.title("🌦️ WRF Térkép Generátor")

col1, col2 = st.columns(2)
with col1:
    s_date = st.date_input("Kezdő dátum", datetime.now())
    s_hour = st.selectbox("Kezdő óra", [f"{h:02d}" for h in range(24)], index=6)
    p_type = st.selectbox("Típus", ["Csapadék", "Hó"])

with col2:
    e_date = st.date_input("Vég dátum", datetime.now() + timedelta(days=1))
    e_hour = st.selectbox("Vég óra", [f"{h:02d}" for h in range(24)], index=6)
    models = st.multiselect("Modellek", list(MODEL_URLS.keys()), default=list(MODEL_URLS.keys()))

show_contours = st.checkbox("Fekete izovonalak", value=False)

if st.button("Térképek Generálása", type="primary"):
    if not models:
        st.error("Válassz legalább egy modellt!")
        st.stop()

    s_dt = datetime.combine(s_date, datetime.strptime(s_hour, "%H").time())
    e_dt = datetime.combine(e_date, datetime.strptime(e_hour, "%H").time())
    
    progress_bar = st.progress(0)
    status_text = st.empty()
    
    # Színskálák
    if p_type == 'Hó':
        cmap = mcolors.ListedColormap(SNOW_HEX)
        levels = [0.1, 0.5, 1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 18, 20, 25, 30, 50, 100, 125, 150]
        unit = "cm"
    else:
        cmap = mcolors.ListedColormap(PRECIP_HEX)
        levels = [1, 2, 3, 5, 8, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 80, 100, 120, 150, 200]
        unit = "mm"
    norm = mcolors.BoundaryNorm(levels, cmap.N, clip=True)

    bp_geom, hu_geom = get_geometries()
    
    ensemble_data = []
    valid_models = []

    for idx, m_name in enumerate(models):
        status_text.text(f"{m_name} letöltése...")
        progress_bar.progress((idx) / len(models))
        
        urls = get_urls(s_dt, e_dt, MODEL_URLS[m_name])
        if len(urls) < 2: continue
        
        # Letöltés
        f_start = download_nc(urls[0])
        f_end = download_nc(urls[-1])
        
        if f_start and f_end:
            try:
                with netCDF4.Dataset(f_start) as nc_s, netCDF4.Dataset(f_end) as nc_e:
                    # Vetület és koordináták kinyerése (csak az elsőnél)
                    if len(ensemble_data) == 0:
                        wrf_var = getvar(nc_s, "T2", timeidx=0)
                        cart_proj = get_cartopy(wrf_var)
                        lats, lons = latlon_coords(wrf_var)
                        lats, lons = to_np(lats), to_np(lons)

                    if p_type == 'Hó':
                        val = to_np(getvar(nc_e, "SNOWH", timeidx=0)) * 100
                    else:
                        r_s = to_np(getvar(nc_s, "RAINNC", timeidx=0) + getvar(nc_s, "RAINC", timeidx=0))
                        r_e = to_np(getvar(nc_e, "RAINNC", timeidx=0) + getvar(nc_e, "RAINC", timeidx=0))
                        val = r_e - r_s
                    
                    val[val < 0] = 0
                    ensemble_data.append(val)
                    valid_models.append(m_name)
            except Exception as e:
                st.error(f"Hiba {m_name}: {e}")
            finally:
                if os.path.exists(f_start): os.remove(f_start)
                if os.path.exists(f_end): os.remove(f_end)
        else:
            st.warning(f"{m_name}: Adat nem elérhető.")

    progress_bar.progress(1.0)
    status_text.empty()
    
    if not ensemble_data:
        st.error("Nem sikerült adatot letölteni.")
        st.stop()

    # Térképek kirajzolása
    st.subheader("Eredmények")
    
    # Tagok és Átlag feldolgozása
    to_plot = list(zip(valid_models, ensemble_data))
    if len(ensemble_data) > 1:
        to_plot.append(("ENSEMBLE ÁTLAG", np.mean(np.stack(ensemble_data), axis=0)))

    for title, data in to_plot:
        # Térkép készítése
        w, h, dpi = 12, 10, 150
        if os.path.exists(DESIGN_PNG_PATH):
            dimg = plt.imread(DESIGN_PNG_PATH)
            w = 16; h = w * (dimg.shape[0]/dimg.shape[1])
            
        fig = plt.figure(figsize=(w, h), dpi=dpi)
        ax = fig.add_axes([0, 0, 1, 1], projection=cart_proj)
        ax.set_extent([16, 23, 45.5, 48.8], crs=ccrs.PlateCarree())
        
        ax.add_feature(cfeature.LAND.with_scale('10m'), facecolor='#FFF9F9')
        
        # Domborzat (egyszerűsítve)
        if os.path.exists(PNG_RELIEF_PATH):
             try:
                 ir = imageio.imread(PNG_RELIEF_PATH)
                 ig = ir[:,:,:3].mean(axis=2) if ir.ndim==3 else ir
                 el = ig.astype(float)*10.0
                 ls = mcolors.LightSource(azdeg=315, altdeg=45)
                 rhs = ls.shade(el, cmap=plt.cm.Greys, vert_exag=200.0, blend_mode='soft')
                 rgba_hill = np.dstack((rhs[:,:,:3], np.where(el>120, 1.0, 0.0)))
                 ax.imshow(rgba_hill, origin='upper', extent=[5.0, 27.0, 40.4175, 51.4175], transform=ccrs.PlateCarree(), zorder=1)
             except: pass

        ax.add_feature(cfeature.LAKES.with_scale('10m'), facecolor='#6baed6', zorder=2)
        ax.add_feature(cfeature.RIVERS.with_scale('10m'), edgecolor='#6baed6', zorder=2)

        # Adat
        if np.any(data > 0) and hu_geom:
             clip = mpatches.PathPatch(geos_to_path(hu_geom)[0], transform=ax.transData)
             cf = ax.contourf(lons, lats, data, levels=levels, cmap=cmap, norm=norm, extend='max', transform=ccrs.PlateCarree(), zorder=4, alpha=0.6)
             cf2 = ax.contourf(lons, lats, data, levels=levels, cmap=cmap, norm=norm, extend='max', transform=ccrs.PlateCarree(), zorder=4, alpha=1.0)
             for col in cf2.collections: col.set_clip_path(clip)
        
        if show_contours:
             clvls = [1, 2, 5, 10, 20, 50]
             cs = ax.contour(lons, lats, data, levels=clvls, colors='black', linewidths=0.8, transform=ccrs.PlateCarree(), zorder=5)
             ax.clabel(cs, inline=True, fmt=f'%d {unit}', fontsize=9)

        ax.add_feature(cfeature.BORDERS.with_scale('10m'), linestyle=':', zorder=6)
        if bp_geom: ax.add_geometries([bp_geom], crs=ccrs.PlateCarree(), edgecolor='black', facecolor='none', zorder=7)
        if hu_geom: ax.add_geometries([hu_geom], crs=ccrs.PlateCarree(), edgecolor='black', facecolor='none', linewidth=2, zorder=8)

        add_media_design(fig, s_dt, e_dt, p_type)
        
        # Megjelenítés
        st.pyplot(fig)
        
        # Letöltés gomb
        fn = f"{p_type}_{title}.png"
        plt.savefig(fn, dpi=150, bbox_inches=None) # Fontos: nincs bbox_inches='tight'!
        with open(fn, "rb") as file:
            st.download_button(label=f"Letöltés: {title}", data=file, file_name=fn, mime="image/png")
        plt.close(fig)
