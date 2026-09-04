"""Generate metadata.csv for CONUS_Forest_History tifs for geeup upload."""
import os
import csv
import calendar
from osgeo import gdal

tif_dir = '/Users/Sayantan.Majumdar@dri.edu/VSCode/pyCropWat_Data/Precip'
output_csv = os.path.join(tif_dir, 'metadata.csv')

tif_files = sorted([f for f in os.listdir(tif_dir) if f.endswith('.tif')])

dtype_map = {
    gdal.GDT_Byte: "Byte", gdal.GDT_UInt16: "UInt16", gdal.GDT_Int16: "Int16",
    gdal.GDT_UInt32: "UInt32", gdal.GDT_Int32: "Int32",
    gdal.GDT_Float32: "Float32", gdal.GDT_Float64: "Float64"
}

rows = []

for tif in tif_files:
    filepath = os.path.join(tif_dir, tif)
    system_index = tif.replace('.tif', '')
    _, year, month = tif[:tif.find('.')].split('_')
    year = int(year)
    month = int(month)
    ds = gdal.Open(filepath)
    if ds:
        xsize, ysize = ds.RasterXSize, ds.RasterYSize
        band = ds.GetRasterBand(1)
        data_type = dtype_map.get(band.DataType, "Byte")
        ci = band.GetColorInterpretation()
        color_interp = "Palette" if ci == gdal.GCI_PaletteIndex else "Gray"
        ds = None
    else:
        print(f"Warning: Could not open {filepath}. Skipping.")
        continue
    num_days_in_month = calendar.monthrange(year, month)[1]
    
    rows.append({
        'system:index': system_index,
        'year': year,
        'month': month,
        'xsize': xsize,
        'ysize': ysize,
        'num_bands': 1,
        'data_type': data_type,
        'color_interpretation': color_interp,
        'inferred_kind': 'continuous',
        'system:time_start': f'{year}-{month:02d}-01',
        'system:time_end': f'{year}-{month:02d}-{num_days_in_month:02d}'
    })

with open(output_csv, 'w', newline='') as f:
    fieldnames = ['system:index', 'year', 'month', 'xsize', 'ysize', 'num_bands', 'data_type', 
                  'color_interpretation', 'inferred_kind', 'system:time_start', 'system:time_end']
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

print(f"Created metadata.csv with {len(rows)} entries")
print(f"Output: {output_csv}")
