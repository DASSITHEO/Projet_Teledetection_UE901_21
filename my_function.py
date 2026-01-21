import numpy as np
from osgeo import gdal, osr, ogr
import pandas as pd

def calculate_ari(band3_path, band5_path, output_path):
    """
    Calcule le Normalized Anthocyanin Reflectance Index (ari) à partir de deux bandes
    et sauvegarde le résultat en TIFF géoréférencé avec float32, projection EPSG:32630,
    résolution 10m et NoData = -9999.
    
    Parameters:
    -----------
    band3_path : str
        Chemin vers la bande 03 (B03).
    band5_path : str
        Chemin vers la bande 05 (B05).
    output_path : str
        Chemin pour sauvegarder le TIFF de sortie.
    """
    
    # Charger les bandes 
    band3_ds = gdal.Open(band3_path)
    band5_ds = gdal.Open(band5_path)
    
    band3 = band3_ds.GetRasterBand(1).ReadAsArray().astype(float)
    band5 = band5_ds.GetRasterBand(1).ReadAsArray().astype(float)
    

    # Calculer le ari 
    ari = (1/band3 - 1/band5) / (1/band3 + 1/band5)
    
    # . Remplacer les valeurs invalides par NoData 
    ari[np.isnan(ari)] = -9999
    ari[np.isinf(ari)] = -9999
    
    #  Créer le nouveau TIFF 
    driver = gdal.GetDriverByName('GTiff')
    rows, cols = ari.shape
    out_ds = driver.Create(output_path, cols, rows, 1, gdal.GDT_Float32)
    
    #  Définir la projection et la transformation
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(32630)
    out_ds.SetProjection(srs.ExportToWkt())
    
    # Récupérer GeoTransform de la bande 3 (ou adapter pour 10m si nécessaire)
    geo_transform = band3_ds.GetGeoTransform()
    out_ds.SetGeoTransform(geo_transform)
    
    #  Écrire les données et définir NoData 
    out_band = out_ds.GetRasterBand(1)
    out_band.WriteArray(ari)
    out_band.SetNoDataValue(-9999)
    
    #  Sauvegarder et fermer 
    out_band.FlushCache()
    out_ds.FlushCache()

    out_band = None
    out_ds = None
    band3_ds = None
    band5_ds = None
    print("Fichier ARI_serie_temp.tif de l'indice ARI a été créé avec succès !")
    


def print_raster_info(raster_path):

    ds = gdal.Open(raster_path)
    if ds is None:
        print("Impossible d’ouvrir le raster")
        return

    # Résolution spatiale
    gt = ds.GetGeoTransform()
    pixel_size_x = gt[1]
    pixel_size_y = abs(gt[5])

    # Type d’encodage
    band = ds.GetRasterBand(1)
    data_type = gdal.GetDataTypeName(band.DataType)

    # NoData
    nodata = band.GetNoDataValue()

    # Projection (optionnel mais utile pour vérifier)
    proj = ds.GetProjection()

    print("---- Informations du raster ----")
    print("Résolution spatiale X :", pixel_size_x, "m")
    print("Résolution spatiale Y :", pixel_size_y, "m")
    print("Type d’encodage :", data_type)
    print("Valeur NoData :", nodata)
    print("Projection :", proj)

    ds = None



def report_from_dict_to_df(dict_report):


    # convert report into dataframe
    report_df = pd.DataFrame.from_dict(dict_report)

    # drop unnecessary rows and columns
    try :
        report_df = report_df.drop(['accuracy', 'macro avg', 'weighted avg'], axis=1)
    except KeyError:
        print(dict_report)
        report_df = report_df.drop(['micro avg', 'macro avg', 'weighted avg'], axis=1)

    report_df = report_df.drop(['support'], axis=0)

    return report_df




def rasterize_shapefile(shapefile, ref_raster, out_raster, field_name):
    """
    Rasterize a vector shapefile to match a reference raster using GDAL.
    """
    # Ouvrir raster de référence
    ref_ds = gdal.Open(ref_raster)
    x_res = ref_ds.RasterXSize
    y_res = ref_ds.RasterYSize
    geotransform = ref_ds.GetGeoTransform()
    proj = ref_ds.GetProjection()

    # Créer raster de sortie
    driver = gdal.GetDriverByName("GTiff")
    out_ds = driver.Create(out_raster, x_res, y_res, 1, gdal.GDT_Byte)
    out_ds.SetGeoTransform(geotransform)
    out_ds.SetProjection(proj)

    # Ouvrir shapefile
    shp_ds = ogr.Open(shapefile)
    shp_layer = shp_ds.GetLayer()

    # Rasteriser
    gdal.RasterizeLayer(out_ds, [1], shp_layer, options=[f"ATTRIBUTE={field_name}"])

    # Fermer fichiers
    out_ds = None
    shp_ds = None

    print(f"Raster créé : {out_raster}")

    import numpy as np



def merge_and_save_multiband(
    images,
    output_path,
    reference_raster_path,
    nodata=None,
    dtype=gdal.GDT_Float32
):
    """
    Fusionne plusieurs images 2D en une image multibande et l'enregistre avec GDAL.

    Parameters
    ----------
    images : list of ndarray (H, W)
        Images à fusionner (mêmes dimensions)
    output_path : str
        Chemin du GeoTIFF de sortie
    reference_raster_path : str
        Raster de référence pour la géoréférence
    nodata : float, optional
        Valeur nodata
    dtype : gdal data type
        Type des données en sortie
    """

    # Vérification des dimensions
    shapes = [img.shape for img in images]
    if len(set(shapes)) != 1:
        raise ValueError("Toutes les images doivent avoir les mêmes dimensions")

    # Empilement → image multibande
    multiband_image = np.stack(images, axis=-1)

    rows, cols, bands = multiband_image.shape

    # Lecture de la géoréférence
    ref_ds = gdal.Open(reference_raster_path)
    geotransform = ref_ds.GetGeoTransform()
    projection = ref_ds.GetProjection()
    ref_ds = None

    # Création du raster de sortie
    driver = gdal.GetDriverByName("GTiff")
    out_ds = driver.Create(
        output_path,
        cols,
        rows,
        bands,
        dtype
    )

    out_ds.SetGeoTransform(geotransform)
    out_ds.SetProjection(projection)

    # Écriture des bandes
    for b in range(bands):
        band = out_ds.GetRasterBand(b + 1)
        band.WriteArray(multiband_image[:, :, b])

        if nodata is not None:
            band.SetNoDataValue(nodata)

        band.FlushCache()

    out_ds = None  # fermeture du fichier


def write_classified_image_nodata0(
    out_filename,
    classified_image,
    ds_ref,
    gdal_dtype
):
    """
    Enregistre une image classifiée en utilisant une image de référence,
    en forçant NoData = 0.

    Parameters
    ----------
    out_filename : str
        Chemin du raster de sortie
    classified_image : np.ndarray
        Image classifiée (2D)
    ds_ref : gdal.Dataset
        Dataset GDAL de référence (CRS, transform, taille)
    gdal_dtype : gdal.DataType
        Type GDAL (ex : gdal.GDT_Byte, gdal.GDT_Int16)
    """

    # Sécurité : forcer les NoData à 0 dans l'array
    classified_image = classified_image.copy()
    classified_image[np.isnan(classified_image)] = 0

    # Création du raster
    driver = gdal.GetDriverByName("GTiff")
    out_ds = driver.Create(
        out_filename,
        ds_ref.RasterXSize,
        ds_ref.RasterYSize,
        1,
        gdal_dtype,
        options=["COMPRESS=LZW"]
    )

    # Copier la géoréférence
    out_ds.SetGeoTransform(ds_ref.GetGeoTransform())
    out_ds.SetProjection(ds_ref.GetProjection())

    # Écriture de la bande
    band = out_ds.GetRasterBand(1)
    band.WriteArray(classified_image)

    # FORÇAGE EXPLICITE DU NoData
    band.SetNoDataValue(0)

    band.FlushCache()
    out_ds = None

    print(f" Image enregistrée avec un type int8 et NoData = 0 : {out_filename}")
