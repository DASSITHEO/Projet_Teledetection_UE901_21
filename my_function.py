import numpy as np
from osgeo import gdal, osr, ogr
import pandas as pd
import os

def calculate_ari(band3_path, band5_path, output_path):
    """
    Calcule l'Anthocyanin Reflectance Index (ARI) pour une série temporelle
    (multi-bandes) et sauvegarde un GeoTIFF multi-bandes.

    ARI = (1/B03 - 1/B05) / (1/B03 + 1/B05)
    """

    # Ouvrir les rasters
    band3_ds = gdal.Open(band3_path)
    band5_ds = gdal.Open(band5_path)

    if band3_ds is None or band5_ds is None:
        raise RuntimeError("Impossible d'ouvrir les fichiers d'entrée")

    # Vérifier le nombre de bandes
    nb_bands = band3_ds.RasterCount
    if band5_ds.RasterCount != nb_bands:
        raise ValueError("B03 et B05 doivent avoir le même nombre de bandes")

    rows = band3_ds.RasterYSize
    cols = band3_ds.RasterXSize

    # Créer le raster de sortie (multi-bandes)
    driver = gdal.GetDriverByName("GTiff")
    out_ds = driver.Create(
        output_path,
        cols,
        rows,
        nb_bands,
        gdal.GDT_Float32,
        options=["COMPRESS=LZW"]
    )

    # Projection et géotransformation
    out_ds.SetGeoTransform(band3_ds.GetGeoTransform())
    out_ds.SetProjection(band3_ds.GetProjection())

    # Calcul ARI bande par bande (temps)
    for i in range(1, nb_bands + 1):
        b3 = band3_ds.GetRasterBand(i).ReadAsArray().astype(np.float32)
        b5 = band5_ds.GetRasterBand(i).ReadAsArray().astype(np.float32)

        # Masque des valeurs invalides
        mask = (b3 <= 0) | (b5 <= 0)

        ari = np.full(b3.shape, -9999, dtype=np.float32)

        ari_valid = (1.0 / b3 - 1.0 / b5) / (1.0 / b3 + 1.0 / b5)
        ari[~mask] = ari_valid[~mask]

        out_band = out_ds.GetRasterBand(i)
        out_band.WriteArray(ari)
        out_band.SetNoDataValue(-9999)
        out_band.FlushCache()

    # Nettoyage
    band3_ds = None
    band5_ds = None
    out_ds = None

    print("✅ Fichier ARI_serie_temp.tif créé avec succès (série temporelle multi-bandes)")



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



def compute_class_statistics(image, mask, classes, nodata=-9999):
    """
    Calcule la moyenne et l'écart-type par classe et par date.

    Parameters
    ----------
    image : ndarray (rows, cols, dates)
        Série temporelle (ex: ARI)
    mask : ndarray (rows, cols)
        Raster des classes
    classes : list
        Liste des labels de classes à analyser
    nodata : float
        Valeur nodata à exclure

    Returns
    -------
    means : ndarray (nb_classes, nb_dates)
    stds : ndarray (nb_classes, nb_dates)
    """

    nb_dates = image.shape[2]
    nb_classes = len(classes)

    means = np.zeros((nb_classes, nb_dates))
    stds = np.zeros((nb_classes, nb_dates))

    for i, cls in enumerate(classes):
        cls_mask = mask == cls

        for d in range(nb_dates):
            values = image[:, :, d][cls_mask]
            values = values[values != nodata]

            if values.size > 0:
                means[i, d] = np.mean(values)
                stds[i, d] = np.std(values)
            else:
                means[i, d] = np.nan
                stds[i, d] = np.nan

    return means, stds



def merge_and_save_multiband(
    raster_paths,
    output_path,
    reference_raster_path,
    nodata=None,
    dtype=gdal.GDT_Float32
):
    """
    Fusionne plusieurs fichiers GeoTIFF multibandes en un seul GeoTIFF multibande.

    Parameters
    ----------
    raster_paths : list of str
        Liste des chemins vers les fichiers GeoTIFF à fusionner
    output_path : str
        Chemin du GeoTIFF de sortie
    reference_raster_path : str
        Raster de référence pour la géoréférence
    nodata : float, optional
        Valeur nodata
    dtype : gdal data type
        Type des données en sortie
    """

    images = []

    # Lecture de toutes les bandes de tous les fichiers
    for path in raster_paths:
        ds = gdal.Open(path)
        if ds is None:
            raise IOError(f"Impossible d'ouvrir le fichier : {path}")

        rows = ds.RasterYSize
        cols = ds.RasterXSize
        nb_bands = ds.RasterCount

        for b in range(1, nb_bands + 1):
            band = ds.GetRasterBand(b)
            arr = band.ReadAsArray()
            images.append(arr)

        ds = None  # fermeture GDAL

    # Vérification des dimensions
    shapes = [img.shape for img in images]
    if len(set(shapes)) != 1:
        raise ValueError("Toutes les bandes doivent avoir les mêmes dimensions")

    # Empilement -> (rows, cols, bands)
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
        out_band = out_ds.GetRasterBand(b + 1)
        out_band.WriteArray(multiband_image[:, :, b])

        if nodata is not None:
            out_band.SetNoDataValue(nodata)

        out_band.FlushCache()

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

def rasterize_vector(vector_path, ref_raster_path, output_path, attribute="strate"):
    ds_ref = gdal.Open(ref_raster_path)
    if ds_ref is None:
        raise IOError("Impossible d'ouvrir le raster de référence")

    gt = ds_ref.GetGeoTransform()
    proj = ds_ref.GetProjection()
    cols = ds_ref.RasterXSize
    rows = ds_ref.RasterYSize

    driver = gdal.GetDriverByName("GTiff")
    ds_out = driver.Create(output_path, cols, rows, 1, gdal.GDT_Byte)
    ds_out.SetGeoTransform(gt)
    ds_out.SetProjection(proj)

    band = ds_out.GetRasterBand(1)
    band.SetNoDataValue(0)
    band.Fill(0)

    vect_ds = ogr.Open(vector_path)
    layer = vect_ds.GetLayer()

    gdal.RasterizeLayer(ds_out, [1], layer, options=[f"ATTRIBUTE={attribute}"])

    ds_out = None
    vect_ds = None
    ds_ref = None

    print(f"✅ Rasterisation terminée : {output_path}")
    return output_path

