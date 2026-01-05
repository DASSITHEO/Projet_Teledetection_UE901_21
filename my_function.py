from osgeo import gdal
import numpy as np

def raster_to_float32_nodata(input_raster, nodata_value=-9999):
    """
    Lit un raster et retourne ses valeurs converties en float32.

    Parameters
    ----------
    input_raster : str
        Chemin du raster d'entrée
    nodata_value : float
        Valeur NoData à utiliser si présente

    Returns
    -------
    array : np.ndarray
        Tableau numpy en float32
    ds : gdal.Dataset
        Dataset GDAL (pour récupérer CRS, résolution, etc.)
    """

    ds = gdal.Open(input_raster)
    band = ds.GetRasterBand(1)

    array = band.ReadAsArray().astype("float32")

    src_nodata = band.GetNoDataValue()
    if src_nodata is not None:
        array[array == src_nodata] = nodata_value

    return array, ds


"""
def image_to_float32_nodata(input_raster, nodata_value=-9999):
    
    #Lit un raster et convertit la bande 1 en float32 et remplace les nodata en -9999.
    #Retourne un seul objet contenant les données et les métadonnées.

    ds = gdal.Open(input_raster)
    band = ds.GetRasterBand(1)

    array = band.ReadAsArray().astype(np.float32)

    src_nodata = band.GetNoDataValue()
    if src_nodata is not None:
        array[array == src_nodata] = nodata_value

    # On attache le tableau converti au dataset
    ds.array = array
    ds.nodata = nodata_value

    return ds """
