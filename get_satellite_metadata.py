"""
Created on August 31, 2019

@author: Yi Wang
"""

from netCDF4 import Dataset

def get_mxd04_metadata(filename):
    """Get Metadata from MxD04 file

    Returns:
    --------
    """

    f = Dataset(filename, 'r')

    # read CoreMetadata.0
    core_md = f.getncattr('CoreMetadata.0')
    core_md = core_md.split('\n')

    # read ArchiveMetada.0
    arc_md = f.getncattr('ArchiveMetadata.0')
    arc_md = arc_md.split('\n')

    f.close()

    # dictionary to save metadata
    metadata = {}

    # GRINGPOINTLONGITUDE
    i = arc_md.index('        OBJECT                 = GRINGPOINTLONGITUDE')
    tmp = arc_md[i+3]
    tmp = tmp.split()
    GRINGPOINTLONGITUDE = (
        float(tmp[2][1:-1]),
        float(tmp[3][:-1]),
        float(tmp[4][:-1]),
        float(tmp[5][:-1]),
        )
    metadata['GRINGPOINTLONGITUDE'] = GRINGPOINTLONGITUDE

    # GRINGPOINTLATITUDE
    i = arc_md.index('        OBJECT                 = GRINGPOINTLATITUDE')
    tmp = arc_md[i+3]
    tmp = tmp.split()
    GRINGPOINTLATITUDE = (
        float(tmp[2][1:-1]),
        float(tmp[3][:-1]),
        float(tmp[4][:-1]),
        float(tmp[5][:-1]),
        )
    metadata['GRINGPOINTLATITUDE'] = GRINGPOINTLATITUDE

    return metadata
