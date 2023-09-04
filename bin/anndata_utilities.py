#! /usr/bin/env python3

import dask.array as da
import zarr
from anndata.experimental import read_dispatched, read_elem


def read_dask_anndata(store: str):
    """
    Read a zarr store as an AnnData object using dask arrays. The normal
    adata.read_zarr method loads everything in memory.

    Adapted from
    https://anndata.readthedocs.io/en/latest/tutorials/notebooks/%7Bread%2Cwrite%7D_dispatched.html
    """
    f = zarr.open(store, mode="r")

    def callback(func, elem_name, elem, iospec):
        if iospec.encoding_type in (
            "dataframe",
            "csr_matrix",
            "csc_matrix",
            "awkward-array",
        ):
            return read_elem(elem)
        elif iospec.encoding_type == "array":
            return da.from_zarr(elem)
        else:
            return func(elem)

    adata = read_dispatched(f, callback=callback)

    return adata
