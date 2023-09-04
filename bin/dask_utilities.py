#! /usr/bin/env python3

import dask.array as da


def slicewise_pearson_r(Y_true, Y_pred, axis: int):
    """
    Pearson correlation among slices of Y_true and Y_pred along an axis.
    """

    assert Y_true.ndim == 2
    assert Y_pred.ndim == 2
    assert axis in [0, 1]
    assert Y_true.shape == Y_pred.shape

    Y_true = Y_true - da.mean(Y_true, axis=axis, keepdims=True)
    Y_pred = Y_pred - da.mean(Y_pred, axis=axis, keepdims=True)

    r = da.sum(Y_true * Y_pred, axis=axis) / da.sqrt(
        da.sum(Y_true**2, axis=axis) * da.sum(Y_pred**2, axis=axis)
    )

    return r
