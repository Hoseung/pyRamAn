
# Also in treeutils
def get_field_view(arr, fields):
    '''
    Returns a view to a 2D numpy array with multiple fields.
    By default, veiw = ndarray['field'] returns a view,
    but ndarray['field1','field2'] COPY data.
    Example:
    >>>
    '''
    import numpy as np
    dtype2 = np.dtype({name:arr.dtype.fields[name] for name in fields})
    return np.ndarray(arr.shape, dtype2, arr, 0, arr.strides)

def remove_field_name(a, name):
    names = list(a.dtype.names)
    if name in names:
        names.remove(name)
    b = a[names]
    return b