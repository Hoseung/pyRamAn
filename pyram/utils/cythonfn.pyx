# cython: language_level = 3

def update_field_vals(field, totalweight, a, b, c, nx, ny, value, average):
    """ This updates the field array (and the totweight array if
    average is True).
    
    The elements to update and their values are inferred from
    a,b,c and value.
    """
    nxny = nx*ny
    indices = a['ind'] + b['ind'] * nx + c['ind'] * nxny
    weights = a['weight'] * b['weight'] * c['weight']
    value = weights * value
    if average:
        for i,ind in enumerate(indices):
            field[ind] += value[i]
            totalweight[ind] += weights[i]
    else:
        for i, ind in enumerate(indices):
            field[ind] += value[i]
