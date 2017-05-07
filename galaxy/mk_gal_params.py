
cluster = dict(
    save=False,
    verbose=False,
    mstar_min=5e9,
    den_lim=1e6,
    den_lim2=5e6,
    rmin=-1,
    Rgal_to_reff=4.0,
    method_com=1,
    follow_bp=None,
    unit_conversion='code',
    debug=False)

HAGN = dict(
    save=False,
    verbose=False,
    mstar_min=1e10,
    den_lim=1e6,
    den_lim2=5e6,
    rmin=-1,
    Rgal_to_reff=4.0,
    method_com=1,
    follow_bp=None,
    unit_conversion='code',
    debug=False)

NH = dict(
    save=False,
    verbose=False,
    mstar_min=1e8,
    den_lim=1e6,
    den_lim2=5e6,
    rmin=-1,
    Rgal_to_reff=5.0,
    method_com=1,
    follow_bp=None,
    unit_conversion='code',
    debug=False)

available = tuple(k for k in locals() if not k.startswith('_'))
