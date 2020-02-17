import numpy as np
import struct

def load_header_double(brick_data):
    nbodies = struct.unpack("i", brick_data[4:8])[0]
    massp = struct.unpack("d", brick_data[16:24])[0]
    aexp = struct.unpack("d", brick_data[32:40])[0]
    omegat = struct.unpack("d", brick_data[48:56])[0]
    age = struct.unpack("d", brick_data[64:72])[0]
    halnum = struct.unpack("i", brick_data[80:84])[0]
    subnum = struct.unpack("i", brick_data[84:88])[0]
    print(nbodies, halnum, subnum, "header done")
    return 96

def load_header_single(brick_data):
    nbodies = struct.unpack("i", brick_data[4:8])[0]
    massp = struct.unpack("f", brick_data[16:20])[0]
    aexp = struct.unpack("f", brick_data[28:32])[0]
    omegat = struct.unpack("f", brick_data[40:44])[0]
    age = struct.unpack("f", brick_data[52:56])[0]
    halnum = struct.unpack("i", brick_data[64:68])[0]
    subnum = struct.unpack("i", brick_data[68:72])[0]
    print(nbodies, halnum, subnum, "header done")
    return 80, halnum, subnum


def load_a_halo_single(brick_data, offset, is_gal=True):
    np = struct.unpack("i", brick_data[offset:offset+4])[0]
    offset += 12  # 12 = 4 + 8
    ids = struct.unpack_from("<{}i".format(np), brick_data[offset:offset+4*np])
    offset += 4*np + 8
    hid = struct.unpack("i", brick_data[offset:offset+4])[0]
    offset += 12
    nstep = struct.unpack("i", brick_data[offset:offset+4])[0]
    offset += 12
    hhost = struct.unpack_from("<5i", brick_data[offset:offset+20])
    offset += 28
    mass = struct.unpack("f", brick_data[offset:offset+4])[0]
    offset += 12
    pos = struct.unpack_from("<3f", brick_data[offset:offset+12])
    offset += 20
    vel = struct.unpack_from("<3f", brick_data[offset:offset+12])
    offset += 20
    ang = struct.unpack_from("<3f", brick_data[offset:offset+12])
    offset += 20
    radius = struct.unpack_from("<4f", brick_data[offset:offset+16])
    offset += 24
    energy = struct.unpack_from("<3f", brick_data[offset:offset+12])
    offset += 20
    spin = struct.unpack("f", brick_data[offset:offset+4])[0]
    offset += 12
    if is_gal:
        gal_props = struct.unpack_from("<3f", brick_data[offset:offset+12])
        offset += 20
    vir = struct.unpack_from("<4f", brick_data[offset:offset+16])
    offset += 24
    profile = struct.unpack_from("<2f", brick_data[offset:offset+8])
    offset += 16
    if is_gal:
        g_nbin = struct.unpack("i", brick_data[offset:offset+4])[0]
        offset += 12
        g_rr = struct.unpack_from("<{}f".format(g_nbin), brick_data[offset:offset+g_nbin*4])
        offset += 8 + g_nbin*4
        g_rho = struct.unpack_from("<{}f".format(g_nbin), brick_data[offset:offset+g_nbin*4])
        offset += 8 + g_nbin*4

    return offset


def load_a_halo_double(brick_data, offset, is_gal=True):
    np = struct.unpack("i", brick_data[offset:offset+4])[0]
    offset += 12  # 12 = 4 + 8
    ids = struct.unpack_from("<{}i".format(np), brick_data[offset:offset+4*np])
    offset += 4*np + 8
    hid = struct.unpack("i", brick_data[offset:offset+4])[0]
    offset += 12
    nstep = struct.unpack("i", brick_data[offset:offset+4])[0]
    offset += 12
    hhost = struct.unpack_from("<5i", brick_data[offset:offset+20])
    offset += 28
    mass = struct.unpack("d", brick_data[offset:offset+8])[0]
    offset += 16
    pos = struct.unpack_from("<3d", brick_data[offset:offset+24])
    offset += 32
    vel = struct.unpack_from("<3d", brick_data[offset:offset+24])
    offset += 32
    ang = struct.unpack_from("<3d", brick_data[offset:offset+24])
    offset += 32
    radius = struct.unpack_from("<4d", brick_data[offset:offset+32])
    offset += 40
    energy = struct.unpack_from("<3d", brick_data[offset:offset+24])
    offset += 32
    spin = struct.unpack("d", brick_data[offset:offset+8])[0]
    offset += 16
    if is_gal:
        gal_props = struct.unpack_from("<3d", brick_data[offset:offset+24])
        offset += 32
    vir = struct.unpack_from("<4d", brick_data[offset:offset+32])
    offset += 40
    profile = struct.unpack_from("<2d", brick_data[offset:offset+16])
    offset += 24
    if is_gal:
        g_nbin = struct.unpack("i", brick_data[offset:offset+4])[0]
        offset += 12
        g_rr = struct.unpack_from("<{}d".format(g_nbin), brick_data[offset:offset+g_nbin*8])
        offset += 8 + g_nbin*8
        g_rho = struct.unpack_from("<{}d".format(g_nbin), brick_data[offset:offset+g_nbin*8])
        offset += 8 + g_nbin*8

    return offset
