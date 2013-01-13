#!/usr/bin/env python

import h5py as h5

def readH5(fname):
    """
    Reads an HDF5 file to a dictionary. 
    """
    try:
        fid = h5.File(fname)
    except:
        print("Can't read file [%s]"%fname)
    groups = []
    for k in fid.keys():
        tmp = {}
        g = fid[k]
        for kk in g.keys():
            print("Reading %s" % kk)
            try:
                tmp[kk] = g[kk].value
            except AttributeError:
                gg = g[kk]
                tmp2 = {}
                for kkk in gg.keys():
                    print("Reading %s" % kkk)
                    try:
                        tmp2[kkk] = gg[kkk].value
                    except AttributeError: 
                        print("Could not read [%s], too deep."%kkk)
                tmp[kk] = tmp2
        groups.append(tmp)
    return groups
