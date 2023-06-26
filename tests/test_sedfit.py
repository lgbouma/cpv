from complexrotators.sedfit import run_SED_analysis

trim_ticids = {
    # ticid: each entry in form (>xmin, <xmax, >ymin, <ymax)
    # e.g. to exclude everything about 1e-9 erg/cm2/s:
    '2234692': [ (None, None, 1e-9, None) ],
}

#ticid = "2234692"
ticid = '142173958'
ticid = '89463560'
ticid = '234295610'
ticid = '177309964'

trimlist = None
if ticid in trim_ticids:
    trimlist = trim_ticids[ticid]

run_SED_analysis(ticid, trimlist=trimlist)
