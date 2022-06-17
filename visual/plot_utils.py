#!/usr/bin/env python
import numpy as np
import pvf_settings as ps

figsizes = [(10, 10.5), (10, 10.5), (10, 6.5), (8, 10.5), (10, 6.5), (14, 6.5)]
figrwcls = [(2, 2), (1, 1), (1, 2), (2, 1), (1, 2), (1, 3)]
figplace = [(3, 2, 0), (0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 0), (1, 0, 0)]
cbsplace = [(1, 1), (-1, -1), (2, 2), (0, 1), (1, 1), (2, 2)]


def plane_in_outputname(figmode, draw2D):
    to_insert = ''
    if figmode == 1 or figmode == 4:
        if draw2D[0]:
            to_insert = 'yz_'
        elif draw2D[1]:
            to_insert = 'xz_'
        elif draw2D[2]:
            to_insert = 'xy_'
    return to_insert


def fsym(vmin, vmax):
    vmx = np.max([np.abs(vmin), np.abs(vmax)])
    vmn = -1.0 * vmx
    return vmn, vmx


def execute_comparison(orig, comp, ctype):
    if ctype == 1:
        return orig - comp
    elif ctype == 2:
        return orig / comp
    elif ctype == 3:
        return (orig / comp) - 1


def scale_manage(sctype, refis, umin, umax, d1, d2, extr):
    symmin = 1.0
    autoscale = False
    d1min, d1max, d2min, d2max, d3min, d3max = min(extr[0]), max(extr[1]), min(extr[2]), max(extr[3]), min(extr[4]), max(extr[5])
    if d2:
        dmin, dmax = d2min, d2max
    else:
        dmin, dmax = d1min, d1max
    if (umin == 0.0 and umax == 0.0):
        vmin, vmax = dmin, dmax
        autoscale = True
    else:
        vmin, vmax = umin, umax

    if (sctype == '1' or sctype == 'symlin'):
        vmin, vmax = fsym(vmin, vmax)

    elif (sctype == '2' or sctype == 'log'):
        if (vmin > 0.0):
            vmin = np.log10(vmin)
        else:
            vmin = np.log10(check_minimum_data(refis))
        if (vmax > 0.0):
            vmax = np.log10(vmax)
        else:
            vmax = ps.fineqv
        if (vmin == np.inf):
            vmin = vmax - ps.fineqv
    elif (sctype == '3' or sctype == 'symlog'):
        if (umin > 0.0 and umax > 0.0):
            symmin = umin
            vmax = np.log10(umax / symmin)
        else:
            if (dmin * dmax > 0.0):
                symmin, smax = sorted([np.abs(dmin), np.abs(dmax)])
            else:
                symmin, smax = check_extremes_absdata(refis)
            if (smax == -np.inf or symmin == np.inf):
                smax = 10.
                symmin = 1.
            vmax = np.log10(smax / symmin)
        vmin = -vmax
        print('SYMMIN value for SYMLOG scaletype: %s' % symmin)

    if vmin == vmax:
        vmin = vmin - ps.fineqv
        vmax = vmax + ps.fineqv

    print('3D data value range: ', d3min, d3max)
    if d2:
        print('Slices  value range: ', d2min, d2max)
    if d1:
        print('1D data value range: ', d1min, d1max)
    print('Plotted value range: ', vmin, vmax)
    return vmin, vmax, symmin, autoscale


def check_minimum_data(refis):
    cmdmin = np.inf
    for blks in refis:
        for bl in blks:
            bxyz, binb = bl[0:2]
            for ncut in range(3):
                if binb[ncut]:
                    cmdmin = min(cmdmin, np.min(bxyz[ncut], initial=np.inf, where=(bxyz[ncut] > 0.0)))
    return cmdmin


def check_extremes_absdata(refis):
    cmdmin, cmdmax = np.inf, -np.inf
    for blks in refis:
        for bl in blks:
            bxyz, binb = bl[0:2]
            for ncut in range(3):
                if binb[ncut]:
                    cmdmin = min(cmdmin, np.min(np.abs(bxyz[ncut]), initial=np.inf, where=(np.abs(bxyz[ncut]) > 0.0)))
                    cmdmax = max(cmdmax, np.max(np.abs(bxyz[ncut]), initial=-np.inf, where=(np.abs(bxyz[ncut]) > 0.0)))
    return cmdmin, cmdmax


def scale_plotarray(pa, sctype, symmin):
    if (sctype == '2' or sctype == 'log'):
        pa = np.log10(pa)
    elif (sctype == '3' or sctype == 'symlog'):
        pa = np.sign(pa) * np.log10(np.maximum(np.abs(pa) / symmin, 1.0))
    return pa


def list3_subtraction(l3, s3):
    return l3[0] - s3[0], l3[1] - s3[1], l3[2] - s3[2]


def list3_division(l3, divisor):
    return l3[0] / divisor, l3[1] / divisor, l3[2] / divisor


def list3_max(l3, r3):
    return [max(l3[0], r3[0]), max(l3[1], r3[1]), max(l3[2], r3[2])]


def list3_min(l3, r3):
    return [min(l3[0], r3[0]), min(l3[1], r3[1]), min(l3[2], r3[2])]


def labelx():
    return lambda var: '$' + str(var)[2:-1].replace('**', '^') + '$'


def labellog(sctype, symmin, cmpr0):
    logname = ''
    compare = ''
    if cmpr0:
        compare = '(compared) '
    if (sctype == '2' or sctype == 'log'):
        logname = 'log '
    elif (sctype == '3' or sctype == 'symlog'):
        logname = '{symmetry level = %8.3e} log ' % symmin
    return logname + compare


def take_nonempty(lst):
    for it in lst:
        if it != []:
            return it
    return []


def colorbar_mode(drawd, drawh, figmode):
    if drawd and drawh and figmode == 3:
        cbar_mode = 'each'
    elif drawd and drawh and figmode != 3:
        cbar_mode = 'none'
    elif drawd or drawh:
        cbar_mode = 'single'
    else:
        cbar_mode = 'none'
    return cbar_mode


def color_axes(wax, color):
    wax.spines['top'].set_color(color)
    wax.spines['bottom'].set_color(color)
    wax.spines['left'].set_color(color)
    wax.spines['right'].set_color(color)
    wax.tick_params(axis='x', colors=color)
    wax.tick_params(axis='y', colors=color)
    return


def detindex(nd, cxyz, smin, smax):
    return int(np.floor(nd * (cxyz - smin) / (smax - smin)))


def ind_limits(nd, cxyz, smin, smax):
    return int(min(nd - 1, max(0, detindex(nd, cxyz, smin, smax))))


def isinbox(cxyz, smin, smax, warn, cc):
    isin = (cxyz >= smin and cxyz <= smax)
    if not isin and warn:
        print('Domain edges %s %s used to plot as the given plot center %s coordinate (%s) is outside the domain.' % (smin, smax, cc, cxyz))
    return isin


def find_indices(nd, cxyz, smin, smax, warn):
    inb = isinbox(cxyz[0], smin[0], smax[0], warn, 'CX'), isinbox(cxyz[1], smin[1], smax[1], warn, 'CY'), isinbox(cxyz[2], smin[2], smax[2], warn, 'CZ')
    icc = ind_limits(nd[0], cxyz[0], smin[0], smax[0]), ind_limits(nd[1], cxyz[1], smin[1], smax[1]), ind_limits(nd[2], cxyz[2], smin[2], smax[2])
    return inb, icc


def check_plotlevels(plotlevels, maxglev, drawa, toplot):
    if plotlevels == '':
        if drawa:
            plotlevels = range(maxglev + 1)
        else:
            plotlevels = 0,
    else:
        npl = []
        for il in plotlevels:
            if il in range(maxglev + 1):
                npl.append(il)
            else:
                print('LEVEL %s not met in the file!' % str(il))
        plotlevels = npl
    if toplot:
        print('Levels to plot: ', plotlevels)
    else:
        print('Compare levels: ', plotlevels)
    return plotlevels


def sanitize_gridlist(gridlist, cgcount):
    if gridlist == '':
        gridlist = range(cgcount)
    else:
        grdl = []
        for ig in gridlist:
            if ig >= 0 and ig < cgcount:
                grdl.append(ig)
            else:
                print('Grid block %s does not exist.' % str(ig))
        gridlist = grdl
        print('GRIDLIST: ', gridlist)
    return gridlist


def reorder_gridcolorlist(gcolor, maxglev, plotlevels):
    gaux1 = gcolor.split(',')
    maxgc = len(gaux1)
    gaux2 = []
    ail = 0
    for il in range(maxglev + 1):
        if il in plotlevels:
            gaux2.append(gaux1[ail % maxgc])
            ail += 1
        else:
            gaux2.append('none')
    return gaux2


def linestyles(markers, maxglev, plotlevels):
    if markers == '':
        markers = ps.plot1d_linestyle
    marknum = len(markers)
    linstyl = []
    for il in range(maxglev + 1):
        if il in plotlevels:
            linstyl.append(markers[-1])
        else:
            linstyl.append('none')
    for im in range(marknum):
        if len(plotlevels) > im:
            linstyl[plotlevels[-(im + 1)]] = markers[im]
    return linstyl


def choose_amr_or_uniform(drawa, drawu, drawd, drawg, drawp, maxglev, gridlist):
    if (drawd or (drawg and not drawp)) and not (drawa or drawu):
        if maxglev == 0 and gridlist == '':
            drawu = True
        else:
            drawa = True
    return drawa, drawu


def check_1D2Ddefaults(axc, n_d, double_cbar):
    draw1D, draw2D = axc
    if not (any(draw1D) or any(draw2D)):
        if sum(n_d) - 2 == n_d[0] * n_d[1] * n_d[2]:
            draw1D = True, True, True
        else:
            draw2D = True, True, True
    figmode = 1
    if any(draw2D):
        if all(draw2D) or (draw2D[0] and draw2D[2] and not draw2D[1]):
            figmode = 0
        elif draw2D[0] and draw2D[1] and not draw2D[2]:
            figmode = 2
        elif draw2D[2] and draw2D[1] and not draw2D[0]:
            figmode = 3
    if double_cbar and (figmode == 1 or figmode == 2):
        figmode = figmode + 3
    return draw1D, draw2D, figmode


def convert_units(infile, toplot):
    au_cm = 1.49597870700e13
    pc_au = 206264.806248712
    pc_cm = pc_au * au_cm
    if infile == 'pc':
        cm = 1.0 / pc_cm
        pc = 1.0
    elif infile == 'au':
        cm = 1.0 / au_cm
        pc = pc_cm * cm
    elif infile == 'kpc':
        cm = 1.0 / (1.0e3 * pc_cm)
        pc = 0.001
    elif infile == 'm':
        cm = 1.0 / 1.0e2
        pc = pc_cm * cm
    else:
        return 1., False

    if toplot == 'cm':
        return cm, True
    elif toplot == 'metr':
        return 1.0e2 * cm, True
    elif toplot == 'km':
        return 1.0e5 * cm, True
    elif toplot == 'au':
        return au_cm * cm, True
    elif toplot == 'pc':
        return pc, True
    elif toplot == 'kpc':
        return 1.0e3 * pc, True
    elif toplot == 'Mpc':
        return 1.0e6 * pc, True
    elif toplot == 'lyr':
        return 9.4605e17 * cm, True
    else:
        return 1., False


def change_units(fromfile, toplot):
    infile = fromfile.decode('utf-8')
    if infile == toplot or toplot == '':
        return 1., fromfile, False
    if toplot == 'k':
        return 1.e3, b"".join([b'k', fromfile]), True
    if toplot == 'M':
        return 1.e6, b"".join([b'M', fromfile]), True
    if toplot == 'G':
        return 1.e9, b"".join([b'G', fromfile]), True
    if toplot == 'm':
        return 1.e-3, b"".join([b'm', fromfile]), True
    if toplot == 'mu':
        return 1.e-6, b"".join([b'mu', fromfile]), True
    conv, chan = convert_units(infile, toplot)
    if chan:
        return conv, bytes(toplot, 'utf-8'), chan
    else:
        return conv, fromfile, chan
