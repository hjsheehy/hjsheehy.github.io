'''Common functions for large projects'''
import os
import string
import numpy as np
import random

def GetCh():
    """Read single character from standard input without echo."""
    import sys, tty, termios
    fd = sys.stdin.fileno()
    old_settings = termios.tcgetattr(fd)
    try:
        tty.setraw(sys.stdin.fileno())
        ch = sys.stdin.read(1)
    finally:
        termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)
    return ch

def YesNo(question):
    c = ""
    print(question + " [Y/n]: ", end = "", flush = True)
    while c not in ("y", "Y", "n", "N"):
        c = GetCh().lower()
    return c == 'y' or c == "Y"

def CreateFolder(folder):
    '''If folder does not exist, create one.'''
    if os.path.exists(folder):
        pass
    else:
        print(f'{folder}\nNo folder found. Created one.')
        os.system(f'mkdir {folder}')

def ClearOutFolder(out_folder, silent=False):
    '''If out_folder does not exists, creates one. Prompts clearing
    of out_folder contents if not-empty.'''
    CreateFolder(out_folder)
    if len(os.listdir(out_folder)) != 0:
        if not silent:
            if YesNo(f'Clear out folder before creating new .conf files located at\n{out_folder}?'):
                pass
            else:
                exit()
        os.system(f'rm -r {out_folder}*')

suffixes = ['B', 'KB', 'MB', 'GB', 'TB', 'PB']

def humansize(nbytes):
    '''Returns human readable memory'''
    i = 0
    while nbytes >= 1024 and i < len(suffixes)-1:
        nbytes /= 1024.
        i += 1
    f = ('%.2f' % nbytes).rstrip('0').rstrip('.')
    return '%s %s' % (f, suffixes[i])

def list_depth(lst):
    return isinstance(lst, list) and max(map(list_depth, lst)) + 1

def tuple_to_tuples(lst):
    """If lst has depth 1, returns [lst]"""
    if list_depth(lst)==1:
        lst=[lst]
    return lst

def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero.

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower offset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax / (vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highest point in the colormap's range.
          Defaults to 1.0 (no upper offset). Should be between
          `midpoint` and 1.0.
    '''
    import matplotlib

    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    matplotlib.pyplot.register_cmap(cmap=newcmap)

    return newcmap

def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])

def points_in_circle(radius, x0=0, y0=0):
    radius-=1
    x_ = np.arange(x0 - radius - 1, x0 + radius + 1, dtype=int)
    y_ = np.arange(y0 - radius - 1, y0 + radius + 1, dtype=int)
    x, y = np.where((x_[:,np.newaxis] - x0)**2 + (y_ - y0)**2 <= radius**2)
    pts = np.array([x,y]).T
    pts = pts - radius - 1
    return pts

def points_on_circle(radius, x0=0, y0=0):
    pts=points_in_circle(radius=radius, x0=x0, y0=y0)
    pts_inside=points_in_circle(radius=radius-1, x0=x0, y0=y0)
    pts=pts[np.all(np.any((pts-pts_inside[:, None]), axis=2), axis=0)]
    return pts
