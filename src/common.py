'''Common functions for large projects'''
import os
import string
import numpy as np

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

def closest_point(points,point):
    points = np.asarray(points)
    deltas = points - point
    dist_2 = np.einsum('ij,ij->i', deltas, deltas)
    return np.argmin(dist_2)

