from math import modf

def dms2dd(dms):
    minsec, degrees = modf(abs(dms))
    seconds, minutes = modf(minsec * 100)
    dd = degrees + minutes / 60 + seconds / 36
    return dd if dms >= 0 else -dd