import numpy as np
import math

from . import errorExit


def MRISP(scale, nfuncs):
  if scale == 0: scale = 1
  cols = int(round(scale * 256))
  rows = 2 * cols
  return np.zeros((cols, rows, nfuncs))


def averageRadius(surface):

  xhi = yhi = zhi = -10000.0
  xlo = ylo = zlo =  10000.0
  
  for (x, y, z) in surface[0]:
    if x > xhi: xhi = x
    if x < xlo: xlo = x
    if y > yhi: yhi = y
    if y < ylo: ylo = y
    if z > zhi: zhi = z
    if z < zlo: zlo = z

  x0 = (xlo + xhi) / 2.0
  y0 = (ylo + yhi) / 2.0
  z0 = (zlo + zhi) / 2.0

  radius = 0.0
  for (x, y, z) in surface[0]:
    dx = x - x0;
    dy = y - y0;
    dz = z - z0;
    radius += math.sqrt(dx**2 + dy**2 + dz**2)

  radius /= len(surface[0])
  return radius


def parameterizeSurface(surface, overlay, fno=0, mrisp=None, scale=1.0):

  # compute the average radius
  a = b = c = averageRadius(surface)

  if mrisp is None:
    mrisp = MRISP(scale, 1)

  mrisp.fill(0.)
  cols = mrisp.shape[0]
  rows = mrisp.shape[1]

  fillstatus = np.zeros((cols, rows))
  # enum for fillstatus
  unfilled = 0
  filling = 1
  filled = 2

  distances = np.zeros((cols, rows))
  coords = []

  # calculate total distances to a point in parameter space
  for (x, y, z) in surface[0]:
    theta = math.atan2(y/b, x/a)
    if theta < 0: theta += 2 * math.pi
    d = c**2 - z**2
    if d < 0: d = 0
    phi = math.atan2(math.sqrt(d), z)

    uf = cols * phi / math.pi
    vf = rows * theta / (2 * math.pi)
    u = int(round(uf))
    v = int(round(vf))

    if u < 0: u = -u
    if u >= cols: u = (2*cols) - (u + 1)
    if v < 0: v += rows
    if v >= rows: v -= rows

    fillstatus[u][v] = filled
    distances[u][v] += 1  # keep track of total # of nodes
    coords.append((u, v))

  # now add in curvatures proportional to their distance from the point
  for i, (u, v) in enumerate(coords):
    totalDistance = distances[u][v]
    if totalDistance > 0.0: mrisp[u][v][fno] += overlay[i] / float(totalDistance)

  # fill in values which were unmapped using soap bubble
  tmpimage = np.zeros((cols, rows, 1))
  npasses = 0
  while True:
    num_unfilled = 0
    for u in range(cols):
      for v in range(rows):
        if fillstatus[u][v] == unfilled:
          total = 0.0
          n = 0
          for uk in [-1, 0, 1]:
            u1 = u + uk
            if u1 < 0: u1 = -u1
            if u1 >= cols: u1 = (2*cols) - (u1 + 1)
            for vk in [-1, 0, 1]:
              v1 = v + vk
              if v1 < 0: v1 += rows
              if v1 >= rows: v1 -= rows
              if fillstatus[u1][v1] == filled:
                total += mrisp[u1][v1][fno]
                n += 1
          if n > 0:
            total /= float(n)  # TODO float might not be necessary
            tmpimage[u][v][fno] = total
            fillstatus[u][v] = filling
          else:
            num_unfilled += 1
        else:
          tmpimage[u][v][fno] = mrisp[u][v][fno]
    fillstatus[fillstatus == filling] = filled
    mrisp[:] = tmpimage
    npasses += 1
    if npasses > 1000: errorExit('could not fill parameterization')
    if num_unfilled == 0: break

  return mrisp
