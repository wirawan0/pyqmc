# Create: 20090112

import math

class avg(object):
  def __init__(self):
    self.clear()
  def clear(self):
    self.N = 0
    self.sum1 = 0.0
    self.sum2 = 0.0
  def add(self, value):
    self.N += 1
    self.sum1 += value
    self.sum2 += value*value
    return self
  def mean(self):
    return self.sum1 / self.N
  def stddev(self):
    N = self.N
    avg1 = self.sum1 / N
    avg2 = self.sum2 / N
    return math.sqrt(abs( (avg2 - avg1**2) * N / (N - 1) ))

# Operator overloading:
avg.__iadd__ = avg.add
avg.__call__ = avg.mean
