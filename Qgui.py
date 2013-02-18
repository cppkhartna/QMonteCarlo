import mpl_toolkits.mplot3d.axes3d as p3
import pylab as p
import os
from ctypes import *

class replica(Structure):
    pass
replica._fields_ = [
        ("x", POINTER(c_double)),
        ("next", POINTER(replica))
        ]

ROOT = os.path.dirname(os.path.abspath(__file__))
qlib = CDLL('%s/qlib.so' % ROOT)
run = qlib.run
run.restype = POINTER(replica)
atom = True
ion = False
replicas = run(c_int(4000), c_int(1000), ion)

fig=p.figure()
ax = p3.Axes3D(fig)
ax.autoscale(True)

xs = []
ys = []
zs = []

rep = replicas.contents
#NULL = POINTER(replica)()
try:
    while True:
        i = rep.x
        xs.append(i[0])
        ys.append(i[1])
        zs.append(i[2])
        rep = rep.next.contents
except ValueError:
    pass

ax.scatter3D(xs, ys, zs)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
p.show()
