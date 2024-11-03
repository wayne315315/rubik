import itertools
from collections import defaultdict

import numpy as np

### List all 54 coordinates
axes = (0,1,2)
entries = (-1,0,1)
coords = np.array([(*coord, n) for n, *coord in itertools.product(axes, entries, entries, entries) if coord[n] != 0])
corner, = np.where(np.all(coords[:,:-1] != 0, axis=-1))
center, = np.where(np.sum(coords[:,:-1] == 0, axis=-1) == 2)
edge = np.asarray([i for i in range(len(coords)) if i not in set(corner)|set(center)])
index = np.asarray([i for i in range(len(coords))])

# block <-> index mapping
b2i = defaultdict(list)
i2b = {}
for i, (*b, n) in enumerate(coords):
    b2i[tuple(b)].append(i)
    i2b[i] = tuple(b)
b2i = dict(b2i)

# coordinate <-> index mapping
c2i = {tuple(coord): i for i, coord in enumerate(coords)}
i2c = {i: tuple(coord) for i, coord in enumerate(coords)}

### Coordination system for Rubik's cube (x,y,z,n) ; x,y,z belongs to {-1,0,1}, normal vector n belongs to {0,1,2}
### Operator system for Rubik's cube rotation (Xi, Yj, Zk)
# axis : 0 -> x axis ; 1 -> y axis ; 2 -> z axis
# level : the level of the plane of rotation, can be -1, 0, or 1
def _rotate(coords, axis, level):
    # setup the unit vector of the rotation axis
    ax = [0,0,0]
    ax[axis] = 1

    # any face within the plane of rotation can be categorized into 2 group based on their normal vector
    # one tangent to the rotation plane, the other perpendicular to the plane
    # cross product the axis vector with the given coordinate to obtain the location after 90 degree rotation
    crosses = np.cross(ax, coords[:,:3])
    crosses[:, axis] = level
    ns = np.array([[n] if n == axis else list({0,1,2} - {n,axis}) for n in coords[:,-1]])
    crosses = np.hstack([crosses, ns])

    # any face out of the plane of rotation remains unchanged
    coords = np.where(coords[:, axis, np.newaxis] != level, coords, crosses)
    return coords

def rotate(coords, axis, level, orient):
    stack = [(axis, level, orient)]
    seq = []
    while stack:
        axis, level, orient = stack.pop()
        if level == 0:
            stack.append((axis, 1, -orient))
            stack.append((axis, -1, -orient))
        elif orient == -1:
            for _ in range(3):
                seq.append((axis, level))
        else:
            seq.append((axis, level))

    seq = seq[::-1]
    for axis, level in seq:
        coords = _rotate(coords, axis, level)
    return coords

### Create Rotation class to customize string representation
class Rotation:
    names = ["".join(item) for item in itertools.product(("x", "y", "z"), ("p","0", "n"), ("p","n"))]
    pairs = list(itertools.product((0,1,2), (1,0,-1), (1,-1)))
    n2p = dict(zip(names, pairs))
    p2n = dict(zip(pairs, names))
   
    global coords, corner, edge, center, index, b2i, i2b, c2i, i2c
    def __init__(self, axis, level, orient):
        self.axis = axis
        self.level = level
        self.orient = orient
        self.indices = tuple([c2i[tuple(c)] for c in rotate(coords, axis, level, orient)])

    
    def __repr__(self):
        name = self.p2n[(self.axis, self.level, self.orient)]
        return name
    
    def __lt__(self, other):
        return str(self) < str(other)

    def __call__(self, item):
        if item.shape != coords.shape and item.shape != (len(coords),):
            raise ValueError("item needs to be either coordinates or indices") 
        indices = np.asarray([c2i[tuple(c)] for c in item]) if item.shape == coords.shape else item.copy()
        indices = np.asarray([self.indices[indices[i]] for i in range(len(indices))])
        item = coords[indices] if item.shape == coords.shape else indices

        return item
    
    @classmethod
    def fromstring(cls, name):
        try:
            assert name in cls.names
        except AssertionError as err:
            raise AssertionError("'%s' is not in  %s" % (name, cls.names))
        axis, level = cls.n2p[name]
        return Rotation(axis, level)


class RotationSequence:
    def __init__(self, seq):
        self.seq = tuple(seq)
        indices = index.copy()
        for r in seq:
            indices = r(indices)
        self.indices = tuple(indices)

    def __repr__(self):
        name = " -> ".join([str(r) for r in self.seq])
        return name

    def __lt__(self, other):
        return str(self) < str(other)

    def __eq__(self, other):
        return self.indices == other.indices

    def __call__(self, item):
        if item.shape != coords.shape and item.shape != (len(coords),):
            raise ValueError("item needs to be either coordinates or indices") 
        indices = np.asarray([c2i[tuple(c)] for c in item]) if item.shape == coords.shape else item.copy()
        indices = np.asarray([self.indices[indices[i]] for i in range(len(indices))])
        item = coords[indices] if item.shape == coords.shape else indices
        return item

### List all 6 rotations (x+, x-, y+, y-, z+, z-)
rs = [Rotation(axis, level, orient) for axis, level, orient in Rotation.pairs]
xpp, xpn, x0p, x0n, xnp, xnn, ypp, ypn, y0p, y0n, ynp, ynn, zpp, zpn, z0p, z0n, znp, znn = rs
