
# FreeSurfer Python Package

The FS python package provides an easier way interact with overlay, image, volume and surface objects. Classes are written in pure python, but some IO and computationally intensive routines are written in c++ and depend on the `utils` freesurfer library. The package can be imported like so:

```python
import freesurfer as fs
```

### ND Arrays

N-Dimensional arrays, between 1 and 3, are represented by `Overlay`, `Image`, and `Volume` classes respectively. These classes all inherit from the `ArrayContainerTemplate` base class, which defines some general array utilities.

`Overlay`, `Image`, and `Volume` objects all contain a `data` member that points to the internal `np.ndarray`. In most cases, `data.ndim` will equal the number dimensions represented by the array type. However, an extra dimension is also allowed - this represents the number of data frames. These objects are constructed from a `np.ndarray` that *does not* get copied by default. To create a volume from an array:

```python
array = np.zeros((256, 256, 256), dtype='uint8')
vol = fs.Volume(array)
vol.data is array  # returns True
```

Array containers can also be loaded from file with the `read()` static member function. As a simple example, to load a `uchar` volume and convert it to a `float` volume with values scaled between 0 and 1:

```python
vol = fs.Volume.read('norm.mgz')
vol.data = vol.data / 255
vol.write('scaled.mgz')
```

### Surfaces

3D meshes are represented by `Surface` objects, and the surface topology is defined entirely by two individual vertex and face arrays:

```python
surf = fs.Surface.read('lh.pial')
surf.vertices  # [nvertices x 3] float array of vertex positions
surf.faces     # [nfaces x 3] int array defining the vertices in each triangle
```

New surfaces can be constructed directly from vertex and face arrays:

```python
surf = fs.Surface(vertices, faces)
surf.write('lh.pial')
```

### Linear Transforms

Linear transforms contain a 4 x 4 affine matrix and can transform a point or series of points between spaces.

```python
xform = fs.LinearTransform.read('reg.lta')
xform.matrix  # [4 x 4] float array 

points = xform.transform(points)
```
