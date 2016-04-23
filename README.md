# Channel DNS Mean Flow Loader

[data source](http://turbulence.ices.utexas.edu/) 

```python
y = np.linspace(0.0, 0.5, 201)
loader = DNSDataLoader(Re, y)
data = loader.data
print data.keys()
# loaded keys
['y', 'y+', 'u+', 'ub+', 'vb+', 'wb+', '-om_z+', 'om_xb+', 'om_yb+', 'om_zb+', 'uvb+', 'uwb+', 'vwb+', 'prb+', 'psb+', 'pstob+', 'pb', 'k+', 'dissip', 'prodoc', 'p-strain', 'p-diff', 't-diff', 'v-diff', 'bal', 'tp-kbal', 'wilcox_y', 'wilcox_unknown_1', 'wilcox_unknown_2', 'wilcox_u+', 'wilcox_uvb+', 'wilcox_k+', 'wilcox_dissip', 'wilcox_produc', 'wilcox_y+']
```