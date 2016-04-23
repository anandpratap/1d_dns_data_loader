# Channel DNS Mean Flow Loader

[data source](http://turbulence.ices.utexas.edu/) 

```python
y = np.linspace(0.0, 0.5, 201)
loader = DNSDataLoader(Re, y)
data = loader.data
print data.keys()
```