using Test
using BasisFunctions
using Winston


## We create spline bases on for the function space [5,15] -> Real with evenly-spaced nodes, and
## apply them to the grid x
x = linspace(5,15,101)
X1 = linearSplineBFE(x, (5,15), 11)
X2 = quadraticSplineBFE(x, (5,15), 11)
X3 = cubicSplineBFE(x, (5,15), 11)

p = FramedPlot()
map(i -> add(p, Curve(x, X1[:,i], "color", "black")), 1:3)
Winston.display(p)

p = FramedPlot()
map(i -> add(p, Curve(x, X2[:,i], "color", "red")), 1:5)
Winston.display(p)

p = FramedPlot()
map(i -> add(p, Curve(x, X3[:,i], "color", "blue")), 1:7)
Winston.display(p)


## Now we use a pre-specified set of nodes
u = [8.05,11.34,12.31,14.53, 14.75]
X1 = splineBFE(1, x, u)
X2 = splineBFE(2, x, u)

p = FramedPlot()
map(i -> add(p, Curve(x, X1[:,i], "color", "black")), 1:3)
Winston.display(p)

p = FramedPlot()
map(i -> add(p, Curve(x, X2[:,i], "color", "black")), 1:2)
Winston.display(p)
