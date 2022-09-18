# Vector-Field-Topolgy-2D

<table border="0">
 <tr>
    <td><b style="font-size:30px">Python version</b></td>
    <td><b style="font-size:30px">Matlab version</b></td>
 </tr>
 <tr>
    <td><p float="left">
<a href = "https://github.com/zaman13/Vector-Field-Topolgy-2D/tree/master/Codes"> <img src="https://img.shields.io/badge/Language-Python-blue" alt="alt text"> </a>
<a href = "https://github.com/zaman13/Vector-Field-Topolgy-2D/blob/master/LICENSE"> <img src="https://img.shields.io/github/license/zaman13/Poisson-solver-2D" alt="alt text"></a>
<a href = "https://github.com/zaman13/Vector-Field-Topolgy-2D/tree/master/Codes"> <img src="https://img.shields.io/badge/version-1.0-red" alt="alt text"> </a>
</p></td>
    <td><p float="left">
<a href = "https://github.com/zaman13/Vector-Field-Topolgy-2D/tree/master/Codes"> <img src="https://img.shields.io/badge/Language-Matlab-blue" alt="alt text"> </a>
<a href = "https://github.com/zaman13/Vector-Field-Topolgy-2D/blob/master/LICENSE"> <img src="https://img.shields.io/github/license/zaman13/Poisson-solver-2D" alt="alt text"></a>
<a href = "https://github.com/zaman13/Vector-Field-Topolgy-2D/tree/master/Codes"> <img src="https://img.shields.io/badge/version-1.6-red" alt="alt text"> </a>
</p></td>
 </tr>
</table>


<img align = "right" src="https://github.com/zaman13/Vector-Field-Topolgy-2D/blob/master/Sample%20output/Python%20output/out3.svg" alt="alt text" width="480">

Vector field topology (VFT) [Ref. 1, 2] is an intuitive visualization technique of vector fields. It is a representation of the global topology based on an analysis of the critical points and their connections. For a complex vector field, it can be difficult to understand the field properties by plotting streamlines or field lines (quiver). By showing only the critical points and key integral lines connecting those points, the VFT technique produces a minimalistic intuitive representation of the field.

VFT methods have been used in many fields including fluid flow analysis, computer vision and photonics.  

The algorithm was initially implemented in MATLAB. A python version has been added since. The two versions are not equivalent. The python implementation contains more upgraded algorithms for identifying the critical points. 

## Usage
Please cite the papers mentioned in the reference section and this github repository if you use the code for your research/work.




## Sample output
### Python version <a href = "https://github.com/zaman13/Vector-Field-Topolgy-2D/tree/master/Codes/Python%20codes/Version%201.0"> <img src="https://img.shields.io/badge/Language-Python-blue" alt="alt text"> </a>

For a test vector filed u<sub>x</sub> = sin(X+Y)exp(-0.5X<sup>2</sup> - 0.3Y<sup>2</sup>), u<sub>y</sub> = cos(X-Y)exp(-0.5X<sup>2</sup> - 0.3Y<sup>2</sup>), the calculated vector field topology is shown below (with quiver lines and without quiver lines):

<p float="left">
<img src="https://github.com/zaman13/Vector-Field-Topolgy-2D/blob/master/Sample%20output/Python%20output/out2.svg" alt="alt text" width="400">
<img src="https://github.com/zaman13/Vector-Field-Topolgy-2D/blob/master/Sample%20output/Python%20output/out1.svg" alt="alt text" width="400">
</p>




The colors on the integral lines correspond to the magnitude of the field.

### Matlab version <a href = "https://github.com/zaman13/Vector-Field-Topolgy-2D/tree/master/Codes/Matlab%20codes/Version%201.6"> <img src="https://img.shields.io/badge/Language-Matlab-blue" alt="alt text"> </a>

For a test vector filed u<sub>x</sub> = -X<sup>2</sup> + Y<sup>2</sup>, u<sub>y</sub> = X<sup>2</sup> + Y<sup>2</sup> - 2, the calculated vector field topology is shown below (with quiver lines and without quiver lines):

<p float="left">
<img src="https://github.com/zaman13/Vector-Field-Topolgy-2D/blob/master/Sample%20output/test_field_1.svg" alt="alt text" width="400">
<img src="https://github.com/zaman13/Vector-Field-Topolgy-2D/blob/master/Sample%20output/test_field_1_no_quiver.svg" alt="alt text" width="400">
</p>

## Theory

### Integral Lines

The integral lines follow the direction of the vector field, denoted here by <b>u</b> = [u<sub>x</sub> u<sub>y</sub>]<sup>T</sup>. For simplicity, we will assume that u is normalized. Let's consider a point on the integral line denoted by the coordinates (x<sub>0</sub>, y<sub>0</sub>). Considering the locus of the integral line is formed by a discrete set of points (x<sub>0</sub>, y<sub>0</sub>), (x<sub>1</sub>, y<sub>1</sub>), (x<sub>2</sub>, y<sub>2</sub>), ..., (x<sub>n</sub>, y<sub>n</sub>). Then,

x<sub>1</sub> = x<sub>0</sub> + u<sub>x</sub>(x<sub>0</sub>,y<sub>0</sub>) Δx

x<sub>2</sub> = x<sub>1</sub> + u<sub>x</sub>(x<sub>1</sub>,y<sub>1</sub>) Δx = x<sub>0</sub> + u<sub>x</sub>(x<sub>0</sub>,y<sub>0</sub>) Δx +  u<sub>x</sub>(x<sub>1</sub>,y<sub>1</sub>) Δx

.
.

x<sub>n</sub> = x<sub>0</sub> + &sum;<sub>i=0</sub><sup>n-1</sup> u<sub>x</sub>(x<sub>i</sub>,y<sub>i</sub>) Δx
<!-- <sub>x<sub>0</sub></sub><sup>x</sup> -->
In the continuous limit,
x = x<sub>0</sub> + ∫ u<sub>x</sub> dx'




## References
1. J. Helman, and L. Hesselink. "Representation and display of vector field topology in fluid flow data sets." Computer,  vol. 22, pp. 27-36, August 1989.
2. J. Helman, and L. Hesselink. "Visualizing vector field topology in fluid flows." IEEE Computer Graphics and Applications, vol. 11, pp. 36-46, May/June 1991.
