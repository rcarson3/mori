# mori
An orientation library built around commonly used orientation representations used in crystallography and engineering applications. It contains conversion, rotation, and data analysis operations for various orientation spaces.

Orientations play a large role in a number of fields ranging from: crystallography, x-ray diffraction, metallurgy, solid mechanics, and the list can go on and on. Therefore, it's important to have a library that easily allows conversions from Eulearian representations to rotation matrices to neo-eulerian representation to quaternions. The initial scope of this library will provide common sets of conversions. In an attempt to have a consistent set of conversions with others in the field, a majority of these conversions have been taken from <sup name="a1">[1](#f1)</sup> Later as it develops the plan is to include the following:

*Crystallographic fundamental region conversions
*Mean orientation calculations
*Misorientation calculations based upon the work of <sup name="a2">[2](#f2)</sup> <sup name="a3">[3](#f3)</sup>


<b id="f1">1</b> [1]:[↩](#a1)
<b id="f2">1</b>[2]:[↩](#a2) Barton N R and Dawson P R 2001 A methodology for determining average lattice orientation and its application to the characterization of grain substructure Metall. Mater. Trans. A 32 1967–75
<b id="f3">1</b>[3]:[↩](#a3) Glez J C and Driver J 2001 Orientation distribution analysis in deformed grains J. Appl. Cryst. 34 280–8
