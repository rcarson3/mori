// This file is a part of the mori - Material Orientation Library in Rust
// Copyright 2018 Robert Carson
//
//    Licensed under the Apache License, Version 2.0 (the "License");
//    you may not use this file except in compliance with the License.
//    You may obtain a copy of the License at
//
//        http://www.apache.org/licenses/LICENSE-2.0
//
//    Unless required by applicable law or agreed to in writing, software
//    distributed under the License is distributed on an "AS IS" BASIS,
//    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//    See the License for the specific language governing permissions and
//    limitations under the License.

//! # mori - Material Orientation Library in Rust
//! An orientation library built around commonly used orientation representations used in crystallography and engineering applications. It contains conversion, rotation, and data analysis operations for various orientation spaces.
//!
//! Orientations play a large role in a number of fields ranging from: crystallography, x-ray diffraction, metallurgy, solid mechanics, and the list can go on and on. Therefore, it's important to have a library that easily allows conversions from Eulearian representations to rotation matrices to neo-eulerian representation to quaternions. The initial scope of this library will provide common sets of conversions. In an attempt to have a consistent set of conversions with others in the field, a majority of these conversions have been taken from <sup name="a1">[1](#f1)</sup> Later as it develops the plan is to include the following:
//!
//! * Crystallographic fundamental region conversions
//! * Mean orientation calculations
//! * Misorientation calculations based upon the work of <sup name="a2">[2](#f2)</sup> <sup name="a3">[3](#f3)</sup>
//!
//!
//! <b id="f1">[1]:</b>[↩](#a1) D Rowenhorst et al 2015 Consistent representations of and conversions between 3D rotations Modelling Simul. Mater. Sci. Eng. 23 083501
//!
//! <b id="f2">[2]:</b>[↩](#a2) Barton N R and Dawson P R 2001 A methodology for determining average lattice orientation and its application to the characterization of grain substructure Metall. Mater. Trans. A 32 1967–75
//!
//! <b id="f3">[3]:</b>[↩](#a3) Glez J C and Driver J 2001 Orientation distribution analysis in deformed grains J. Appl. Cryst. 34 280–8
//!
//! # Evolution of crate
//! Besides the listed above capabilities the following features are still being looked at to extend it:
//!
//! * Parallel capabilities
//! * In-place conversion processes of orientations conversions so conversions can go into an already existing memory structure.
//! * Various quaternion operations that aren't necessarily related to typical orientations operations but are useful.
//!

#[macro_use]
extern crate ndarray;
use ndarray::prelude::*;

///Contains orientation conversions from one to another, rotations of vectors and tensor data, and finally relative orientation operations
pub mod orientations;
//Future module for us to use. This will require us to use 
//pub mod par_orientations;