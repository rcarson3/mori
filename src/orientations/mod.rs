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

//! # Module: orientation
//!
//! This submodule contains methods related to the conversion of orientations from one to another.
//! It also contains methods directly related to orientations that are able to rotate vectors and 2nd order tensors.
//! It allows easy addition of new orientations in which one only has to add at least one conversion to an already existing orientation method.
//! Afterwards, all of the other orientations representations can be reached from one of the already existing orientation representations.

use super::*;
///It holds the necessary tools to create/view an array of Bunge angles. Also, the conversion processes to other orientation representations.
pub mod bunge;
///It holds the necessary tools to create/view an array of rotation matrices. Also, the conversion processes to other orientation representations.
///Finally, it holds the tools necessary to rotate matrices and 2nd order tensors.
pub mod rmat;
///It holds the necessary tools to create/view an array of axis-angle representations. Also, the conversion processes to other orientation representations.
pub mod ang_axis;
///It holds the necessary tools to create/view an array of Rodrigues vectors. Also, the conversion processes to other orientation representations.
pub mod rod_vec;
///It holds the necessary tools to create/view an array of unit quaterions. Also, the conversion processes to other orientation representations.
pub mod quat;

pub use self::bunge::*;
pub use self::rmat::*;
pub use self::ang_axis::*;
pub use self::rod_vec::*;
pub use self::quat::*;


#[derive(Clone, Debug)]
pub struct Homochoric{
    ori: Array2<f64>,
}

///A set of generic orientation conversions from one to another orientation representation 
pub trait OriConv{
    fn to_bunge(&self) -> Bunge;
    fn to_rmat(&self) -> RMat;
    fn to_ang_axis(&self) -> AngAxis;
    fn to_ang_axis_comp(&self) -> AngAxisComp;
    fn to_rod_vec(&self) -> RodVec;
    fn to_rod_vec_comp(&self) -> RodVecComp;
    fn to_quat(&self) -> Quat;
    fn to_homochoric(&self) -> Homochoric;
}
///A set of methods that allow for rotations of supplied vector data (3x1) dim
pub trait RotVector{
    fn rot_vector(&self, vec: ArrayView2<f64>) -> Array2<f64>;
    fn rot_vector_mut(&self, vec: ArrayView2<f64>, rvec: ArrayViewMut2<f64>);
    fn rot_vector_inplace(&self, vec: ArrayViewMut2<f64>);
}
///A set of methods that allow for rotations of supplied 2nd order tensor data (3x3) dim
pub trait RotTensor{
    fn rot_tensor(&self, tensor: ArrayView3<f64>) -> Array3<f64>;
    fn rot_tensor_mut(&self, tensor: ArrayView3<f64>, rtensor: ArrayViewMut3<f64>);
    fn rot_tensor_inplace(&self, tensor: ArrayViewMut3<f64>);
}


