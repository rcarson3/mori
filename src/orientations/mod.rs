// This file is a part of the mori - Material Orientation Library in Rust
// Copyright 2018 Robert Carson

//    Licensed under the Apache License, Version 2.0 (the "License");
//    you may not use this file except in compliance with the License.
//    You may obtain a copy of the License at

//        http://www.apache.org/licenses/LICENSE-2.0

//    Unless required by applicable law or agreed to in writing, software
//    distributed under the License is distributed on an "AS IS" BASIS,
//    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//    See the License for the specific language governing permissions and
// limitations under the License.

use super::*;

pub mod bunge;
pub mod rmat;
pub mod ang_axis;
pub mod rod_vec;
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

pub trait RotVector{
    fn rot_vector(&self, vec: ArrayView2<f64>) -> Array2<f64>;
    fn rot_vector_mut(&self, vec: ArrayView2<f64>, rvec: ArrayViewMut2<f64>);
    fn rot_vector_inplace(&self, vec: ArrayViewMut2<f64>);
}
pub trait RotTensor{
    fn rot_tensor(&self, tensor: ArrayView3<f64>) -> Array3<f64>;
    fn rot_tensor_mut(&self, tensor: ArrayView3<f64>, rtensor: ArrayViewMut3<f64>);
    fn rot_tensor_inplace(&self, tensor: ArrayViewMut3<f64>);
}


