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

use super::*;

///A structure that holds an array of axis-angle representation of a rotation
#[derive(Clone, Debug)]
pub struct AngAxis{
    pub ori: Array2<f64>,
}
///A structure that holds an array of compact axis-angle representation of a rotation
#[derive(Clone, Debug)]
pub struct AngAxisComp{
    ori: Array2<f64>
}

impl AngAxis{

    ///Creates an array of zeros for the initial axis-angle parameterization when data is not fed into it
    pub fn new(size: usize) -> AngAxis{
        assert!(size > 0, "Size inputted: {}, was not greater than 0", size);

        let mut ori = Array2::<f64>::zeros((4, size).f());

        azip!(mut angaxis (ori.axis_iter_mut(Axis(1))) in {angaxis[3] = 1.0_f64});

        AngAxis{
            ori,
        }
    }//End of new

    ///Creates an axis-angle parameterization  type with the supplied data as long as the supplied data is in the following format
    ///shape (3, nelems), memory order = fortran/column major.
    ///If it doesn't fit those standards it will fail.
    pub fn new_init(ori: Array2<f64>) -> AngAxis{

        let nrow = ori.rows();

        assert!(nrow == 4, "Number of rows of array was: {}, which is not equal to 4", nrow);
        //We need to deal with a borrowing of ori here, so we need to have strides dropped at one point.
        {
            let strides = ori.strides();

            assert!(strides[0] == 1, "The memory stride is not column major (f order)");
        }

        AngAxis{
            ori,
        }
    }//End of new_init

    ///Return a ndarray view of the orientation data
    pub fn ori_view(&self) -> ArrayView2<f64>{
        self.ori.view()
    }

    //Return a ndarray mutable view of the orientation data
    pub fn ori_view_mut(&mut self) -> ArrayViewMut2<f64>{
        self.ori.view_mut()
    }
}//End of AngAxis impl

///The orientation conversions of a series of axis-angle representation to a number of varying different orientation
///representations commonly used in material orientation processing. 
impl OriConv for AngAxis{
    ///Converts the axis-angle representation over to a rotation matrix which has the following properties
    ///shape (3, 3, nelems), memory order = fortran/column major.
    fn to_rmat(&self) -> RMat{

        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array3::<f64>::zeros((3, 3, nelems).f());

        azip!(mut rmat (ori.axis_iter_mut(Axis(2))), ref angaxis (self.ori.axis_iter(Axis(1))) in {
            let c = angaxis[3].cos();
            let s = angaxis[3].sin();

            rmat[[0, 0]] = c + (1.0_f64 - c) * (angaxis[0] * angaxis[0]);
            rmat[[1, 0]] = (1.0_f64 - c) * (angaxis[0] * angaxis[1]) + s * angaxis[2];
            rmat[[2, 0]] = (1.0_f64 - c) * (angaxis[0] * angaxis[2]) - s * angaxis[1];

            rmat[[0, 1]] = (1.0_f64 - c) * (angaxis[0] * angaxis[1]) - s * angaxis[2];
            rmat[[1, 1]] = c + (1.0_f64 - c) * (angaxis[1] * angaxis[1]);
            rmat[[2, 1]] = (1.0_f64 - c) * (angaxis[1] * angaxis[2]) + s * angaxis[0];

            rmat[[0, 2]] = (1.0_f64 - c) * (angaxis[0] * angaxis[2]) + s * angaxis[1];
            rmat[[1, 2]] = (1.0_f64 - c) * (angaxis[1] * angaxis[2]) - s * angaxis[0];
            rmat[[2, 2]] = c + (1.0_f64 - c) * (angaxis[2] * angaxis[2]);
        });

        RMat::new_init(ori)
    }//End of to_rmat

    ///Converts the axis-angle representation over to Bunge angles which has the following properties
    ///shape (3, nelems), memory order = fortran/column major.
    fn to_bunge(&self) -> Bunge{
        let rmat = self.to_rmat();
        //When a pure conversion doesn't exist we just use the already existing ones in other orientation
        //representations 
        rmat.to_bunge()
    }//End of to_bunge

    ///Returns a clone of the axis-angle structure
    fn to_ang_axis(&self) -> AngAxis{
        self.clone()
    }

    ///Converts the axis-angle representation over to a compact angle-axis representation which has the following properties
    ///shape (3, nelems), memory order = fortran/column major.
    fn to_ang_axis_comp(&self) -> AngAxisComp{
        //We first convert to a axis-angle representation. Then we scale our normal vector by our the rotation
        //angle which is the fourth component of our axis-angle vector.
        // let ang_axis = self.to_ang_axis();

        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((3, nelems).f());

        azip!(mut angaxis_comp (ori.axis_iter_mut(Axis(1))), ref angaxis (self.ori.axis_iter(Axis(1))) in {
            angaxis_comp[0] = angaxis[0] * angaxis[3];
            angaxis_comp[1] = angaxis[1] * angaxis[3];
            angaxis_comp[2] = angaxis[2] * angaxis[3];
        });

        AngAxisComp::new_init(ori)
    }//End of to_ang_axis_comp

    ///Converts the axis-angle representation over to a Rodrigues vector representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    fn to_rod_vec(&self) -> RodVec{
        //We first convert to a axis-angle representation. Then we just need to change the last component
        //of our axis-angle representation to be tan(phi/2) instead of phi
        // let ang_axis = self.to_ang_axis();

        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((4, nelems).f());

        let inv2 = 1.0_f64/2.0_f64;

        azip!(mut rodvec (ori.axis_iter_mut(Axis(1))), ref angaxis (self.ori.axis_iter(Axis(1))) in {
            rodvec[0] = angaxis[0];
            rodvec[1] = angaxis[1];
            rodvec[2] = angaxis[2];
            rodvec[3] = f64::tan(inv2 * angaxis[3]);
        });

        RodVec::new_init(ori)
    }//End of to_rod_vec

    ///Converts the axis-angle representation over to a compact Rodrigues vector representation which has the following properties
    ///shape (3, nelems), memory order = fortran/column major.
    fn to_rod_vec_comp(&self) -> RodVecComp{
        //We first convert to a Rodrigues vector representation. Then we scale our normal vector by our the rotation
        //angle which is the fourth component of our axis-angle vector.
        //If we want to be more efficient about this in the future with out as many copies used we can reuse a lot of the code
        //used in the to_ang_axis code. However, we will end up with a lot of similar/repeated code then. We could put that
        //code in a helper function that isn't seen.
        let rod_vec = self.to_rod_vec();
        rod_vec.to_rod_vec_comp()
    }//End of to_rod_vec_comp

    ///Converts the axis-angle representation over to a unit quaternion representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    fn to_quat(&self) -> Quat{

        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((3, nelems).f());

        let inv2 = 1.0_f64 / 2.0_f64;

        azip!(mut quat (ori.axis_iter_mut(Axis(1))), ref angaxis (self.ori.axis_iter(Axis(1))) in {
            let s = f64::sin(inv2 * angaxis[3]); 

            quat[0] = f64::cos(inv2 * angaxis[3]);
            quat[1] = s * angaxis[0];
            quat[2] = s * angaxis[1];
            quat[3] = s * angaxis[2];
        });

        Quat::new_init(ori)
    }//End of to_quat

    ///Converts the axis-angle representation over to a homochoric representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    fn to_homochoric(&self) ->Homochoric{

        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((4, nelems).f());

        let inv3  = 1.0_f64 / 3.0_f64;
        let inv34 = 3.0_f64 / 4.0_f64;

        azip!(mut homoch (ori.axis_iter_mut(Axis(1))), ref angaxis (self.ori.axis_iter(Axis(1))) in {
            let pow_term    = inv34 * (angaxis[3] - angaxis[3].sin()); 

            homoch[0] = angaxis[0];
            homoch[1] = angaxis[1];
            homoch[2] = angaxis[2];
            homoch[3] = pow_term.powf(inv3);
        });

        Homochoric{
            ori,
        }
    }//End of to_homochoric
}//End of Impl OriConv for AngAxis

impl AngAxisComp{

    ///Creates an array of zeros for the initial compact axis-angle parameterization when data is not fed into it
    pub fn new(size: usize) -> AngAxisComp{
        assert!(size > 0, "Size inputted: {}, was not greater than 0", size);

        let mut ori = Array2::<f64>::zeros((4, size).f());

        azip!(mut angaxis (ori.axis_iter_mut(Axis(1))) in {angaxis[3] = 1.0_f64});

        AngAxisComp{
            ori,
        }
    }//End of new

    ///Creates a compact axis-angle parameterization  type with the supplied data as long as the supplied data is in the following format
    ///shape (3, nelems), memory order = fortran/column major.
    ///If it doesn't fit those standards it will fail.
    pub fn new_init(ori: Array2<f64>) -> AngAxisComp{

        let nrow = ori.rows();

        assert!(nrow == 4, "Number of rows of array was: {}, which is not equal to 4", nrow);
        //We need to deal with a borrowing of ori here, so we need to have strides dropped at one point.
        {
            let strides = ori.strides();

            assert!(strides[0] == 1, "The memory stride is not column major (f order)");
        }

        AngAxisComp{
            ori,
        }
    }//End of new_init

    ///Return a ndarray view of the orientation data
    pub fn ori_view(&self) -> ArrayView2<f64>{
        self.ori.view()
    }

    ///Return a ndarray mutable view of the orientation data
    pub fn ori_view_mut(&mut self) -> ArrayViewMut2<f64>{
        self.ori.view_mut()
    }
}//End of AngAxisComp impl

///The orientation conversions of a series of compact axis-angle representation to a number of varying different orientation
///representations commonly used in material orientation processing. 
impl OriConv for AngAxisComp{
    ///Converts the compact axis-angle representation over to a rotation matrix which has the following properties
    ///shape (3, 3, nelems), memory order = fortran/column major.
    fn to_rmat(&self) -> RMat{
        let ang_axis = self.to_ang_axis();
        ang_axis.to_rmat()
    }//End of to_rmat

    ///Converts the compact axis-angle representation over to Bunge angles which has the following properties
    ///shape (3, nelems), memory order = fortran/column major.
    fn to_bunge(&self) -> Bunge{
        let rmat = self.to_rmat();
        rmat.to_bunge()
    }//End of to_bunge

    ///Converts the compact axis-angle representation over to an angle-axis representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    fn to_ang_axis(&self) -> AngAxis{
        //We first convert to a axis-angle representation. Then we scale our normal vector by our the rotation
        //angle which is the fourth component of our axis-angle vector.
        // let ang_axis = self.to_ang_axis();

        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((4, nelems).f());

        azip!(mut angaxis (ori.axis_iter_mut(Axis(1))), ref angaxis_comp (self.ori.axis_iter(Axis(1))) in {
            let norm_angaxis = f64::sqrt({
                angaxis_comp[0] * angaxis_comp[0] 
                + angaxis_comp[1] * angaxis_comp[1] 
                + angaxis_comp[2] * angaxis_comp[2]
                });
            let inv_norm_angaxis = 1.0_f64 / norm_angaxis;

            angaxis[0] = angaxis_comp[0] * inv_norm_angaxis;
            angaxis[1] = angaxis_comp[1] * inv_norm_angaxis;
            angaxis[2] = angaxis_comp[2] * inv_norm_angaxis;
            angaxis[3] = norm_angaxis;
        });

        AngAxis::new_init(ori)
    }

    ///Returns a clone of the compact axis-angle structure
    fn to_ang_axis_comp(&self) -> AngAxisComp{
        self.clone()
    }//End of to_ang_axis_comp

    ///Converts the compact axis-angle representation over to a Rodrigues vector representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    fn to_rod_vec(&self) -> RodVec{
        let ang_axis = self.to_ang_axis();
        ang_axis.to_rod_vec()
    }//End of to_rod_vec

    ///Converts the compact axis-angle representation over to a compact Rodrigues vector representation which has the following properties
    ///shape (3, nelems), memory order = fortran/column major.
    fn to_rod_vec_comp(&self) -> RodVecComp{
        //We first convert to a Rodrigues vector representation. Then we scale our normal vector by our the rotation
        //angle which is the fourth component of our axis-angle vector.
        //If we want to be more efficient about this in the future with out as many copies used we can reuse a lot of the code
        //used in the to_ang_axis code. However, we will end up with a lot of similar/repeated code then. We could put that
        //code in a helper function that isn't seen.
        let rod_vec = self.to_rod_vec();
        rod_vec.to_rod_vec_comp()
    }//End of to_rod_vec_comp

    ///Converts the compact axis-angle representation over to a unit quaternion representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    fn to_quat(&self) -> Quat{
        let ang_axis = self.to_ang_axis();
        ang_axis.to_quat()
    }//End of to_quat

    ///Converts the compact axis-angle representation over to a homochoric representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    fn to_homochoric(&self) ->Homochoric{
        let ang_axis = self.to_ang_axis();
        ang_axis.to_homochoric()
    }//End of to_homochoric
}//End of Impl OriConv for AngAxisComp