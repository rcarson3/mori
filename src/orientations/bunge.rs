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

///A struct that holds an array of Bunge angles.
#[derive(Clone, Debug)]
pub struct Bunge{
    ori: Array2<f64>,
}

impl Bunge{

    ///Creates an array of zeros for the initial Bunge angles when data is not fed into it
    pub fn new(size: usize) -> Bunge{
        assert!(size > 0, "Size inputted: {}, was not greater than 0", size);

        let ori = Array2::<f64>::zeros((3, size).f());

        Bunge{
            ori,
        }
    }//End of new

    ///Creates a Bunge type with the supplied data as long as the supplied data is in the following format
    ///shape (3, nelems), memory order = fortran/column major.
    ///If it doesn't fit those standards it will fail.
    pub fn new_init(ori: Array2<f64>) -> Bunge{

        let nrow = ori.rows();

        assert!(nrow == 3, "Number of rows of array was: {}, which is not equal to 3", nrow);
        //We need to deal with a borrowing of ori here, so we need to have strides dropped at one point.
        {
            let strides = ori.strides();

            assert!(strides[0] == 1, "The memory stride is not column major (f order)");
        }

        Bunge{
            ori,
        }
    }//End of new_init

    ///Return a ndarray view of the orientation data
    pub fn ori_view(&self) -> ArrayView2<f64>{
        self.ori.view()
    }

    ///Return a ndarray mutable view of the orienation data
    pub fn ori_view_mut(&mut self) -> ArrayViewMut2<f64>{
        self.ori.view_mut()
    }
}//End of Bunge impl

///The orientation conversions of a series of Bunge angles to a number of varying different orientation
///representations commonly used in material orientation processing. 
impl OriConv for Bunge{
    ///The conversion from Bunge angles to Bungle angles just returns a copy of the original data structure.
    fn to_bunge(&self) -> Bunge{
        self.clone()
    }//End of to_bunge

    ///Converts the Bunge angles over to a rotation matrix which has the following properties
    ///shape (3, 3, nelems), memory order = fortran/column major.
    fn to_rmat(&self) -> RMat{

        let nelems = self.ori.len_of(Axis(1));
        
        let mut ori = Array3::<f64>::zeros((3, 3, nelems).f());

        azip!(mut rmat (ori.axis_iter_mut(Axis(2))), ref bunge (self.ori.axis_iter(Axis(1))) in {
            let s1 = bunge[0].sin();
            let c1 = bunge[0].cos();
            let s2 = bunge[1].sin();
            let c2 = bunge[1].cos();
            let s3 = bunge[2].sin();
            let c3 = bunge[2].cos();

            rmat[[0, 0]] = c1 * c3 - s1 * s3 * c2;
            rmat[[1, 0]] = -c1 * s3 - s1 * c2 * c3;
            rmat[[2, 0]] = s1 * s2;

            rmat[[0, 1]] = s1 * c3 + c1 * c2 * s3;
            rmat[[1, 1]] = -s1 * s3 + c1 * c2 * c3;
            rmat[[2, 1]] = -c1 * s2;

            rmat[[0, 2]] = s2 * s3;
            rmat[[1, 2]] = s2 * c3;
            rmat[[2, 2]] = c2;
        });

        RMat::new_init(ori)
    }//End of to_rmat

    ///Converts the Bunge angles over to an angle-axis representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    fn to_ang_axis(&self) -> AngAxis{
        let rmat = self.to_rmat();
        rmat.to_ang_axis()

        //This should also work but it doesn't appear to so we're using the above
        //instead
        // let nelems = self.ori.len_of(Axis(1));
        
        // let mut ori = Array2::<f64>::zeros((4, nelems).f());

        // let inv2 = 1.0_f64/2.0_f64;

        // azip!(mut angaxis (ori.axis_iter_mut(Axis(1))), ref bunge (self.ori.axis_iter(Axis(1))) in {
        //     let t       = f64::tan(bunge[1] * inv2);
        //     let sigma   = inv2 * (bunge[0] + bunge[2]);
        //     let csigma  = sigma.cos();
        //     let ssigma  = sigma.sin();
        //     let delta   = inv2 * (bunge[0] - bunge[2]);
        //     let tau     = f64::sqrt(t * t + ssigma * ssigma);
        //     let mut alpha   = 2.0_f64 * f64::atan(tau / csigma);
        //     let itau = 1.0_f64 / tau;
        //     let mut p:f64;

        //     //If alpha is greater than pi we need to set p equal to 1 and
        //     //set alpha equal to 2*pi - alpha. If it isn't then p is set
        //     //to -1. Afterwards everything else is the same.
        //     if alpha > std::f64::consts::PI{
        //         p = 1.0_f64;
        //         alpha = 2.0_f64 * std::f64::consts::PI - alpha;
        //     }else{
        //         p = -1.0_f64;
        //     }

        //     angaxis[0] = p * itau * t * delta.cos();
        //     angaxis[1] = p * itau * t * delta.sin();
        //     angaxis[2] = p * itau * ssigma;
        //     angaxis[3] = alpha;
        // });

        // AngAxis::new_init(ori)
    }//End of to_ang_axis

    ///Converts the Bunge angles over to a compact angle-axis representation which has the following properties
    ///shape (3, nelems), memory order = fortran/column major.
    fn to_ang_axis_comp(&self) -> AngAxisComp{
        //We first convert to a axis-angle representation. Then we scale our normal vector by our the rotation
        //angle which is the fourth component of our axis-angle vector.
        let ang_axis = self.to_ang_axis();
        ang_axis.to_ang_axis_comp()
    }//End of to_ang_axis_comp

    ///Converts the Bunge angles over to a Rodrigues vector representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    fn to_rod_vec(&self) -> RodVec{
        //We first convert to a axis-angle representation. Then we just need to change the last component
        //of our axis-angle representation to be tan(phi/2) instead of phi
        let ang_axis = self.to_ang_axis();
        ang_axis.to_rod_vec()
    }//End of to_rod_vec

    ///Converts the Bunge angles over to a compact Rodrigues vector representation which has the following properties
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

    ///Converts the Bunge angles over to a unit quaternion representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    fn to_quat(&self) -> Quat{

        let ang_axis = self.to_ang_axis();
        ang_axis.to_quat()
        //The below should also work but it appears to be following short so
        //we're going to do the following...
        // let nelems = self.ori.len_of(Axis(1));
        
        // let mut ori = Array2::<f64>::zeros((4, nelems).f());
        
        // //This is a constant factor used in the loop.
        // let inv2 = 1.0_f64/2.0_f64;

        // azip!(mut quat (ori.axis_iter_mut(Axis(1))), ref bunge (self.ori.axis_iter(Axis(1))) in {
        //     let sigma   = inv2 * (bunge[0] + bunge[2]);
        //     let delta   = inv2 * (bunge[0] - bunge[2]);
        //     let c       = f64::cos(inv2 * bunge[1]);
        //     let s       = f64::sin(inv2 * bunge[1]);
        //     let q0      = c * sigma.cos();
        //     let p       = if q0 < 0.0_f64{-1.0_f64} else{1.0_f64};

        //     quat[0] = p * q0;
        //     quat[1] = -p * s * delta.cos();
        //     quat[2] = -p * s * delta.sin();
        //     quat[3] = -p * c * sigma.cos();
        // });

        // Quat::new_init(ori)           
    }//End of to_quat

    ///Converts Bunge angles over to a homochoric representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    fn to_homochoric(&self) -> Homochoric{
        let ang_axis = self.to_ang_axis();
        ang_axis.to_homochoric()
    }//End of to_homochoric    

}//End of Impl Ori_Conv for Bunge