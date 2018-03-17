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

#[derive(Clone, Debug)]
pub struct RMat{
    pub ori: Array3<f64>,
} 


impl RMat{

    //Creates a series of identity matrices for the initial rotation matrix when the data is not fed into it
    pub fn new(size: usize) -> RMat{
        assert!(size > 0, "Size inputted: {}, was not greater than 0", size);
        let mut ori = Array3::<f64>::zeros((3, 3, size).f());
        let eye = Array2::<f64>::eye(3);

        azip!(mut rmat (ori.axis_iter_mut(Axis(2))) in {rmat.assign(&eye);});

        RMat{
            ori,
        }
    }//End of new

    //Creates a rotation matrix type with the supplied data as long as the supplied data is in the following format
    //shape (3, 3, nelems), memory order = fortran/column major.
    //If it doesn't fit those standards it will fail.
    //The data is also being assummed to be orthogonal, but it's currently not checked to see if that is the case.
    //I might add a check later on to make sure that's the case.
    pub fn new_init(ori: Array3<f64>) -> RMat{

        let nrow = ori.len_of(Axis(0));
        let ncol = ori.len_of(Axis(1));

        assert!(nrow == 3 && ncol == 3, "Number of rows and cols of array was: {:?}, which is not equal to 3", [nrow, ncol]);
        //We need to deal with a borrowing of ori here, so we need to have strides dropped at one point.
        {
            let strides = ori.strides();

            assert!(strides[0] == 1, "The memory stride is not column major (f order)");
        }

        RMat{
            ori,
        }
    }//End of new_init

    pub fn to_bunge(&self) -> Bunge{
        
        let nelems = self.ori.len_of(Axis(2));

        let mut ori = Array2::<f64>::zeros((3, nelems));

        //We need to check the R_33 component to see if it's near 1.0 
        let tol = std::f64::EPSILON;

        azip!(mut bunge (ori.axis_iter_mut(Axis(1))), ref rmat (self.ori.axis_iter(Axis(2))) in {
            if f64::abs(rmat[[2, 2]] - 1.0_f64) < tol{
                bunge[0] = f64::atan2(rmat[[0, 1]], rmat[[0, 0]]);
                bunge[1] = std::f64::consts::FRAC_PI_2 * (1.0_f64 - rmat[[2, 2]]);
                bunge[2] = 0.0_f64;
            }else{
                let eta  = 1.0_f64 / f64::sqrt(1.0_f64 - rmat[[2, 2]] * rmat[[2, 2]]);
                bunge[0] = f64::atan2(rmat[[2, 0]] * eta, -rmat[[2, 1]] * eta);
                bunge[1] = rmat[[2, 2]].acos();
                bunge[2] = f64::atan2(rmat[[0, 2]] * eta, rmat[[1, 2]] * eta);
            }
        });

        Bunge{
            ori,
        }
    }//End of to_bunge

    //If to_rmat is called it just returns self
    pub fn to_rmat(&self) -> RMat{
        self.clone()
    }//End of to_rmat

    //Converts the rotation matrix over to an angle-axis representation which has the following properties
    //shape (4, nelems), memory order = fortran/column major.
    pub fn to_ang_axis(&self) -> AngAxis{

        let nelems = self.ori.len_of(Axis(2));

        let mut ori = Array2::<f64>::zeros((4, nelems));

        let inv2 = 1.0_f64/2.0_f64;

        azip!(mut angaxis (ori.axis_iter_mut(Axis(1))), ref rmat (self.ori.axis_iter(Axis(2))) in {
            //The trace of Rmat
            let tr_r = rmat[[0, 0]] + rmat[[1, 1]] + rmat[[2, 2]];
            //This is the angle of rotation about our normal axis.
            let phi  = f64::acos(inv2 * (tr_r - 1.0_f64));
            //The first case is if there is no rotation of axis
            if phi.abs() < std::f64::EPSILON{
                angaxis[0] = 0.0_f64;
                angaxis[1] = 0.0_f64;
                angaxis[2] = 1.0_f64;
                angaxis[3] = 0.0_f64;
            }else{
                let inv_sin = 1.0_f64 / phi.sin();
                //The first three terms are the axial vector of RMat times (1/(2*sin(phi)))
                angaxis[0] = inv_sin * inv2 * (rmat[[2, 1]] - rmat[[1, 2]]);
                angaxis[1] = inv_sin * inv2 * (rmat[[0, 2]] - rmat[[2, 0]]);
                angaxis[2] = inv_sin * inv2 * (rmat[[1, 0]] - rmat[[0, 1]]);
                angaxis[3] = phi;
            }
        });

        AngAxis{
            ori,
        }
    }//End of to_ang_axis

    //Converts the rotation matrix over to a compact angle-axis representation which has the following properties
    //shape (3, nelems), memory order = fortran/column major.
    pub fn to_ang_axis_comp(&self) -> AngAxisComp{
        //We first convert to a angle axis representation. Then we scale our normal vector by our the rotation
        //angle which is the fourth component of our angle axis vector.
        let ang_axis = self.to_ang_axis();
        ang_axis.to_ang_axis_comp()
    }//End of to_ang_axis_comp

    //Converts the rotation matrix over to a rodrigues vector representation which has the following properties
    //shape (4, nelems), memory order = fortran/column major.
    pub fn to_rod_vec(&self) -> RodVec{
        //We first convert to a angle axis representation. Then we just need to change the last component
        //of our angle axis representation to be tan(phi/2) instead of phi
        let ang_axis = self.to_ang_axis();
        ang_axis.to_rod_vec()
    }//End of to_rod_vec

    //Converts the rotation matrix over to a compact rodrigues vector representation which has the following properties
    //shape (3, nelems), memory order = fortran/column major.
    pub fn to_rod_vec_comp(&self) -> RodVecComp{
        //We first convert to a rodrigues vector representation. Then we scale our normal vector by our the rotation
        //angle which is the fourth component of our angle axis vector.
        //If we want to be more efficient about this in the future with out as many copies used we can reuse a lot of the code
        //used in the to_ang_axis code. However, we will end up with a lot of similar/repeated code then. We could put that
        //code in a helper function that isn't seen.
        let rod_vec = self.to_rod_vec();
        rod_vec.to_rod_vec_comp()
    }//End of to_rod_vec_comp

    //Converts the rotation matrix over to a unit quaternion representation which has the following properties
    //shape (4, nelems), memory order = fortran/column major.
    pub fn to_quat(&self) -> Quat{

        let nelems = self.ori.len_of(Axis(2));

        let mut ori = Array2::<f64>::zeros((4, nelems));

        let inv2 = 1.0_f64/2.0_f64;

        azip!(mut quat (ori.axis_iter_mut(Axis(1))), ref rmat (self.ori.axis_iter(Axis(2))) in {
            let q0 = inv2 * f64::sqrt(1.0_f64 + rmat[[0, 0]] + rmat[[1, 1]] + rmat[[2, 2]]);
            let mut q1 = inv2 * f64::sqrt(1.0_f64 + rmat[[0, 0]] - rmat[[1, 1]] - rmat[[2, 2]]);
            let mut q2 = inv2 * f64::sqrt(1.0_f64 - rmat[[0, 0]] + rmat[[1, 1]] - rmat[[2, 2]]);
            let mut q3 = inv2 * f64::sqrt(1.0_f64 - rmat[[0, 0]] - rmat[[1, 1]] + rmat[[2, 2]]);

            if rmat[[1, 2]] > rmat[[2, 1]]{
                q1 *= -1.0_f64;
            }
            if rmat[[2, 0]] > rmat[[0, 2]]{
                q2 *= -1.0_f64;
            }
            if rmat[[0, 1]] > rmat[[1, 0]]{
                q3 *= -1.0_f64;
            }

            //We need to normalize our quaternion when we store it.
            let inv_norm = 1.0_f64 / f64::sqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
            //Once the inv_norm is calculated we apply it to each individual component of our unit
            //quaternion and obtain our 
            quat[0] = q0 * inv_norm;
            quat[1] = q1 * inv_norm;
            quat[2] = q2 * inv_norm;
            quat[3] = q3 * inv_norm;
        });

        Quat{
            ori,
        }
    }//End of to_quat
}//End of Impl RMat
