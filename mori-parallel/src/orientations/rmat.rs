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
use std::cmp;

///A structure that holds an array of rotation/orientation matrices
#[derive(Clone, Debug)]
pub struct RMat{
    ori: Array3<f64>,
} 


impl RMat{

    ///Creates a series of identity matrices for the initial rotation matrix when the data is not fed into it
    pub fn new(size: usize) -> RMat{
        assert!(size > 0, "Size inputted: {}, was not greater than 0", size);
        let mut ori = Array3::<f64>::zeros((3, 3, size).f());
        let eye = Array2::<f64>::eye(3);

        azip!(mut rmat (ori.axis_iter_mut(Axis(2))) in {rmat.assign(&eye);});

        RMat{
            ori,
        }
    }//End of new

    ///Creates a rotation matrix type with the supplied data as long as the supplied data is in the following format
    ///shape (3, 3, nelems), memory order = fortran/column major.
    ///If it doesn't fit those standards it will fail.
    ///The data is also being assummed to be orthogonal, but it's currently not checked to see if that is the case.
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

    ///Return a ndarray view of the orientation data
    pub fn ori_view(&self) -> ArrayView3<f64>{
        self.ori.view()
    }

    ///Return a ndarray mutable view of the orientation data
    pub fn ori_view_mut(&mut self) -> ArrayViewMut3<f64>{
        self.ori.view_mut()
    }

    ///Returns the  transpose of the rotation matrix.
    pub fn par_transpose(&self) -> RMat{
        let nelems = self.ori.len_of(Axis(2));
        let mut ori = Array3::<f64>::zeros((3, 3, nelems).f());
        
        par_azip!(mut rmat_t (ori.axis_iter_mut(Axis(2))), ref rmat (self.ori.axis_iter(Axis(2))) in {
            rmat_t.assign(&rmat.t());
        });

        RMat::new_init(ori)
    }

    //Returns the transpose of the rotation matrix in place
    pub fn par_transpose_inplace(&mut self){
        par_azip!(mut rmat (self.ori.axis_iter_mut(Axis(2))) in {
            rmat.swap([0, 1], [1, 0]);
            rmat.swap([0, 2], [2, 0]);
            rmat.swap([2, 1], [1, 2]);
        });
    }

}//End of RMat impl

///The orientation conversions of a series of rotation matrices to a number of varying different orientation
///representations commonly used in material orientation processing. 
impl ParallelOriConv for RMat{

    ///Converts the rotation matrices to the equivalent bunge angles which has the following properties
    ///shape (3, nelems), memory order = fortran/column major.
    fn par_to_bunge(&self) -> Bunge{
        
        let nelems = self.ori.len_of(Axis(2));

        let mut ori = Array2::<f64>::zeros((3, nelems).f());

        //We need to check the R_33 component to see if it's near 1.0 
        let tol = f64::sqrt(std::f64::EPSILON);

        par_azip!(mut bunge (ori.axis_iter_mut(Axis(1))), ref rmat (self.ori.axis_iter(Axis(2))) in {
            if f64::abs(rmat[[2, 2]]) > (1.0_f64 - tol){
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

        Bunge::new_init(ori)
    }//End of par_to_bunge

    ///Returns a copy of the initial rotation matrix data structure
    fn par_to_rmat(&self) -> RMat{
        self.clone()
    }//End of par_to_rmat

    ///Converts the rotation matrix over to an angle-axis representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    fn par_to_ang_axis(&self) -> AngAxis{

        let nelems = self.ori.len_of(Axis(2));

        let mut ori = Array2::<f64>::zeros((4, nelems).f());

        let inv2 = 1.0_f64/2.0_f64;

        par_azip!(mut angaxis (ori.axis_iter_mut(Axis(1))), ref rmat (self.ori.axis_iter(Axis(2))) in {
            //The trace of Rmat
            let tr_r = rmat[[0, 0]] + rmat[[1, 1]] + rmat[[2, 2]];
            //This is the angle of rotation about our normal axis.
            let phi  = f64::acos(inv2 * (tr_r - 1.0_f64));
            //The first case is if there is no rotation of axis
            if phi.abs() < std::f64::EPSILON{
                angaxis[2] = 1.0_f64;
            }else{
                let inv_sin = 1.0_f64 / phi.sin();
                //The first three terms are the axial vector of RMat times (1/(2*sin(phi)))
                angaxis[0] = inv_sin * inv2 * (rmat[[2, 1]] - rmat[[1, 2]]);
                angaxis[1] = inv_sin * inv2 * (rmat[[0, 2]] - rmat[[2, 0]]);
                angaxis[2] = inv_sin * inv2 * (rmat[[1, 0]] - rmat[[0, 1]]);
                angaxis[3] = phi;
            }
        });

        AngAxis::new_init(ori)
    }//End of par_to_ang_axis

    ///Converts the rotation matrix over to a compact angle-axis representation which has the following properties
    ///shape (3, nelems), memory order = fortran/column major.
    fn par_to_ang_axis_comp(&self) -> AngAxisComp{
        let nelems = self.ori.len_of(Axis(2));

        let mut ori = Array2::<f64>::zeros((3, nelems).f());

        let inv2 = 1.0_f64/2.0_f64;

        par_azip!(mut angaxis (ori.axis_iter_mut(Axis(1))), ref rmat (self.ori.axis_iter(Axis(2))) in {
            //The trace of Rmat
            let tr_r = rmat[[0, 0]] + rmat[[1, 1]] + rmat[[2, 2]];
            //This is the angle of rotation about our normal axis.
            let phi  = f64::acos(inv2 * (tr_r - 1.0_f64));
            if phi.abs() > std::f64::EPSILON{
                let inv_sin = 1.0_f64 / phi.sin();
                //The first three terms are the axial vector of RMat times (1/(2*sin(phi)))
                angaxis[0] = inv_sin * inv2 * (rmat[[2, 1]] - rmat[[1, 2]]) * phi;
                angaxis[1] = inv_sin * inv2 * (rmat[[0, 2]] - rmat[[2, 0]]) * phi;
                angaxis[2] = inv_sin * inv2 * (rmat[[1, 0]] - rmat[[0, 1]]) * phi;
            }
        });

        AngAxisComp::new_init(ori)
    }//End of par_to_ang_axis_comp

    ///Converts the rotation matrix over to a Rodrigues vector representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    fn par_to_rod_vec(&self) -> RodVec{
        let nelems = self.ori.len_of(Axis(2));

        let mut ori = Array2::<f64>::zeros((4, nelems).f());

        let inv2 = 1.0_f64/2.0_f64;

        par_azip!(mut rod_vec (ori.axis_iter_mut(Axis(1))), ref rmat (self.ori.axis_iter(Axis(2))) in {
            //The trace of Rmat
            let tr_r = rmat[[0, 0]] + rmat[[1, 1]] + rmat[[2, 2]];
            //This is the angle of rotation about our normal axis.
            let phi  = f64::acos(inv2 * (tr_r - 1.0_f64));
            //The first case is if there is no rotation of axis
            if phi.abs() < std::f64::EPSILON{
                rod_vec[2] = 1.0_f64;
            }else{
                let inv_sin = 1.0_f64 / phi.sin();
                //The first three terms are the axial vector of RMat times (1/(2*sin(phi)))
                rod_vec[0] = inv_sin * inv2 * (rmat[[2, 1]] - rmat[[1, 2]]);
                rod_vec[1] = inv_sin * inv2 * (rmat[[0, 2]] - rmat[[2, 0]]);
                rod_vec[2] = inv_sin * inv2 * (rmat[[1, 0]] - rmat[[0, 1]]);
                rod_vec[3] = f64::tan(phi * inv2);
            }
        });

        RodVec::new_init(ori)
    }//End of par_to_rod_vec

    ///Converts the rotation matrix over to a compact Rodrigues vector representation which has the following properties
    ///shape (3, nelems), memory order = fortran/column major.
    fn par_to_rod_vec_comp(&self) -> RodVecComp{
        let nelems = self.ori.len_of(Axis(2));

        let mut ori = Array2::<f64>::zeros((4, nelems).f());

        let inv2 = 1.0_f64/2.0_f64;

        par_azip!(mut rod_vec (ori.axis_iter_mut(Axis(1))), ref rmat (self.ori.axis_iter(Axis(2))) in {
            //The trace of Rmat
            let tr_r = rmat[[0, 0]] + rmat[[1, 1]] + rmat[[2, 2]];
            //This is the angle of rotation about our normal axis.
            let phi  = f64::acos(inv2 * (tr_r - 1.0_f64));
            //The first case is if there is no rotation of axis
            if phi.abs() > std::f64::EPSILON{
                let inv_sin = 1.0_f64 / phi.sin();
                let tan2 = f64::tan(phi * inv2);
                //The first three terms are the axial vector of RMat times (1/(2*sin(phi)))
                rod_vec[0] = inv_sin * inv2 * (rmat[[2, 1]] - rmat[[1, 2]]) * tan2;
                rod_vec[1] = inv_sin * inv2 * (rmat[[0, 2]] - rmat[[2, 0]]) * tan2;
                rod_vec[2] = inv_sin * inv2 * (rmat[[1, 0]] - rmat[[0, 1]]) * tan2;
            }
        });

        RodVecComp::new_init(ori)
    }//End of par_to_rod_vec_comp

    ///Converts the rotation matrix over to a unit quaternion representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    fn par_to_quat(&self) -> Quat{
        let nelems = self.ori.len_of(Axis(2));

        let mut ori = Array2::<f64>::zeros((4, nelems).f());

        let inv2 = 1.0_f64/2.0_f64;

        par_azip!(mut quat (ori.axis_iter_mut(Axis(1))), ref rmat (self.ori.axis_iter(Axis(2))) in {
            let mut angaxis = Array1::<f64>::zeros((4).f());
            //The trace of Rmat
            let tr_r = rmat[[0, 0]] + rmat[[1, 1]] + rmat[[2, 2]];
            //This is the angle of rotation about our normal axis.
            let phi  = f64::acos(inv2 * (tr_r - 1.0_f64));
            //The first case is if there is no rotation of axis
            if phi.abs() < std::f64::EPSILON{
                angaxis[2] = 1.0_f64;
            }else{
                let inv_sin = 1.0_f64 / phi.sin();
                //The first three terms are the axial vector of RMat times (1/(2*sin(phi)))
                angaxis[0] = inv_sin * inv2 * (rmat[[2, 1]] - rmat[[1, 2]]);
                angaxis[1] = inv_sin * inv2 * (rmat[[0, 2]] - rmat[[2, 0]]);
                angaxis[2] = inv_sin * inv2 * (rmat[[1, 0]] - rmat[[0, 1]]);
                angaxis[3] = phi;
            }

            let s = f64::sin(inv2 * angaxis[3]); 

            quat[0] = f64::cos(inv2 * angaxis[3]);
            quat[1] = s * angaxis[0];
            quat[2] = s * angaxis[1];
            quat[3] = s * angaxis[2];
        });

        Quat::new_init(ori)
    }//End of par_to_quat

    ///Converts the rotation matrix representation over to a homochoric representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    fn par_to_homochoric(&self) -> Homochoric{
        let ang_axis = self.par_to_ang_axis();
        ang_axis.par_to_homochoric()
    }//End of par_to_homochoric

    ///Converts the rotation matrices to the equivalent bunge angles which has the following properties
    ///shape (3, nelems), memory order = fortran/column major.
    ///This operation is done inplace and does not create a new structure
    fn par_to_bunge_inplace(&self, bunge: &mut Bunge){
        let mut ori = bunge.ori_view_mut();

        let new_nelem = ori.len_of(Axis(1));
        let nelem = self.ori.len_of(Axis(2));

        assert!(new_nelem == nelem, 
        "The number of elements in the original ori field do no match up with the new field.
        The old field had {} elements, and the new field has {} elements",
        nelem, new_nelem);

        //We need to check the R_33 component to see if it's near 1.0 
        let tol = f64::sqrt(std::f64::EPSILON);

        par_azip!(mut bunge (ori.axis_iter_mut(Axis(1))), ref rmat (self.ori.axis_iter(Axis(2))) in {
            if f64::abs(rmat[[2, 2]]) > (1.0_f64 - tol){
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

    }

    ///Returns a copy of the initial rotation matrix data structure
    ///This operation is done inplace and does not create a new structure
    fn par_to_rmat_inplace(&self, rmat: &mut RMat){
        let mut ori = rmat.ori_view_mut();

        let new_nelem = ori.len_of(Axis(2));
        let nelem = self.ori.len_of(Axis(2));

        assert!(new_nelem == nelem, 
        "The number of elements in the original ori field do no match up with the new field.
        The old field had {} elements, and the new field has {} elements",
        nelem, new_nelem);

        ori.assign(&self.ori);
    }

    ///Converts the rotation matrix over to an angle-axis representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    ///This operation is done inplace and does not create a new structure
    fn par_to_ang_axis_inplace(&self, ang_axis: &mut AngAxis){

        let mut ori = ang_axis.ori_view_mut();

        let new_nelem = ori.len_of(Axis(1));
        let nelem = self.ori.len_of(Axis(2));

        assert!(new_nelem == nelem, 
        "The number of elements in the original ori field do no match up with the new field.
        The old field had {} elements, and the new field has {} elements",
        nelem, new_nelem);

        let inv2 = 1.0_f64/2.0_f64;

        par_azip!(mut angaxis (ori.axis_iter_mut(Axis(1))), ref rmat (self.ori.axis_iter(Axis(2))) in {
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
    }

    ///Converts the rotation matrix over to a compact angle-axis representation which has the following properties
    ///shape (3, nelems), memory order = fortran/column major.
    ///This operation is done inplace and does not create a new structure
    fn par_to_ang_axis_comp_inplace(&self, ang_axis_comp: &mut AngAxisComp){

        let mut ori = ang_axis_comp.ori_view_mut();

        let new_nelem = ori.len_of(Axis(1));
        let nelem = self.ori.len_of(Axis(2));

        assert!(new_nelem == nelem, 
        "The number of elements in the original ori field do no match up with the new field.
        The old field had {} elements, and the new field has {} elements",
        nelem, new_nelem);

        let inv2 = 1.0_f64/2.0_f64;

        par_azip!(mut angaxis (ori.axis_iter_mut(Axis(1))), ref rmat (self.ori.axis_iter(Axis(2))) in {
            //The trace of Rmat
            let tr_r = rmat[[0, 0]] + rmat[[1, 1]] + rmat[[2, 2]];
            //This is the angle of rotation about our normal axis.
            let phi  = f64::acos(inv2 * (tr_r - 1.0_f64));
            if phi.abs() > std::f64::EPSILON{
                let inv_sin = 1.0_f64 / phi.sin();
                //The first three terms are the axial vector of RMat times (1/(2*sin(phi)))
                angaxis[0] = inv_sin * inv2 * (rmat[[2, 1]] - rmat[[1, 2]]) * phi;
                angaxis[1] = inv_sin * inv2 * (rmat[[0, 2]] - rmat[[2, 0]]) * phi;
                angaxis[2] = inv_sin * inv2 * (rmat[[1, 0]] - rmat[[0, 1]]) * phi;
            }
        });
    }

    ///Converts the rotation matrix over to a Rodrigues vector representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    ///This operation is done inplace and does not create a new structure
    fn par_to_rod_vec_inplace(&self, rod_vec: &mut RodVec){

        let mut ori = rod_vec.ori_view_mut();

        let new_nelem = ori.len_of(Axis(1));
        let nelem = self.ori.len_of(Axis(2));

        assert!(new_nelem == nelem, 
        "The number of elements in the original ori field do no match up with the new field.
        The old field had {} elements, and the new field has {} elements",
        nelem, new_nelem);

        let inv2 = 1.0_f64/2.0_f64;

        par_azip!(mut rod_vec (ori.axis_iter_mut(Axis(1))), ref rmat (self.ori.axis_iter(Axis(2))) in {
            //The trace of Rmat
            let tr_r = rmat[[0, 0]] + rmat[[1, 1]] + rmat[[2, 2]];
            //This is the angle of rotation about our normal axis.
            let phi  = f64::acos(inv2 * (tr_r - 1.0_f64));
            //The first case is if there is no rotation of axis
            if phi.abs() < std::f64::EPSILON{
                rod_vec[2] = 1.0_f64;
            }else{
                let inv_sin = 1.0_f64 / phi.sin();
                //The first three terms are the axial vector of RMat times (1/(2*sin(phi)))
                rod_vec[0] = inv_sin * inv2 * (rmat[[2, 1]] - rmat[[1, 2]]);
                rod_vec[1] = inv_sin * inv2 * (rmat[[0, 2]] - rmat[[2, 0]]);
                rod_vec[2] = inv_sin * inv2 * (rmat[[1, 0]] - rmat[[0, 1]]);
                rod_vec[3] = f64::tan(phi * inv2);
            }
        });
    }

    ///Converts the rotation matrix over to a compact Rodrigues vector representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    ///This operation is done inplace and does not create a new structure
    fn par_to_rod_vec_comp_inplace(&self, rod_vec_comp: &mut RodVecComp){

        let mut ori = rod_vec_comp.ori_view_mut();

        let new_nelem = ori.len_of(Axis(1));
        let nelem = self.ori.len_of(Axis(2));

        assert!(new_nelem == nelem, 
        "The number of elements in the original ori field do no match up with the new field.
        The old field had {} elements, and the new field has {} elements",
        nelem, new_nelem);

        let inv2 = 1.0_f64/2.0_f64;

        par_azip!(mut rod_vec (ori.axis_iter_mut(Axis(1))), ref rmat (self.ori.axis_iter(Axis(2))) in {
            //The trace of Rmat
            let tr_r = rmat[[0, 0]] + rmat[[1, 1]] + rmat[[2, 2]];
            //This is the angle of rotation about our normal axis.
            let phi  = f64::acos(inv2 * (tr_r - 1.0_f64));
            //The first case is if there is no rotation of axis
            if phi.abs() > std::f64::EPSILON{
                let inv_sin = 1.0_f64 / phi.sin();
                let tan2 = f64::tan(phi * inv2);
                //The first three terms are the axial vector of RMat times (1/(2*sin(phi)))
                rod_vec[0] = inv_sin * inv2 * (rmat[[2, 1]] - rmat[[1, 2]]) * tan2;
                rod_vec[1] = inv_sin * inv2 * (rmat[[0, 2]] - rmat[[2, 0]]) * tan2;
                rod_vec[2] = inv_sin * inv2 * (rmat[[1, 0]] - rmat[[0, 1]]) * tan2;
            }
        });
    }
    ///This operation is done inplace and does not create a new structure
    fn par_to_quat_inplace(&self, quat: &mut Quat){

        let mut ori = quat.ori_view_mut();

        let new_nelem = ori.len_of(Axis(1));
        let nelem = self.ori.len_of(Axis(2));

        assert!(new_nelem == nelem, 
        "The number of elements in the original ori field do no match up with the new field.
        The old field had {} elements, and the new field has {} elements",
        nelem, new_nelem);

        let inv2 = 1.0_f64/2.0_f64;

        par_azip!(mut quat (ori.axis_iter_mut(Axis(1))), ref rmat (self.ori.axis_iter(Axis(2))) in {
            let mut angaxis = Array1::<f64>::zeros((4).f());
            //The trace of Rmat
            let tr_r = rmat[[0, 0]] + rmat[[1, 1]] + rmat[[2, 2]];
            //This is the angle of rotation about our normal axis.
            let phi  = f64::acos(inv2 * (tr_r - 1.0_f64));
            //The first case is if there is no rotation of axis
            if phi.abs() < std::f64::EPSILON{
                angaxis[2] = 1.0_f64;
            }else{
                let inv_sin = 1.0_f64 / phi.sin();
                //The first three terms are the axial vector of RMat times (1/(2*sin(phi)))
                angaxis[0] = inv_sin * inv2 * (rmat[[2, 1]] - rmat[[1, 2]]);
                angaxis[1] = inv_sin * inv2 * (rmat[[0, 2]] - rmat[[2, 0]]);
                angaxis[2] = inv_sin * inv2 * (rmat[[1, 0]] - rmat[[0, 1]]);
                angaxis[3] = phi;
            }

            let s = f64::sin(inv2 * angaxis[3]); 

            quat[0] = f64::cos(inv2 * angaxis[3]);
            quat[1] = s * angaxis[0];
            quat[2] = s * angaxis[1];
            quat[3] = s * angaxis[2];
        });
    }
    ///This operation is done inplace and does not create a new structure
    fn par_to_homochoric_inplace(&self, homochoric: &mut Homochoric){
        let ang_axis = self.par_to_ang_axis();
        ang_axis.par_to_homochoric_inplace(homochoric);
    }



}//End of Impl Ori_Conv for RMat

///A series of commonly used operations to rotate vector data by a given rotation
impl ParallelRotVector for RMat{

    ///rot_vector takes in a 2D array view of a series of vectors. It then rotates these vectors using the
    ///given rotation matrices. The newly rotated vectors are then returned. This function requires the
    ///number of elements in the rotation matrix to be either 1 or nelems.
    ///It also requires that the number of elements in the unrotated vector be 1 or nelems.
    ///If this condition is not met the function will error out.
    ///vec - the vector to be rotated must have dimensions 3xnelems or 3x1.
    ///Output - the rotated vector and has dimensions 3xnelems.
    fn par_rot_vector(&self, vec: ArrayView2<f64>) -> Array2<f64>{

        let nelems = vec.len_of(Axis(1));
        let rnelems = self.ori.len_of(Axis(2));

        let rows  = vec.len_of(Axis(0));
        assert!((rows == 3), "The number of rows must be 3. The number of rows provided is {}", rows); 

        assert!( (nelems == rnelems) | (rnelems == 1) | (nelems == 1), 
        "The number of elements in the vector field must be equal to the number of elements in the
        Rotation Matix structure, or their must only be one element in Rotation Matrix. The final case is
        that there must only be one element in the vector field. There are
        currently {} elements in vector and {} elements in Rotation Matrix",
        nelems, rnelems);

        let mnelems = cmp::max(rnelems, nelems);
        let mut rvec = Array2::<f64>::zeros((3, mnelems).f());

        //We need to see if we have more than one rotation matrix that we're multiplying by
        if rnelems == nelems {
            //Here we're iterating through each vector, rotation matrix, and rotated vector value
            //and assigning R*v to the rotated vector.
            par_azip!(mut rvec (rvec.axis_iter_mut(Axis(1))), ref vec (vec.axis_iter(Axis(1))), 
            ref rmat (self.ori.axis_iter(Axis(2))) in {
                rvec.assign({&rmat.dot(&vec)});
            });
        }else if rnelems == 1{
            //We just have one rmat so we can multiple that by our vec and get all of our
            //rvec values all at once.
            //The subview is necessary because we represent RMat as an Array3 and we need an Array1 or Array2
            //to use the dot function
            rvec.assign({&self.ori.subview(Axis(2), 0).dot(&vec)});
        }else{
            //We now need to look at the case where we only have one vector to rotate but several rotation matrices
            let vec = vec.subview(Axis(1), 0);
            par_azip!(mut rvec (rvec.axis_iter_mut(Axis(1))), ref rmat (self.ori.axis_iter(Axis(2))) in {
                rvec.assign({&rmat.dot(&vec)});
            });
        }//End if-else
        //Now we just need to return the rvec value
        rvec
    }//End of rot_vector

    ///rot_vector_mut takes in a 2D array view of a series of vectors and a mutable 2D ArrayView of the 
    ///rotated vector. It then rotates these vectors using the given rotation matrices. The newly rotated
    /// vectors are assigned to the supplied rotated vector, rvec. This function requires the
    ///number of elements in the rotation matrix to be either 1 or nelems where rvec has nelems in it.
    ///It also requires that the number of elements in the unrotated vector be 1 or nelems.
    ///If these conditions are not met the function will error out.
    ///vec - the vector to be rotated must have dimensions 3xnelems or 3x1.
    ///rvec - the rotated vector and has dimensions 3xnelems.
    fn par_rot_vector_mut(&self, vec: ArrayView2<f64>, mut rvec: ArrayViewMut2<f64>) {

        let nelems = vec.len_of(Axis(1));
        let rvnelems = rvec.len_of(Axis(1));
        let rnelems = self.ori.len_of(Axis(2));
        let mnelems = cmp::max(rnelems, nelems);

        let rows  = vec.len_of(Axis(0));
        assert!((rows == 3), "The number of rows must be 3. The number of rows provided is {}", rows); 

        assert!((mnelems == rvnelems),
        "The number of elements in the unrotated vector/rotation matrix field must be equal to the number of elements
        in the supplied rotated vector field. There are currently {} elements in the unrotated vector
        field/rotated matrix field and {} elements in the rotated vector field", 
        mnelems, rvnelems);

        assert!( (nelems == rnelems) | (rnelems == 1) | (nelems == 1), 
        "The number of elements in the vector field must be equal to the number of elements in the
        Rotation Matix structure, or their must only be one element in Rotation Matrix. The final case is
        that there must only be one element in the vector field. There are
        currently {} elements in vector and {} elements in Rotation Matrix",
        nelems, rnelems);



        //We need to see if we have more than one rotation matrix that we're multiplying by
        if rnelems == nelems {
            //Here we're iterating through each vector, rotation matrix, and rotated vector value
            //and assigning R*v to the rotated vector.
            par_azip!(mut rvec (rvec.axis_iter_mut(Axis(1))), ref vec (vec.axis_iter(Axis(1))), 
            ref rmat (self.ori.axis_iter(Axis(2))) in {
                rvec.assign({&rmat.dot(&vec)});
            });
        } else if rnelems == 1{
            //We just have one rmat so we can multiple that by our vec and get all of our
            //rvec values all at once.
            //The subview is necessary because we represent RMat as an Array3 and we need an Array1 or Array2
            //to use the dot function
            rvec.assign({&self.ori.subview(Axis(2), 0).dot(&vec)});
        }else{
            //We now need to look at the case where we only have one vector to rotate but several rotation matrices
            let vec = vec.subview(Axis(1), 0);
            par_azip!(mut rvec (rvec.axis_iter_mut(Axis(1))), ref rmat (self.ori.axis_iter(Axis(2))) in {
                rvec.assign({&rmat.dot(&vec)});
            });
        }//End if-else
    }//End of rot_vector_mut

    ///rot_vector_inplace takes in a mutable 2D array view of a series of vectors. It then rotates these vectors using the
    ///given rotation matrices. The newly rotated vectors are assigned to original vector. This function requires the
    ///number of elements in the rotation matrix to be either 1 or nelems where vec has nelems in it.
    ///If this condition is not met the function will error out.
    ///vec - the vector to be rotated must have dimensions 3xnelems.
    fn par_rot_vector_inplace(&self, mut vec: ArrayViewMut2<f64>){

        let nelems = vec.len_of(Axis(1));
        let rnelems = self.ori.len_of(Axis(2));

        let rows  = vec.len_of(Axis(0));
        assert!((rows == 3), "The number of rows must be 3. The number of rows provided is {}", rows); 

        assert!( (nelems == rnelems) | (rnelems == 1), 
        "The number of elements in the vector field must be equal to the number of elements in the
        Rotation Matix structure, or their must only be one element in Rotation Matrix. There are
        currently {} elements in vector and {} elements in Rotation Matrix",
        nelems, rnelems);

        //We need to see if we have more than one rotation matrix that we're multiplying by
        if rnelems == nelems {
            //Here we're iterating through each vector, rotation matrix, and rotated vector value
            //and assigning R*v to the vector.
            par_azip!(mut vec (vec.axis_iter_mut(Axis(1))), ref rmat (self.ori.axis_iter(Axis(2))) in {
                //A cleaner way needs to exists to perform this operation.
                let mut rvec = Array1::<f64>::zeros((3).f());
                rvec.assign({&rmat.dot(&vec)});
                vec.assign({&rvec});
            });
        } else{
            //We just have one rmat so we can multiple that by our vec and get all of our
            //rvec values all at once.
            //The subview is necessary because we represent RMat as an Array3 and we need an Array1 or Array2
            //to use the dot function
            par_azip!(mut vec (vec.axis_iter_mut(Axis(1))) in{
                //A cleaner way needs to exists to perform this operation.
                let mut rvec = Array1::<f64>::zeros((3).f());
                rvec.assign({&self.ori.subview(Axis(2), 0).dot(&vec)});
                vec.assign({&rvec});
            });
        }//End if-else
    }//End of rot_vector_inplace
}//Endo of Impl RotVector

///A series of commonly used operations to rotate 2nd order 3x3 tensor data by a given rotation
impl ParallelRotTensor for RMat{

    ///rot_tensor takes in a 3D array view of a series of tensors. It then rotates these tensors using the
    ///given rotation matrices. The newly rotated tensors are then returned. This function requires the
    ///number of elements in the rotation matrix to be either 1 or nelems where tensor has nelems in it.
    ///If this condition is not met the function will error out.
    ///tensor - the tensors to be rotated must have dimensions 3x3xnelems.
    ///Output - the rotated tensors and has dimensions 3x3xnelems.
    fn par_rot_tensor(&self, tensor: ArrayView3<f64>) -> Array3<f64>{

        let nelems  = tensor.len_of(Axis(2));
        let rnelems = self.ori.len_of(Axis(2));
        let rows    = tensor.len_of(Axis(0));
        let cols    = tensor.len_of(Axis(1));

        assert!((cols == 3) & (rows == 3), 
        "The number of columns and rows must be 3. The number of rows and cols provided is {}
        and {} respectively", rows, cols); 

        assert!( (nelems == rnelems) | (rnelems == 1), 
        "The number of elements in the tensor field must be equal to the number of elements in the
        Rotation Matix structure, or their must only be one element in Rotation Matrix. There are
        currently {} elements in tensor and {} elements in Rotation Matrix",
        nelems, rnelems);

        let mut rtensor = Array3::<f64>::zeros((3, 3, nelems).f());

        //We need to see if we have more than one rotation matrix that we're multiplying by
        if rnelems == nelems {
            //Here we're iterating through each tensor, rotation matrix, and rotated tensor value
            //and assigning R*T*R^T to the rotated tensor.
            par_azip!(mut rtensor (rtensor.axis_iter_mut(Axis(2))), ref tensor (tensor.axis_iter(Axis(2))), 
            ref rmat (self.ori.axis_iter(Axis(2))) in {
                rtensor.assign({&rmat.dot(&tensor.dot(&rmat.t()))});
            });
        } else{
            //We just have one rmat so we can multiple that by our tensors.
            //The subview is necessary because we represent RMat as an Array3 and we need an Array1 or Array2
            //to use the dot function
            par_azip!(mut rtensor (rtensor.axis_iter_mut(Axis(2))), ref tensor (tensor.axis_iter(Axis(2))) in {
            rtensor.assign({
                &self.ori.subview(Axis(2), 0).dot(&tensor.dot(&self.ori.subview(Axis(2), 0).t()))
                });
            });
        }

        //Now we just need to return the rtensor value
        rtensor
    }//End of rot_tensor

    ///rot_tensor_mut takes in a 3D array view of a series of tensors and a mutable 3D ArrayView of the 
    ///rotated tensors. It then rotates these tensors using the given rotation matrices. The newly rotated
    ///tensors are assigned to the supplied rotated tensors, rtensor. This function requires the
    ///number of elements in the rotation matrix to be either 1 or nelems where tensor has nelems in it.
    ///It also requires the number of elements in rtensor and tensor to be equal.
    ///tensor  - the tensors to be rotated must have dimensions 3x3xnelems.
    ///rtensor - the rotated tensors and has dimensions 3x3xnelems.
    fn par_rot_tensor_mut(&self, tensor: ArrayView3<f64>, mut rtensor: ArrayViewMut3<f64>) {

        let nelems  = tensor.len_of(Axis(2));
        let rtnelems  = rtensor.len_of(Axis(2));
        let rnelems = self.ori.len_of(Axis(2));
        let rows    = tensor.len_of(Axis(0));
        let cols    = tensor.len_of(Axis(1));

        assert!((cols == 3) & (rows == 3), 
        "The number of columns and rows must be 3. The number of rows and cols provided is {}
        and {} respectively", rows, cols); 

        assert!((nelems == rtnelems),
        "The number of elements in the unrotated tensor field must be equal to the number of elements
        in the supplied rotated tensor field. There are currently {} elements in the unrotated tensor
        field and {} elements in the rotated tensor field", 
        nelems, rtnelems);

        assert!( (nelems == rnelems) | (rnelems == 1), 
        "The number of elements in the tensor field must be equal to the number of elements in the
        Rotation Matix structure, or their must only be one element in Rotation Matrix. There are
        currently {} elements in tensor and {} elements in Rotation Matrix",
        nelems, rnelems);

        //We need to see if we have more than one rotation matrix that we're multiplying by
        if rnelems == nelems {
            //Here we're iterating through each tensor, rotation matrix, and rotated tensor value
            //and assigning R*T*R^T to the rotated tensor.
            par_azip!(mut rtensor (rtensor.axis_iter_mut(Axis(2))), ref tensor (tensor.axis_iter(Axis(2))), 
            ref rmat (self.ori.axis_iter(Axis(2))) in {
                rtensor.assign({&rmat.dot(&tensor.dot(&rmat.t()))});
            });
        } else{
            //We just have one rmat so we can multiple that by our tensors.
            //The subview is necessary because we represent RMat as an Array3 and we need an Array1 or Array2
            //to use the dot function
            par_azip!(mut rtensor (rtensor.axis_iter_mut(Axis(2))), ref tensor (tensor.axis_iter(Axis(2))) in {
                rtensor.assign({
                    &self.ori.subview(Axis(2), 0).dot(&tensor.dot(&self.ori.subview(Axis(2), 0).t()))
                    });
            });
        }
    }//End of rot_tensor_mut

    ///rot_tensor_inplace takes in a mutable 3D array view of a series of tensors. 
    ///It then rotates these tensors using the given rotation matrices. The newly rotated
    ///tensors are assigned in place of the supplied variable tensor. This function requires the
    ///number of elements in the rotation matrix to be either 1 or nelems where tensor has nelems in it.
    ///tensor - the tensors to be rotated must have dimensions 3x3xnelems.
    fn par_rot_tensor_inplace(&self, mut tensor: ArrayViewMut3<f64>) {

        let nelems  = tensor.len_of(Axis(2));
        let rnelems = self.ori.len_of(Axis(2));
        let rows    = tensor.len_of(Axis(0));
        let cols    = tensor.len_of(Axis(1));

        assert!((cols == 3) & (rows == 3), 
        "The number of columns and rows must be 3. The number of rows and cols provided is {}
        and {} respectively", rows, cols); 

        assert!( (nelems == rnelems) | (rnelems == 1), 
        "The number of elements in the tensor field must be equal to the number of elements in the
        Rotation Matix structure, or their must only be one element in Rotation Matrix. There are
        currently {} elements in tensor and {} elements in Rotation Matrix",
        nelems, rnelems);

        //We need to see if we have more than one rotation matrix that we're multiplying by
        if rnelems == nelems {
            //Here we're iterating through each tensor, rotation matrix, and rotated tensor value
            //and assigning R*T*R^T to the rotated tensor.
            par_azip!(mut tensor (tensor.axis_iter_mut(Axis(2))), ref rmat (self.ori.axis_iter(Axis(2))) in {
                //A cleaner way needs to exists to perform this operation.
                let mut rtensor = Array2::<f64>::zeros((3, 3).f());
                rtensor.assign({&rmat.dot(&tensor.dot(&rmat.t()))});
                tensor.assign({&rtensor});
            });
        } else{
            //We just have one rmat so we can multiple that by our tensors.
            //The subview is necessary because we represent RMat as an Array3 and we need an Array1 or Array2
            //to use the dot function
            par_azip!(mut  tensor (tensor.axis_iter_mut(Axis(2))) in {
                //A cleaner way needs to exists to perform this operation.
                let mut rtensor = Array2::<f64>::zeros((3, 3).f());
                rtensor.assign({&self.ori.subview(Axis(2), 0).dot(&tensor.dot(&self.ori.subview(Axis(2), 0).t()))});
                tensor.assign({&rtensor});
            });
        }
    }//End of rot_tensor_inplace

}//End of Impl RotTensor for RMat
