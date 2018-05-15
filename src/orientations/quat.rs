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

///A structure that holds an array of unit quaternions
#[derive(Clone, Debug)]
pub struct Quat{
    ori: Array2<f64>,
}

impl Quat{
    ///Creates an array of zeros for the initial unit quaternion when data is not fed into it
    pub fn new(size: usize) -> Quat{
        assert!(size > 0, "Size inputted: {}, was not greater than 0", size);

        let mut ori = Array2::<f64>::zeros((4, size).f());

        azip!(mut quat (ori.axis_iter_mut(Axis(1))) in {quat[0] = 1.0_f64});

        Quat{
            ori,
        }
    }//End of new

    ///Creates a unit quaternion type with the supplied data as long as the supplied data is in the following format
    ///shape (4, nelems), memory order = fortran/column major.
    ///If it doesn't fit those standards it will fail.
    pub fn new_init(ori: Array2<f64>) -> Quat{

        let nrow = ori.rows();

        assert!(nrow == 4, "Number of rows of array was: {}, which is not equal to 4", nrow);
        //We need to deal with a borrowing of ori here, so we need to have strides dropped at one point.
        {
            let strides = ori.strides();

            assert!(strides[0] == 1, "The memory stride is not column major (f order)");
        }

        Quat{
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

    ///Returns a new Quat that is equal to the conjugate/inverse of the unit quaternion which is simply the negative
    ///of the vector portions of the unit quaternion.
    pub fn conjugate(&self) -> Quat{
        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((4, nelems).f());
        
        azip!(mut quat_c (ori.axis_iter_mut(Axis(1))), ref quat (self.ori.axis_iter(Axis(1))) in {
            quat_c[0] = quat[0];
            quat_c[1] = -1.0_f64 * quat[1];
            quat_c[2] = -1.0_f64 * quat[2];
            quat_c[3] = -1.0_f64 * quat[3];
        });

        Quat::new_init(ori)
    }

    ///Performs in place the conjugate/inverse of the unit quaternion which is simply the negative
    ///of the vector portions of the unit quaternion.
    pub fn conjugate_inplace(&mut self){
        azip!(mut quat_c (self.ori.axis_iter_mut(Axis(1))) in {
            quat_c[1] *= -1.0_f64;
            quat_c[2] *= -1.0_f64;
            quat_c[3] *= -1.0_f64;
        });
    }
}//End of Impl of Quat

///The orientation conversions of a series of unit quaternions to a number of varying different orientation
///representations commonly used in material orientation processing. 
impl OriConv for Quat{
    ///Converts the unit quaternion representation over to Bunge angles which has the following properties
    ///shape (3, nelems), memory order = fortran/column major.
    fn to_bunge(&self) -> Bunge{

        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((3, nelems).f());

        let tol = f64::sqrt(std::f64::EPSILON);

        azip!(mut bunge (ori.axis_iter_mut(Axis(1))), ref quat (self.ori.axis_iter(Axis(1))) in {
            let q03 = quat[0] * quat[0] + quat[3] * quat[3];
            let q12 = quat[1] * quat[1] + quat[2] * quat[2];
            let xi = f64::sqrt(q03 * q12);
            //We get to now go through all of the different cases that this might break down into
            if xi.abs() < tol && q12.abs() < tol {
                bunge[0] = f64::atan2(-2.0_f64 * quat[0] * quat[3], quat[0] * quat[0] - quat[3] * quat[3]);
                //All of the other values are zero
            }else if xi.abs() < tol && q03.abs() < tol{
                bunge[0] = f64::atan2(2.0_f64 * quat[1] * quat[2], quat[1] * quat[1] - quat[2] * quat[2]);
                bunge[1] = std::f64::consts::PI;
                //The other value is zero
            }else{
                let inv_xi = 1.0_f64 / xi;
                //The atan2 terms are pretty long so we're breaking it down into a couple of temp terms
                let t1 = inv_xi * (quat[1] * quat[3] - quat[0] * quat[2]);
                let t2 = inv_xi * (-quat[0] * quat[1] - quat[2] * quat[3]);
                //We can now assign the first two bunge angles
                bunge[0] = t1.atan2(t2);
                bunge[1] = f64::atan2(2.0_f64 * xi, q03 - q12);
                //Once again these terms going into the atan2 term are pretty long
                let t1 = inv_xi * (quat[0] * quat[2] + quat[1] * quat[3]);
                let t2 = inv_xi * (quat[2] * quat[3] - quat[0] * quat[1]);
                //We can finally find the final bunge angle
                bunge[2] = t1.atan2(t2);
            }
        });

        Bunge::new_init(ori)
    }//End of to_bunge

    ///Converts the unit quaternion representation over to rotation matrix which has the following properties
    ///shape (3, 3, nelems), memory order = fortran/column major.
    fn to_rmat(&self) -> RMat{

        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array3::<f64>::zeros((3, 3, nelems).f());

        azip!(mut rmat (ori.axis_iter_mut(Axis(2))), ref quat (self.ori.axis_iter(Axis(1))) in {
            let qbar =  quat[0] * quat[0] - (quat[1] * quat[1] + quat[2] * quat[2] + quat[3] * quat[3]);

            rmat[[0, 0]] = qbar + 2.0_f64 * quat[1] * quat[1];
            rmat[[1, 0]] = 2.0_f64 * (quat[1] * quat[2] + quat[0] * quat[3]);
            rmat[[2, 0]] = 2.0_f64 * (quat[1] * quat[3] - quat[0] * quat[2]);

            rmat[[0, 1]] = 2.0_f64 * (quat[1] * quat[2] - quat[0] * quat[3]);
            rmat[[1, 1]] = qbar + 2.0_f64 * quat[2] * quat[2];
            rmat[[2, 1]] = 2.0_f64 * (quat[2] * quat[3] + quat[0] * quat[1]);

            rmat[[0, 2]] = 2.0_f64 * (quat[1] * quat[3] + quat[0] * quat[2]);
            rmat[[1, 2]] = 2.0_f64 * (quat[2] * quat[3] - quat[0] * quat[1]);
            rmat[[2, 2]] = qbar + 2.0_f64 * quat[3] * quat[3];
        });

        RMat::new_init(ori)
    }//End of to_rmat

    ///Converts the unit quaternion representation over to angle-axis representation which has the following properties
    ///shape (3, nelems), memory order = fortran/column major.
    fn to_ang_axis(&self) -> AngAxis{

        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((4, nelems).f());

        let tol = std::f64::EPSILON;

        azip!(mut angaxis (ori.axis_iter_mut(Axis(1))), ref quat (self.ori.axis_iter(Axis(1))) in {
            let phi = 2.0_f64 * quat[0].acos();
            if quat[0].abs() < tol{
                angaxis[0] = quat[1];
                angaxis[1] = quat[2];
                angaxis[2] = quat[3];
                angaxis[3] = std::f64::consts::PI;
            }else if phi.abs() < tol{
                angaxis[2] = 1.0_f64;
            }else{
                let s   = quat[0].signum() / f64::sqrt(quat[1] * quat[1] + quat[2] * quat[2] + quat[3] * quat[3]);

                angaxis[0] = s * quat[1];
                angaxis[1] = s * quat[2];
                angaxis[2] = s * quat[3];
                angaxis[3] = phi;
            }
        });

        AngAxis::new_init(ori)
    }//End of to_ang_axis

    ///Converts the unit quaternion over to a compact angle-axis representation which has the following properties
    ///shape (3, nelems), memory order = fortran/column major.
    fn to_ang_axis_comp(&self) -> AngAxisComp{

        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((3, nelems).f());

        let tol = std::f64::EPSILON;

        azip!(mut angaxis (ori.axis_iter_mut(Axis(1))), ref quat (self.ori.axis_iter(Axis(1))) in {
            let phi = 2.0_f64 * quat[0].acos();
            if quat[0].abs() < tol{
                angaxis[0] = quat[1] * std::f64::consts::PI;
                angaxis[1] = quat[2] * std::f64::consts::PI;
                angaxis[2] = quat[3] * std::f64::consts::PI;
            }else{
                let s   = quat[0].signum() / f64::sqrt(quat[1] * quat[1] + quat[2] * quat[2] + quat[3] * quat[3]);

                angaxis[0] = s * quat[1] * phi;
                angaxis[1] = s * quat[2] * phi;
                angaxis[2] = s * quat[3] * phi; 
            }
        });

        AngAxisComp::new_init(ori)
    }//End of to_ang_axis_comp

    ///Converts the unit quaternion over to a Rodrigues vector representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    fn to_rod_vec(&self) -> RodVec{

        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((4, nelems).f());

        let tol = std::f64::EPSILON;

        azip!(mut rod_vec (ori.axis_iter_mut(Axis(1))), ref quat (self.ori.axis_iter(Axis(1))) in {
            let phi = quat[0].acos();
            if quat[0].abs() < tol{
                rod_vec[0] = quat[1];
                rod_vec[1] = quat[2];
                rod_vec[2] = quat[3];
                rod_vec[3] = std::f64::INFINITY;
            }else if phi.abs() < tol{
                rod_vec[2] = 1.0_f64;
            }else{
                let s   = quat[0].signum() / f64::sqrt(quat[1] * quat[1] + quat[2] * quat[2] + quat[3] * quat[3]);

                rod_vec[0] = s * quat[1];
                rod_vec[1] = s * quat[2];
                rod_vec[2] = s * quat[3];
                rod_vec[3] = phi.tan();
            }
        });

        RodVec::new_init(ori)
    }//End of to_rod_vec

    ///Converts the unit quaternion over to a compact Rodrigues vector representation which has the following properties
    ///shape (3, nelems), memory order = fortran/column major.
    fn to_rod_vec_comp(&self) -> RodVecComp{
        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((3, nelems).f());
        let tol = std::f64::EPSILON;

        azip!(mut rod_vec_comp (ori.axis_iter_mut(Axis(1))), ref quat (self.ori.axis_iter(Axis(1))) in {
            let tan_phi = f64::tan(quat[0].acos());
            //This case will not allow for anything to be retrievable later on...
            if quat[0].abs() < tol{
                rod_vec_comp[0] = std::f64::INFINITY;
                rod_vec_comp[1] = std::f64::INFINITY;
                rod_vec_comp[2] = std::f64::INFINITY;
            }else{
                let s   = quat[0].signum() / f64::sqrt(quat[1] * quat[1] + quat[2] * quat[2] + quat[3] * quat[3]);

                rod_vec_comp[0] = s * quat[1] * tan_phi;
                rod_vec_comp[1] = s * quat[2] * tan_phi;
                rod_vec_comp[2] = s * quat[3] * tan_phi; 
            }
        });

        RodVecComp::new_init(ori)
    }//End of to_rod_vec_comp

    ///This returns a clone of the original unit quaternion structure
    fn to_quat(&self) -> Quat{
        self.clone()
    }//End of to_quat

    ///Converts the quaternion representation over to a homochoric representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    fn to_homochoric(&self) -> Homochoric{
        let ang_axis = self.to_ang_axis();
        ang_axis.to_homochoric()
    }//End of to_homochoric

    ///Converts the unit quaternion representation over to Bunge angles which has the following properties
    ///shape (3, nelems), memory order = fortran/column major.
    ///This operation is done inplace and does not create a new structure
    fn to_bunge_inplace(&self, bunge: &mut Bunge){
        let mut ori = bunge.ori_view_mut();

        let new_nelem = ori.len_of(Axis(1));
        let nelem = self.ori.len_of(Axis(1));

        assert!(new_nelem == nelem, 
        "The number of elements in the original ori field do no match up with the new field.
        The old field had {} elements, and the new field has {} elements",
        nelem, new_nelem);

        let tol = f64::sqrt(std::f64::EPSILON);

        azip!(mut bunge (ori.axis_iter_mut(Axis(1))), ref quat (self.ori.axis_iter(Axis(1))) in {
            let q03 = quat[0] * quat[0] + quat[3] * quat[3];
            let q12 = quat[1] * quat[1] + quat[2] * quat[2];
            let xi = f64::sqrt(q03 * q12);
            //We get to now go through all of the different cases that this might break down into
            if xi.abs() < tol && q12.abs() < tol {
                bunge[0] = f64::atan2(-2.0_f64 * quat[0] * quat[3], quat[0] * quat[0] - quat[3] * quat[3]);
                //All of the other values are zero
            }else if xi.abs() < tol && q03.abs() < tol{
                bunge[0] = f64::atan2(2.0_f64 * quat[1] * quat[2], quat[1] * quat[1] - quat[2] * quat[2]);
                bunge[1] = std::f64::consts::PI;
                //The other value is zero
            }else{
                let inv_xi = 1.0_f64 / xi;
                //The atan2 terms are pretty long so we're breaking it down into a couple of temp terms
                let t1 = inv_xi * (quat[1] * quat[3] - quat[0] * quat[2]);
                let t2 = inv_xi * (-quat[0] * quat[1] - quat[2] * quat[3]);
                //We can now assign the first two bunge angles
                bunge[0] = t1.atan2(t2);
                bunge[1] = f64::atan2(2.0_f64 * xi, q03 - q12);
                //Once again these terms going into the atan2 term are pretty long
                let t1 = inv_xi * (quat[0] * quat[2] + quat[1] * quat[3]);
                let t2 = inv_xi * (quat[2] * quat[3] - quat[0] * quat[1]);
                //We can finally find the final bunge angle
                bunge[2] = t1.atan2(t2);
            }
        });

    }

    ///Converts the unit quaternion representation over to rotation matrix which has the following properties
    ///shape (3, 3, nelems), memory order = fortran/column major.
    ///This operation is done inplace and does not create a new structure
    fn to_rmat_inplace(&self, rmat: &mut RMat){
        let mut ori = rmat.ori_view_mut();

        let new_nelem = ori.len_of(Axis(2));
        let nelem = self.ori.len_of(Axis(1));

        assert!(new_nelem == nelem, 
        "The number of elements in the original ori field do no match up with the new field.
        The old field had {} elements, and the new field has {} elements",
        nelem, new_nelem);

        azip!(mut rmat (ori.axis_iter_mut(Axis(2))), ref quat (self.ori.axis_iter(Axis(1))) in {
            let qbar =  quat[0] * quat[0] - (quat[1] * quat[1] + quat[2] * quat[2] + quat[3] * quat[3]);

            rmat[[0, 0]] = qbar + 2.0_f64 * quat[1] * quat[1];
            rmat[[1, 0]] = 2.0_f64 * (quat[1] * quat[2] + quat[0] * quat[3]);
            rmat[[2, 0]] = 2.0_f64 * (quat[1] * quat[3] - quat[0] * quat[2]);

            rmat[[0, 1]] = 2.0_f64 * (quat[1] * quat[2] - quat[0] * quat[3]);
            rmat[[1, 1]] = qbar + 2.0_f64 * quat[2] * quat[2];
            rmat[[2, 1]] = 2.0_f64 * (quat[2] * quat[3] + quat[0] * quat[1]);

            rmat[[0, 2]] = 2.0_f64 * (quat[1] * quat[3] + quat[0] * quat[2]);
            rmat[[1, 2]] = 2.0_f64 * (quat[2] * quat[3] - quat[0] * quat[1]);
            rmat[[2, 2]] = qbar + 2.0_f64 * quat[3] * quat[3];
        });
    }

    ///Converts the unit quaternion over to a angle-axis representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    ///This operation is done inplace and does not create a new structure
    fn to_ang_axis_inplace(&self, ang_axis: &mut AngAxis){
        let mut ori = ang_axis.ori_view_mut();

        let new_nelem = ori.len_of(Axis(1));
        let nelem = self.ori.len_of(Axis(1));

        assert!(new_nelem == nelem, 
        "The number of elements in the original ori field do no match up with the new field.
        The old field had {} elements, and the new field has {} elements",
        nelem, new_nelem);

        let tol = std::f64::EPSILON;

        azip!(mut angaxis (ori.axis_iter_mut(Axis(1))), ref quat (self.ori.axis_iter(Axis(1))) in {
            let phi = 2.0_f64 * quat[0].acos();
            if quat[0].abs() < tol{
                angaxis[0] = quat[1];
                angaxis[1] = quat[2];
                angaxis[2] = quat[3];
                angaxis[3] = std::f64::consts::PI;
            }else if phi.abs() < tol{
                angaxis[2] = 1.0_f64;
            }else{
                let s   = quat[0].signum() / f64::sqrt(quat[1] * quat[1] + quat[2] * quat[2] + quat[3] * quat[3]);

                angaxis[0] = s * quat[1];
                angaxis[1] = s * quat[2];
                angaxis[2] = s * quat[3];
                angaxis[3] = phi;
            }
        });
    }

    ///Converts the unit quaternion over to a compact angle-axis representation which has the following properties
    ///shape (3, nelems), memory order = fortran/column major.
    ///This operation is done inplace and does not create a new structure
    fn to_ang_axis_comp_inplace(&self, ang_axis_comp: &mut AngAxisComp){
        let mut ori = ang_axis_comp.ori_view_mut();

        let new_nelem = ori.len_of(Axis(1));
        let nelem = self.ori.len_of(Axis(1));

        assert!(new_nelem == nelem, 
        "The number of elements in the original ori field do no match up with the new field.
        The old field had {} elements, and the new field has {} elements",
        nelem, new_nelem);

        let tol = std::f64::EPSILON;

        azip!(mut angaxis (ori.axis_iter_mut(Axis(1))), ref quat (self.ori.axis_iter(Axis(1))) in {
            let phi = 2.0_f64 * quat[0].acos();
            if quat[0].abs() < tol{
                angaxis[0] = quat[1] * std::f64::consts::PI;
                angaxis[1] = quat[2] * std::f64::consts::PI;
                angaxis[2] = quat[3] * std::f64::consts::PI;
            }else{
                let s   = quat[0].signum() / f64::sqrt(quat[1] * quat[1] + quat[2] * quat[2] + quat[3] * quat[3]);

                angaxis[0] = s * quat[1] * phi;
                angaxis[1] = s * quat[2] * phi;
                angaxis[2] = s * quat[3] * phi; 
            }
        });
    }

    ///Converts the unit quaternion over to a Rodrigues vector representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    ///This operation is done inplace and does not create a new structure
    fn to_rod_vec_inplace(&self, rod_vec: &mut RodVec){
        let mut ori = rod_vec.ori_view_mut();

        let new_nelem = ori.len_of(Axis(1));
        let nelem = self.ori.len_of(Axis(1));

        assert!(new_nelem == nelem, 
        "The number of elements in the original ori field do no match up with the new field.
        The old field had {} elements, and the new field has {} elements",
        nelem, new_nelem);

        let tol = std::f64::EPSILON;

        azip!(mut rod_vec (ori.axis_iter_mut(Axis(1))), ref quat (self.ori.axis_iter(Axis(1))) in {
            let phi = quat[0].acos();
            if quat[0].abs() < tol{
                rod_vec[0] = quat[1];
                rod_vec[1] = quat[2];
                rod_vec[2] = quat[3];
                rod_vec[3] = std::f64::INFINITY;
            }else if phi.abs() < tol{
                rod_vec[2] = 1.0_f64;
            }else{
                let s   = quat[0].signum() / f64::sqrt(quat[1] * quat[1] + quat[2] * quat[2] + quat[3] * quat[3]);

                rod_vec[0] = s * quat[1];
                rod_vec[1] = s * quat[2];
                rod_vec[2] = s * quat[3];
                rod_vec[3] = phi.tan();
            }
        });
    }

    ///Converts the unit quaternion over to a compact Rodrigues vector representation which has the following properties
    ///shape (3, nelems), memory order = fortran/column major.
    ///This operation is done inplace and does not create a new structure
    fn to_rod_vec_comp_inplace(&self, rod_vec_comp: &mut RodVecComp){
        let mut ori = rod_vec_comp.ori_view_mut();

        let new_nelem = ori.len_of(Axis(1));
        let nelem = self.ori.len_of(Axis(1));

        assert!(new_nelem == nelem, 
        "The number of elements in the original ori field do no match up with the new field.
        The old field had {} elements, and the new field has {} elements",
        nelem, new_nelem);

        let tol = std::f64::EPSILON;

        azip!(mut rod_vec_comp (ori.axis_iter_mut(Axis(1))), ref quat (self.ori.axis_iter(Axis(1))) in {
            let tan_phi = f64::tan(quat[0].acos());
            //This case will not allow for anything to be retrievable later on...
            if quat[0].abs() < tol{
                rod_vec_comp[0] = std::f64::INFINITY;
                rod_vec_comp[1] = std::f64::INFINITY;
                rod_vec_comp[2] = std::f64::INFINITY;
            }else{
                let s   = quat[0].signum() / f64::sqrt(quat[1] * quat[1] + quat[2] * quat[2] + quat[3] * quat[3]);

                rod_vec_comp[0] = s * quat[1] * tan_phi;
                rod_vec_comp[1] = s * quat[2] * tan_phi;
                rod_vec_comp[2] = s * quat[3] * tan_phi; 
            }
        });
    }

    ///This returns a clone of the original unit quaternion structure
    ///This operation is done inplace and does not create a new structure
    fn to_quat_inplace(&self, quat: &mut Quat){
        let mut ori = quat.ori_view_mut();

        let new_nelem = ori.len_of(Axis(1));
        let nelem = self.ori.len_of(Axis(1));

        assert!(new_nelem == nelem, 
        "The number of elements in the original ori field do no match up with the new field.
        The old field had {} elements, and the new field has {} elements",
        nelem, new_nelem);

        ori.assign(&self.ori);
    }

    ///Converts the quaternion representation over to a homochoric representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    ///This operation is done inplace and does not create a new structure
    fn to_homochoric_inplace(&self, homochoric: &mut Homochoric){
        let ang_axis = self.to_ang_axis();
        ang_axis.to_homochoric_inplace(homochoric);    
    }


}//End of impl of unit Quaternion

///A series of commonly used operations to rotate vector data by a given rotation
impl RotVector for Quat{

    ///rot_vector takes in a 2D array view of a series of vectors. It then rotates these vectors using the
    ///given Quaternion. The newly rotated vectors are then returned. This function requires the
    ///number of elements in the Quaternion to be either 1.
    ///The unrotated vector might also contain either 1 or nelems number of elements.
    ///If this condition is not met the function will error out.
    ///vec - the vector to be rotated must have dimensions 3xnelems or 3x1.
    ///Output - the rotated vector and has dimensions 3xnelems.
    fn rot_vector(&self, vec: ArrayView2<f64>) -> Array2<f64>{

        let nelems = vec.len_of(Axis(1));
        let rnelems = self.ori.len_of(Axis(1));

        let rows  = vec.len_of(Axis(0));
        assert!((rows == 3), "The number of rows must be 3. The number of rows provided is {}", rows); 

        assert!( (nelems == rnelems) | (rnelems == 1) | (nelems == 1), 
        "The number of elements in the vector field must be equal to the number of elements in the
        Quaternion structure, or their must only be one element in Quaternion. The final case is
        that there must only be one element in the vector field. There are
        currently {} elements in vector and {} elements in Quaternion",
        nelems, rnelems);

        let mnelems = cmp::max(rnelems, nelems);
        let mut rvec = Array2::<f64>::zeros((3, mnelems).f());

        //We need to see if we have more than one Quaternion that we're multiplying by
        if rnelems == nelems {
            //The rotations here can be given by reference 1  equation 24 in the README.
            azip!(mut rvec (rvec.axis_iter_mut(Axis(1))), ref vec (vec.axis_iter(Axis(1))), 
            ref quat (self.ori.axis_iter(Axis(1))) in {
                quat_rot_vec(&quat, &vec, rvec);     
            });
        } else if rnelems == 1{
            //We just have one Quaternion so perform pretty much the above to get all of our values
            let quat = self.ori.subview(Axis(1), 0);

            azip!(mut rvec (rvec.axis_iter_mut(Axis(1))), ref vec (vec.axis_iter(Axis(1))) in {  
                quat_rot_vec(&quat, &vec, rvec);      
            });
        }else{
            //We just have one vector so perform pretty much the above to get all of our values
            let vec = vec.subview(Axis(1), 0);

            azip!(mut rvec (rvec.axis_iter_mut(Axis(1))), ref quat (self.ori.axis_iter(Axis(1))) in {  
                quat_rot_vec(&quat, &vec, rvec);  
            });
        }//End of if-else
        //Now we just need to return the rvec value
        rvec
    }//End of rot_vector

    ///rot_vector_mut takes in a 2D array view of a series of vectors and a mutable 2D ArrayView of the 
    ///rotated vector. It then rotates these vectors using the given Quaternion. The newly rotated
    /// vectors are assigned to the supplied rotated vector, rvec. This function requires the
    ///number of elements in the Quaternion to be either 1 or nelems.
    ///The unrotated vector might also contain either 1 or nelems number of elements.
    ///It also requires the number of elements in rvec and vec to be equal.
    ///If these conditions are not met the function will error out.
    ///vec - the vector to be rotated must have dimensions 3xnelems or 3x1.
    ///rvec - the rotated vector and has dimensions 3xnelems.
    fn rot_vector_mut(&self, vec: ArrayView2<f64>, mut rvec: ArrayViewMut2<f64>) {

        let nelems = vec.len_of(Axis(1));
        let rvnelems = rvec.len_of(Axis(1));
        let rnelems = self.ori.len_of(Axis(1));
        let mnelems = cmp::max(rnelems, nelems);

        let rows  = vec.len_of(Axis(0));
        assert!((rows == 3), "The number of rows must be 3. The number of rows provided is {}", rows); 

        assert!((mnelems == rvnelems),
        "The number of elements in the unrotated vector or quaternion field must be equal to the number of elements
        in the supplied rotated vector field. There are currently {} elements in the unrotated vector or quaternion
        field and {} elements in the rotated vector field", 
        mnelems, rvnelems);

        assert!( (nelems == rnelems) | (rnelems == 1) | (nelems == 1), 
        "The number of elements in the vector field must be equal to the number of elements in the
        Quaternion structure, or their must only be one element in Quaternion. The final case is
        that there must only be one element in the vector field. There are
        currently {} elements in vector and {} elements in Quaternion",
        nelems, rnelems);

        //We need to see if we have more than one Quaternion that we're multiplying by
        if rnelems == nelems {
            //The rotations here can be given by reference 1  equation 24 in the README.
            azip!(mut rvec (rvec.axis_iter_mut(Axis(1))), ref vec (vec.axis_iter(Axis(1))), 
            ref quat (self.ori.axis_iter(Axis(1))) in {
                quat_rot_vec(&quat, &vec, rvec);         
            });
        } else if rnelems == 1{
            //We just have one Quaternion so perform pretty much the above to get all of our values
            let quat = self.ori.subview(Axis(1), 0);

            azip!(mut rvec (rvec.axis_iter_mut(Axis(1))), ref vec (vec.axis_iter(Axis(1))) in {  
                quat_rot_vec(&quat, &vec, rvec);  
            });
        } else{
            //We just have one vector so perform pretty much the above to get all of our values
            let vec = vec.subview(Axis(1), 0);

            azip!(mut rvec (rvec.axis_iter_mut(Axis(1))), ref quat (self.ori.axis_iter(Axis(1))) in {  
                quat_rot_vec(&quat, &vec, rvec);  
            });
        }//End of if-else
    }//End of rot_vector_mut

    ///rot_vector_inplace takes in a mutable 2D array view of a series of vectors. It then rotates these vectors using the
    ///given Quaternion. The newly rotated vectors are assigned to original vector. This function requires the
    ///number of elements in the Quaternion to be either 1 or nelems where vec has nelems in it.
    ///If this condition is not met the function will error out.
    ///vec - the vector to be rotated must have dimensions 3xnelems.
    fn rot_vector_inplace(&self, mut vec: ArrayViewMut2<f64>){

        let nelems = vec.len_of(Axis(1));
        let rnelems = self.ori.len_of(Axis(1));

        let rows  = vec.len_of(Axis(0));
        assert!((rows == 3), "The number of rows must be 3. The number of rows provided is {}", rows); 

        assert!( (nelems == rnelems) | (rnelems == 1), 
        "The number of elements in the vector field must be equal to the number of elements in the
        Quaternion structure, or their must only be one element in Quaternion. There are
        currently {} elements in vector and {} elements in Quaternion",
        nelems, rnelems);

        //We need to see if we have more than one Quaternion that we're multiplying by
        if rnelems == nelems {
            //The rotations here can be given by reference 1  equation 24 in the README.
            azip!(mut vec (vec.axis_iter_mut(Axis(1))), ref quat (self.ori.axis_iter(Axis(1))) in {
                let mut rvec = Array1::<f64>::zeros((3).f());
                quat_rot_vec(&quat, &vec.view(), rvec.view_mut());
                vec.assign({&rvec});    
            });
        } else{
            //We just have one Quaternion so perform pretty much the above to get all of our values
            let quat = self.ori.subview(Axis(1), 0);

            azip!(mut vec (vec.axis_iter_mut(Axis(1))) in {
                let mut rvec = Array1::<f64>::zeros((3).f()); 
                quat_rot_vec(&quat, &vec.view(), rvec.view_mut());
                vec.assign({&rvec});  
            });
        }//End of if-else
    }//End of rot_vector_inplace
}//Endo of Impl RotVector

///All of the quaternion vector rotation operations can be described by using the below series of functions.
///This also reduces the amount of repetive code that existed earlier within rot_vector. 
fn quat_rot_vec(quat: &ArrayView1<f64>, vec: &ArrayView1<f64>, mut rvec: ArrayViewMut1<f64>){
    let q02 = 2.0_f64 * quat[0];
    //(q_0^2 - ||q||^2)
    let q02_m_nq = quat[0] * quat[0] - (quat[1] * quat[1] + quat[2] * quat[2] + quat[3] * quat[3]);
    let dot_prod2 = 2.0_f64 * (quat[1] * vec[0] + quat[2] * vec[1] + quat[3] * vec[2]);
    let mut cross_prod = Array1::<f64>::zeros((3).f());

    cross_prod[0] = -quat[3] * vec[1] + quat[2] * vec[2];
    cross_prod[1] = quat[3] * vec[0] - quat[1] * vec[2];
    cross_prod[2] = -quat[2] * vec[0] + quat[1] * vec[1];

    rvec[0] = vec[0] * q02_m_nq + cross_prod[0] * q02 + quat[1] * dot_prod2;
    rvec[1] = vec[1] * q02_m_nq + cross_prod[1] * q02 + quat[2] * dot_prod2;
    rvec[2] = vec[2] * q02_m_nq + cross_prod[2] * q02 + quat[3] * dot_prod2;  
}