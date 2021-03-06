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

        #[cfg(feature = "parallel")]
        par_azip!((mut quat in ori.axis_iter_mut(Axis(1))) {quat[0] = 1.0_f64});

        #[cfg(not(feature = "parallel"))]
        azip!((mut quat in ori.axis_iter_mut(Axis(1))) {quat[0] = 1.0_f64});

        Quat{
            ori,
        }
    }//End of new

    ///Creates a unit quaternion type with the supplied data as long as the supplied data is in the following format
    ///shape (4, nelems), memory order = fortran/column major.
    ///If it doesn't fit those standards it will fail.
    pub fn new_init(ori: Array2<f64>) -> Quat{

        let nrow = ori.nrows();

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

        #[cfg(feature = "parallel")]
        par_azip!((mut quat_c in ori.axis_iter_mut(Axis(1)), ref quat in self.ori.axis_iter(Axis(1))) {
            quat_c[0] = quat[0];
            quat_c[1] = -1.0_f64 * quat[1];
            quat_c[2] = -1.0_f64 * quat[2];
            quat_c[3] = -1.0_f64 * quat[3];
        });

        #[cfg(not(feature = "parallel"))]
        azip!((mut quat_c in ori.axis_iter_mut(Axis(1)), ref quat in self.ori.axis_iter(Axis(1))) {
            quat_c[0] = quat[0];
            quat_c[1] = -1.0_f64 * quat[1];
            quat_c[2] = -1.0_f64 * quat[2];
            quat_c[3] = -1.0_f64 * quat[3];
        });

        Quat::new_init(ori)
    }//End of conjugate

    ///Performs in place the conjugate/inverse of the unit quaternion which is simply the negative
    ///of the vector portions of the unit quaternion. The inverse is said to be the same here because
    ///for unit quaternions that is the case. If we didn't have unit quaternions that would not be the case.
    pub fn conjugate_inplace(&mut self){

        #[cfg(feature = "parallel")]
        par_azip!((mut quat_c in self.ori.axis_iter_mut(Axis(1))) {
            quat_c[1] *= -1.0_f64;
            quat_c[2] *= -1.0_f64;
            quat_c[3] *= -1.0_f64;
        });

        #[cfg(not(feature = "parallel"))]
        azip!((mut quat_c in self.ori.axis_iter_mut(Axis(1))) {
            quat_c[1] *= -1.0_f64;
            quat_c[2] *= -1.0_f64;
            quat_c[3] *= -1.0_f64;
        });
    }//End of conjugate_inplace

    ///Performs a quaternion product operation between two unit quaternions ->
    ///q_new = (q_01*q_02 - q_1.q_2, q_01*q_2 +q_02*q_1 + q_1 x q_2) where q_1 and q_2 are the vector components of the
    ///quaternion. The result returned is a new quaternion. Also, if the conjugate/inverse is used for the first quaternion
    ///one can obtain the relative orientation/rotation between quaternion 1 and quaternion 2.
    ///This function requires the number of elements in self to be either 1.
    ///The quat2 field might also contain either 1 or nelems number of elements.
    ///If this condition is not met the function will error out.
    ///quat2 - the quaternion to be rotated must have dimensions 4xnelems or 4x1.
    ///Output - the quaternion product and has dimensions 4xnelems.
    pub fn product(&self, quat2: &Quat) -> Quat{
        
        let ori_quat2 = quat2.ori_view();
        let nelems = ori_quat2.len_of(Axis(1));
        let rnelems = self.ori.len_of(Axis(1));

        assert!( (nelems == rnelems) | (rnelems == 1) | (nelems == 1), 
        "The number of elements in quat2 field must be equal to the number of elements in the
        Quaternion structure, or their must only be one element in Quaternion. The final case is
        that there must only be one element in the quat2 field. There are
        currently {} elements in quat2 and {} elements in Quaternion",
        nelems, rnelems);

        let mnelems = cmp::max(rnelems, nelems);
        let mut quat_prod = Array2::<f64>::zeros((4, mnelems).f());

                //We need to see if we have more than one Quaternion that we're multiplying by
        if rnelems == nelems {
            //The rotations here can be given by reference 1  equation 24 in the README.

            #[cfg(feature = "parallel")]
            par_azip!((quat_prod in quat_prod.axis_iter_mut(Axis(1)), ref quat2 in ori_quat2.axis_iter(Axis(1)), 
            ref quat1 in self.ori.axis_iter(Axis(1))) {
                quat_product(&quat1, &quat2, quat_prod);     
            });

            #[cfg(not(feature = "parallel"))]
            azip!((quat_prod in quat_prod.axis_iter_mut(Axis(1)), ref quat2 in ori_quat2.axis_iter(Axis(1)), 
            ref quat1 in self.ori.axis_iter(Axis(1))) {
                quat_product(&quat1, &quat2, quat_prod);     
            });
        } else if rnelems == 1{
            //We just have one Quaternion so perform pretty much the above to get all of our values
            let quat1 = self.ori.index_axis(Axis(1), 0);

            #[cfg(feature = "parallel")]
            par_azip!((quat_prod in quat_prod.axis_iter_mut(Axis(1)), ref quat2 in ori_quat2.axis_iter(Axis(1))) {  
                quat_product(&quat1, &quat2, quat_prod);      
            });

            #[cfg(not(feature = "parallel"))]
            azip!((quat_prod in quat_prod.axis_iter_mut(Axis(1)), ref quat2 in ori_quat2.axis_iter(Axis(1))) {  
                quat_product(&quat1, &quat2, quat_prod);      
            });
        }else{
            //We just have one vector so perform pretty much the above to get all of our values
            let quat2 = ori_quat2.index_axis(Axis(1), 0);

            #[cfg(feature = "parallel")]
            par_azip!((quat_prod in quat_prod.axis_iter_mut(Axis(1)), ref quat1 in self.ori.axis_iter(Axis(1))) {  
                quat_product(&quat1, &quat2, quat_prod);  
            });

            #[cfg(not(feature = "parallel"))]
            azip!((quat_prod in quat_prod.axis_iter_mut(Axis(1)), ref quat1 in self.ori.axis_iter(Axis(1))) {  
                quat_product(&quat1, &quat2, quat_prod);  
            });
        }//End of if-else

        Quat::new_init(quat_prod)
    }//End of product

    ///Performs a quaternion product operation between two unit quaternions ->
    ///q_new = (q_01*q_02 - q_1.q_2, q_01*q_2 +q_02*q_1 + q_1 x q_2) where q_1 and q_2 are the vector components of the
    ///quaternion. The result is stored in a supplied quaternion field. 
    ///Also, if the conjugate/inverse is used for the first quaternion
    ///one can obtain the relative orientation/rotation between quaternion 1 and quaternion 2.
    ///This function requires the number of elements in self to be either 1.
    ///The quat2 field might also contain either 1 or nelems number of elements.
    ///The quat_prod field must contain nelems number of elements.
    ///If this condition is not met the function will error out.
    ///quat2 - the quaternion to be rotated must have dimensions 4xnelems or 4x1.
    ///quat_prod - the quaternion product that was supplied that we are going to store data in must have dims 4xnelems
    pub fn product_mut(&self, quat2: &Quat, quat_prod: &mut Quat){
        
        let ori_quat2 = quat2.ori_view();
        let mut ori_quat_prod = quat_prod.ori_view_mut();

        let nelems = ori_quat2.len_of(Axis(1));
        let rvnelems = ori_quat_prod.len_of(Axis(1));
        let rnelems = self.ori.len_of(Axis(1));

        let mnelems = cmp::max(rnelems, nelems);

        assert!((mnelems == rvnelems),
        "The number of elements in the quat2 or Quaternion field must be equal to the number of elements
        in the supplied quat_prod field. There are currently {} elements in the quat2 or Quaternion
        field and {} elements in the quat_prod field", 
        mnelems, rvnelems);

        assert!( (nelems == rnelems) | (rnelems == 1) | (nelems == 1), 
        "The number of elements in quat2 field must be equal to the number of elements in the
        Quaternion structure, or their must only be one element in Quaternion. The final case is
        that there must only be one element in the quat2 field. There are
        currently {} elements in quat2 and {} elements in Quaternion",
        nelems, rnelems);

        //We need to see if we have more than one Quaternion that we're multiplying by
        if rnelems == nelems {
            //The rotations here can be given by reference 1  equation 23 in the README.
            #[cfg(feature = "parallel")]
            par_azip!((quat_prod in ori_quat_prod.axis_iter_mut(Axis(1)), ref quat2 in ori_quat2.axis_iter(Axis(1)), 
            ref quat1 in self.ori.axis_iter(Axis(1))) {
                quat_product(&quat1, &quat2, quat_prod);     
            });

            #[cfg(not(feature = "parallel"))]
            azip!((quat_prod in ori_quat_prod.axis_iter_mut(Axis(1)), ref quat2 in ori_quat2.axis_iter(Axis(1)), 
            ref quat1 in self.ori.axis_iter(Axis(1))) {
                quat_product(&quat1, &quat2, quat_prod);     
            });
        } else if rnelems == 1{
            //We just have one Quaternion so perform pretty much the above to get all of our values
            let quat1 = self.ori.index_axis(Axis(1), 0);

            #[cfg(feature = "parallel")]
            par_azip!((quat_prod in ori_quat_prod.axis_iter_mut(Axis(1)), ref quat2 in ori_quat2.axis_iter(Axis(1))) {  
                quat_product(&quat1, &quat2, quat_prod);      
            });

            #[cfg(not(feature = "parallel"))]
            azip!((quat_prod in ori_quat_prod.axis_iter_mut(Axis(1)), ref quat2 in ori_quat2.axis_iter(Axis(1))) {  
                quat_product(&quat1, &quat2, quat_prod);      
            });
        }else{
            //We just have one vector so perform pretty much the above to get all of our values
            let quat2 = ori_quat2.index_axis(Axis(1), 0);

            #[cfg(feature = "parallel")]
            par_azip!((quat_prod in ori_quat_prod.axis_iter_mut(Axis(1)), ref quat1 in self.ori.axis_iter(Axis(1))) {  
                quat_product(&quat1, &quat2, quat_prod);  
            });

            #[cfg(not(feature = "parallel"))]
            azip!((quat_prod in ori_quat_prod.axis_iter_mut(Axis(1)), ref quat1 in self.ori.axis_iter(Axis(1))) {  
                quat_product(&quat1, &quat2, quat_prod);  
            });
        }//End of if-else
    }//End of product_mut

}//End of Impl of Quat

//A helper function for Impl of Quat

///All of the quaternion product operations can be described by using the below series of functions.
///q_new = (q_01*q_02 - q_1.q_2, q_01*q_2 +q_02*q_1 + q_1 x q_2) where q_1 and q_2 are the vector components of the
///quaternion. 
fn quat_product(quat1: &ArrayView1<f64>, quat2: &ArrayView1<f64>, mut quat_prod: ArrayViewMut1<f64>){
    let q01q02 = quat1[0] * quat2[0];
    //(q_0^2 - ||q||^2)
    let q01q02_qd = q01q02 - (quat1[1] * quat2[1] + quat1[2] * quat2[2] + quat1[3] * quat2[3]);
    let mut cross_prod = Array1::<f64>::zeros((3).f());

    cross_prod[0] = -quat1[3] * quat2[2] + quat1[2] * quat2[3];
    cross_prod[1] = quat1[3] * quat2[1] - quat1[1] * quat2[3];
    cross_prod[2] = -quat1[2] * quat2[1] + quat1[1] * quat2[2];

    quat_prod[0] = q01q02_qd;
    quat_prod[1] = quat1[0] * quat2[1] + quat2[0] * quat1[1] + cross_prod[0];
    quat_prod[2] = quat1[0] * quat2[2] + quat2[0] * quat1[2] + cross_prod[0];
    quat_prod[3] = quat1[0] * quat2[3] + quat2[0] * quat1[3] + cross_prod[0];
}//End of quat_product

///The orientation conversions of a series of unit quaternions to a number of varying different orientation
///representations commonly used in material orientation processing. 
impl OriConv for Quat{
    ///Converts the unit quaternion representation over to Bunge angles which has the following properties
    ///shape (3, nelems), memory order = fortran/column major.
    fn to_bunge(&self) -> Bunge{

        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((3, nelems).f());

        let tol = f64::sqrt(std::f64::EPSILON);

        let f = |mut bunge: ArrayViewMut1::<f64>, ref quat: ArrayView1::<f64>| {
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
        };

        #[cfg(feature = "parallel")]
        par_azip!((bunge in ori.axis_iter_mut(Axis(1)), quat in self.ori.axis_iter(Axis(1))) {
            f(bunge, quat);
        });

        #[cfg(not(feature = "parallel"))]
        azip!((bunge in ori.axis_iter_mut(Axis(1)), quat in self.ori.axis_iter(Axis(1))) {
            f(bunge, quat);
        });

        Bunge::new_init(ori)
    }//End of to_bunge

    ///Converts the unit quaternion representation over to rotation matrix which has the following properties
    ///shape (3, 3, nelems), memory order = fortran/column major.
    fn to_rmat(&self) -> RMat{

        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array3::<f64>::zeros((3, 3, nelems).f());

        let f = |mut rmat: ArrayViewMut2::<f64>, ref quat: ArrayView1::<f64>| {
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
        };

        #[cfg(feature = "parallel")]
        par_azip!((rmat in ori.axis_iter_mut(Axis(2)), quat in self.ori.axis_iter(Axis(1))) {
            f(rmat, quat);
        });

        #[cfg(not(feature = "parallel"))]
        azip!((rmat in ori.axis_iter_mut(Axis(2)), quat in self.ori.axis_iter(Axis(1))) {
            f(rmat, quat);
        });

        RMat::new_init(ori)
    }//End of to_rmat

    ///Converts the unit quaternion representation over to angle-axis representation which has the following properties
    ///shape (3, nelems), memory order = fortran/column major.
    fn to_ang_axis(&self) -> AngAxis{

        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((4, nelems).f());

        let tol = std::f64::EPSILON;

        let f = |mut ang_axis: ArrayViewMut1::<f64>, ref quat: ArrayView1::<f64>| {
            let phi = 2.0_f64 * quat[0].acos();
            if quat[0].abs() < tol{
                ang_axis[0] = quat[1];
                ang_axis[1] = quat[2];
                ang_axis[2] = quat[3];
                ang_axis[3] = std::f64::consts::PI;
            }else if phi.abs() < tol{
                ang_axis[2] = 1.0_f64;
            }else{
                let s   = quat[0].signum() / f64::sqrt(quat[1] * quat[1] + quat[2] * quat[2] + quat[3] * quat[3]);

                ang_axis[0] = s * quat[1];
                ang_axis[1] = s * quat[2];
                ang_axis[2] = s * quat[3];
                ang_axis[3] = phi;
            }
        };

        #[cfg(feature = "parallel")]
        par_azip!((ang_axis in ori.axis_iter_mut(Axis(1)), quat in self.ori.axis_iter(Axis(1))) {
            f(ang_axis, quat);
        });

        #[cfg(not(feature = "parallel"))]
        azip!((ang_axis in ori.axis_iter_mut(Axis(1)), quat in self.ori.axis_iter(Axis(1))) {
            f(ang_axis, quat);
        });

        AngAxis::new_init(ori)
    }//End of to_ang_axis

    ///Converts the unit quaternion over to a compact angle-axis representation which has the following properties
    ///shape (3, nelems), memory order = fortran/column major.
    fn to_ang_axis_comp(&self) -> AngAxisComp{

        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((3, nelems).f());

        let tol = std::f64::EPSILON;

        let f = |mut ang_axis: ArrayViewMut1::<f64>, ref quat: ArrayView1::<f64>| {
            let phi = 2.0_f64 * quat[0].acos();
            if quat[0].abs() < tol{
                ang_axis[0] = quat[1] * std::f64::consts::PI;
                ang_axis[1] = quat[2] * std::f64::consts::PI;
                ang_axis[2] = quat[3] * std::f64::consts::PI;
            }else{
                let s   = quat[0].signum() / f64::sqrt(quat[1] * quat[1] + quat[2] * quat[2] + quat[3] * quat[3]);

                ang_axis[0] = s * quat[1] * phi;
                ang_axis[1] = s * quat[2] * phi;
                ang_axis[2] = s * quat[3] * phi; 
            }
        };

        #[cfg(feature = "parallel")]
        par_azip!((ang_axis in ori.axis_iter_mut(Axis(1)), quat in self.ori.axis_iter(Axis(1))) {
            f(ang_axis, quat);
        });

        #[cfg(not(feature = "parallel"))]
        azip!((ang_axis in ori.axis_iter_mut(Axis(1)), quat in self.ori.axis_iter(Axis(1))) {
            f(ang_axis, quat);
        });

        AngAxisComp::new_init(ori)
    }//End of to_ang_axis_comp

    ///Converts the unit quaternion over to a Rodrigues vector representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    fn to_rod_vec(&self) -> RodVec{

        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((4, nelems).f());

        let tol = std::f64::EPSILON;

        let f = |mut rod_vec: ArrayViewMut1::<f64>, ref quat: ArrayView1::<f64>| {
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
        };

        #[cfg(feature = "parallel")]
        par_azip!((rod_vec in ori.axis_iter_mut(Axis(1)), quat in self.ori.axis_iter(Axis(1))) {
            f(rod_vec, quat);
        });

        #[cfg(not(feature = "parallel"))]
        azip!((rod_vec in ori.axis_iter_mut(Axis(1)), quat in self.ori.axis_iter(Axis(1))) {
            f(rod_vec, quat);
        });

        RodVec::new_init(ori)
    }//End of to_rod_vec

    ///Converts the unit quaternion over to a compact Rodrigues vector representation which has the following properties
    ///shape (3, nelems), memory order = fortran/column major.
    fn to_rod_vec_comp(&self) -> RodVecComp{
        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((3, nelems).f());
        let tol = std::f64::EPSILON;

        let f = |mut rod_vec_comp: ArrayViewMut1::<f64>, ref quat: ArrayView1::<f64>| {
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
        };

        #[cfg(feature = "parallel")]
        par_azip!((rod_vec_comp in ori.axis_iter_mut(Axis(1)), quat in self.ori.axis_iter(Axis(1))) {
            f(rod_vec_comp, quat);
        });

        #[cfg(not(feature = "parallel"))]
        azip!((rod_vec_comp in ori.axis_iter_mut(Axis(1)), quat in self.ori.axis_iter(Axis(1))) {
            f(rod_vec_comp, quat);
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

        let f = |mut bunge: ArrayViewMut1::<f64>, ref quat: ArrayView1::<f64>| {
            let q03 = quat[0] * quat[0] + quat[3] * quat[3];
            let q12 = quat[1] * quat[1] + quat[2] * quat[2];
            let xi = f64::sqrt(q03 * q12);
            //We get to now go through all of the different cases that this might break down into
            if xi.abs() < tol && q12.abs() < tol {
                bunge[0] = f64::atan2(-2.0_f64 * quat[0] * quat[3], quat[0] * quat[0] - quat[3] * quat[3]);
                bunge[1] = 0.0_f64;
                bunge[2] = 0.0_f64;
            }else if xi.abs() < tol && q03.abs() < tol{
                bunge[0] = f64::atan2(2.0_f64 * quat[1] * quat[2], quat[1] * quat[1] - quat[2] * quat[2]);
                bunge[1] = std::f64::consts::PI;
                bunge[2] = 0.0_f64;
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
        };

        #[cfg(feature = "parallel")]
        par_azip!((bunge in ori.axis_iter_mut(Axis(1)), quat in self.ori.axis_iter(Axis(1))) {
            f(bunge, quat);
        });

        #[cfg(not(feature = "parallel"))]
        azip!((bunge in ori.axis_iter_mut(Axis(1)), quat in self.ori.axis_iter(Axis(1))) {
            f(bunge, quat);
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

        let f = |mut rmat: ArrayViewMut2::<f64>, ref quat: ArrayView1::<f64>| {
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
        };

        #[cfg(feature = "parallel")]
        par_azip!((rmat in ori.axis_iter_mut(Axis(2)), quat in self.ori.axis_iter(Axis(1))) {
            f(rmat, quat);
        });

        #[cfg(not(feature = "parallel"))]
        azip!((rmat in ori.axis_iter_mut(Axis(2)), quat in self.ori.axis_iter(Axis(1))) {
            f(rmat, quat);
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

        let f = |mut ang_axis: ArrayViewMut1::<f64>, ref quat: ArrayView1::<f64>| {
            let phi = 2.0_f64 * quat[0].acos();
            if quat[0].abs() < tol{
                ang_axis[0] = quat[1];
                ang_axis[1] = quat[2];
                ang_axis[2] = quat[3];
                ang_axis[3] = std::f64::consts::PI;
            }else if phi.abs() < tol{
                ang_axis[0] = 0.0_f64;
                ang_axis[1] = 0.0_f64;
                ang_axis[2] = 1.0_f64;
                ang_axis[3] = 0.0_f64;
            }else{
                let s   = quat[0].signum() / f64::sqrt(quat[1] * quat[1] + quat[2] * quat[2] + quat[3] * quat[3]);

                ang_axis[0] = s * quat[1];
                ang_axis[1] = s * quat[2];
                ang_axis[2] = s * quat[3];
                ang_axis[3] = phi;
            }
        };

        #[cfg(feature = "parallel")]
        par_azip!((ang_axis in ori.axis_iter_mut(Axis(1)), quat in self.ori.axis_iter(Axis(1))) {
            f(ang_axis, quat);
        });

        #[cfg(not(feature = "parallel"))]
        azip!((ang_axis in ori.axis_iter_mut(Axis(1)), quat in self.ori.axis_iter(Axis(1))) {
            f(ang_axis, quat);
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

        let f = |mut ang_axis: ArrayViewMut1::<f64>, ref quat: ArrayView1::<f64>| {
            let phi = 2.0_f64 * quat[0].acos();
            if quat[0].abs() < tol{
                ang_axis[0] = quat[1] * std::f64::consts::PI;
                ang_axis[1] = quat[2] * std::f64::consts::PI;
                ang_axis[2] = quat[3] * std::f64::consts::PI;
            }else{
                let s   = quat[0].signum() / f64::sqrt(quat[1] * quat[1] + quat[2] * quat[2] + quat[3] * quat[3]);

                ang_axis[0] = s * quat[1] * phi;
                ang_axis[1] = s * quat[2] * phi;
                ang_axis[2] = s * quat[3] * phi; 
            }
        };

        #[cfg(feature = "parallel")]
        par_azip!((ang_axis in ori.axis_iter_mut(Axis(1)), quat in self.ori.axis_iter(Axis(1))) {
            f(ang_axis, quat);
        });

        #[cfg(not(feature = "parallel"))]
        azip!((ang_axis in ori.axis_iter_mut(Axis(1)), quat in self.ori.axis_iter(Axis(1))) {
            f(ang_axis, quat);
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

        let f = |mut rod_vec: ArrayViewMut1::<f64>, ref quat: ArrayView1::<f64>| {
            let phi = quat[0].acos();
            if quat[0].abs() < tol{
                rod_vec[0] = quat[1];
                rod_vec[1] = quat[2];
                rod_vec[2] = quat[3];
                rod_vec[3] = std::f64::INFINITY;
            }else if phi.abs() < tol{
                rod_vec[0] = 0.0_f64;
                rod_vec[1] = 0.0_f64;
                rod_vec[2] = 1.0_f64;
                rod_vec[3] = 0.0_f64;
            }else{
                let s   = quat[0].signum() / f64::sqrt(quat[1] * quat[1] + quat[2] * quat[2] + quat[3] * quat[3]);

                rod_vec[0] = s * quat[1];
                rod_vec[1] = s * quat[2];
                rod_vec[2] = s * quat[3];
                rod_vec[3] = phi.tan();
            }
        };

        #[cfg(feature = "parallel")]
        par_azip!((rod_vec in ori.axis_iter_mut(Axis(1)), quat in self.ori.axis_iter(Axis(1))) {
            f(rod_vec, quat);
        });

        #[cfg(not(feature = "parallel"))]
        azip!((rod_vec in ori.axis_iter_mut(Axis(1)), quat in self.ori.axis_iter(Axis(1))) {
            f(rod_vec, quat);
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

        let f = |mut rod_vec_comp: ArrayViewMut1::<f64>, ref quat: ArrayView1::<f64>| {
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
        };

        #[cfg(feature = "parallel")]
        par_azip!((rod_vec_comp in ori.axis_iter_mut(Axis(1)), quat in self.ori.axis_iter(Axis(1))) {
            f(rod_vec_comp, quat);
        });

        #[cfg(not(feature = "parallel"))]
        azip!((rod_vec_comp in ori.axis_iter_mut(Axis(1)), quat in self.ori.axis_iter(Axis(1))) {
            f(rod_vec_comp, quat);
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

            #[cfg(feature = "parallel")]
            par_azip!((rvec in rvec.axis_iter_mut(Axis(1)), ref vec in vec.axis_iter(Axis(1)), 
            ref quat in self.ori.axis_iter(Axis(1))) {
                quat_rot_vec(&quat, &vec, rvec);     
            });

            #[cfg(not(feature = "parallel"))]
            azip!((rvec in rvec.axis_iter_mut(Axis(1)), ref vec in vec.axis_iter(Axis(1)), 
            ref quat in self.ori.axis_iter(Axis(1))) {
                quat_rot_vec(&quat, &vec, rvec);     
            });
        } else if rnelems == 1{
            //We just have one Quaternion so perform pretty much the above to get all of our values
            let quat = self.ori.index_axis(Axis(1), 0);

            #[cfg(feature = "parallel")]
            par_azip!((rvec in rvec.axis_iter_mut(Axis(1)), ref vec in vec.axis_iter(Axis(1))) {  
                quat_rot_vec(&quat, &vec, rvec);      
            });

            #[cfg(not(feature = "parallel"))]
            azip!((rvec in rvec.axis_iter_mut(Axis(1)), ref vec in vec.axis_iter(Axis(1))) {  
                quat_rot_vec(&quat, &vec, rvec);      
            });
        }else{
            //We just have one vector so perform pretty much the above to get all of our values
            let vec = vec.index_axis(Axis(1), 0);

            #[cfg(feature = "parallel")]
            par_azip!((rvec in rvec.axis_iter_mut(Axis(1)), ref quat in self.ori.axis_iter(Axis(1))) {  
                quat_rot_vec(&quat, &vec, rvec);  
            });

            #[cfg(not(feature = "parallel"))]
            azip!((rvec in rvec.axis_iter_mut(Axis(1)), ref quat in self.ori.axis_iter(Axis(1))) {  
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

            #[cfg(feature = "parallel")]
            par_azip!((rvec in rvec.axis_iter_mut(Axis(1)), ref vec in vec.axis_iter(Axis(1)), 
            ref quat in self.ori.axis_iter(Axis(1))) {
                quat_rot_vec(&quat, &vec, rvec);         
            });

            #[cfg(not(feature = "parallel"))]
            azip!((rvec in rvec.axis_iter_mut(Axis(1)), ref vec in vec.axis_iter(Axis(1)), 
            ref quat in self.ori.axis_iter(Axis(1))) {
                quat_rot_vec(&quat, &vec, rvec);         
            });
        } else if rnelems == 1{
            //We just have one Quaternion so perform pretty much the above to get all of our values
            let quat = self.ori.index_axis(Axis(1), 0);

            #[cfg(feature = "parallel")]
            par_azip!((rvec in rvec.axis_iter_mut(Axis(1)), ref vec in vec.axis_iter(Axis(1))) {  
                quat_rot_vec(&quat, &vec, rvec);  
            });

            #[cfg(not(feature = "parallel"))]
            azip!((rvec in rvec.axis_iter_mut(Axis(1)), ref vec in vec.axis_iter(Axis(1))) {  
                quat_rot_vec(&quat, &vec, rvec);  
            });
        } else{
            //We just have one vector so perform pretty much the above to get all of our values
            let vec = vec.index_axis(Axis(1), 0);

            #[cfg(feature = "parallel")]
            par_azip!((rvec in rvec.axis_iter_mut(Axis(1)), ref quat in self.ori.axis_iter(Axis(1))) {  
                quat_rot_vec(&quat, &vec, rvec);  
            });

            #[cfg(not(feature = "parallel"))]
            azip!((rvec in rvec.axis_iter_mut(Axis(1)), ref quat in self.ori.axis_iter(Axis(1))) {  
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

            #[cfg(feature = "parallel")]
            par_azip!((mut vec in vec.axis_iter_mut(Axis(1)), ref quat in self.ori.axis_iter(Axis(1))) {
                let mut rvec = Array1::<f64>::zeros((3).f());
                quat_rot_vec(&quat, &vec.view(), rvec.view_mut());
                vec.assign(&rvec);    
            });

            #[cfg(not(feature = "parallel"))]
            azip!((mut vec in vec.axis_iter_mut(Axis(1)), ref quat in self.ori.axis_iter(Axis(1))) {
                let mut rvec = Array1::<f64>::zeros((3).f());
                quat_rot_vec(&quat, &vec.view(), rvec.view_mut());
                vec.assign(&rvec);    
            });
        } else{
            //We just have one Quaternion so perform pretty much the above to get all of our values
            let quat = self.ori.index_axis(Axis(1), 0);

            #[cfg(feature = "parallel")]
            par_azip!((mut vec in vec.axis_iter_mut(Axis(1))) {
                let mut rvec = Array1::<f64>::zeros((3).f()); 
                quat_rot_vec(&quat, &vec.view(), rvec.view_mut());
                vec.assign(&rvec);  
            });

            #[cfg(not(feature = "parallel"))]
            azip!((mut vec in vec.axis_iter_mut(Axis(1))) {
                let mut rvec = Array1::<f64>::zeros((3).f()); 
                quat_rot_vec(&quat, &vec.view(), rvec.view_mut());
                vec.assign(&rvec);  
            });
        }//End of if-else
    }//End of rot_vector_inplace
}//Endo of Impl RotVector

//A helper function for Impl RotVector for Quat

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