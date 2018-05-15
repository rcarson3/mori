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

///A structure that holds an array of axis-angle representation of a rotation
#[derive(Clone, Debug)]
pub struct AngAxis{
    ori: Array2<f64>,
}


impl AngAxis{

    ///Creates an array of zeros for the initial axis-angle parameterization when data is not fed into it
    pub fn new(size: usize) -> AngAxis{
        assert!(size > 0, "Size inputted: {}, was not greater than 0", size);

        let mut ori = Array2::<f64>::zeros((4, size).f());

        azip!(mut angaxis (ori.axis_iter_mut(Axis(1))) in {angaxis[2] = 1.0_f64});

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

    ///Converts the axis-angle representation over to Bunge angles which has the following properties
    ///shape (3, nelems), memory order = fortran/column major.
    fn to_bunge(&self) -> Bunge{

        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((3, nelems).f());

        //We need to check the R_33 component to see if it's near 1.0 
        let tol = f64::sqrt(std::f64::EPSILON);

        azip!(mut bunge (ori.axis_iter_mut(Axis(1))), ref angaxis (self.ori.axis_iter(Axis(1))) in {
            let c = angaxis[3].cos();
            let s = angaxis[3].sin();

            let mut rmat = Array2::<f64>::zeros((3, 3).f());

            rmat[[0, 0]] = c + (1.0_f64 - c) * (angaxis[0] * angaxis[0]);
            rmat[[1, 0]] = (1.0_f64 - c) * (angaxis[0] * angaxis[1]) + s * angaxis[2];
            rmat[[2, 0]] = (1.0_f64 - c) * (angaxis[0] * angaxis[2]) - s * angaxis[1];

            rmat[[0, 1]] = (1.0_f64 - c) * (angaxis[0] * angaxis[1]) - s * angaxis[2];
            rmat[[1, 1]] = c + (1.0_f64 - c) * (angaxis[1] * angaxis[1]);
            rmat[[2, 1]] = (1.0_f64 - c) * (angaxis[1] * angaxis[2]) + s * angaxis[0];

            rmat[[0, 2]] = (1.0_f64 - c) * (angaxis[0] * angaxis[2]) + s * angaxis[1];
            rmat[[1, 2]] = (1.0_f64 - c) * (angaxis[1] * angaxis[2]) - s * angaxis[0];
            rmat[[2, 2]] = c + (1.0_f64 - c) * (angaxis[2] * angaxis[2]);

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

    }//End of to_bunge

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

        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((3, nelems).f());

        let inv2 = 1.0_f64/2.0_f64;
        azip!(mut rodvec (ori.axis_iter_mut(Axis(1))), ref angaxis (self.ori.axis_iter(Axis(1))) in {
            let tan2 = f64::tan(inv2 * angaxis[3]);
            rodvec[0] = angaxis[0] * tan2;
            rodvec[1] = angaxis[1] * tan2;
            rodvec[2] = angaxis[2] * tan2;
        });

        RodVecComp::new_init(ori)
    }//End of to_rod_vec_comp

    ///Converts the axis-angle representation over to a unit quaternion representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    fn to_quat(&self) -> Quat{

        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((4, nelems).f());

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

    ///Converts the axis-angle representation over to Bunge angles which has the following properties
    ///shape (3, nelems), memory order = fortran/column major.
    ///This operation is done inplace and does not create a new structure
    fn to_bunge_inplace(&self, bunge: &mut Bunge){
        let mut ori = bunge.ori_view_mut();

        let new_nelem = ori.len_of(Axis(2));
        let nelem = self.ori.len_of(Axis(1));

        assert!(new_nelem == nelem, 
        "The number of elements in the original ori field do no match up with the new field.
        The old field had {} elements, and the new field has {} elements",
        nelem, new_nelem);

        //We need to check the R_33 component to see if it's near 1.0 
        let tol = f64::sqrt(std::f64::EPSILON);

        azip!(mut bunge (ori.axis_iter_mut(Axis(1))), ref angaxis (self.ori.axis_iter(Axis(1))) in {
            let c = angaxis[3].cos();
            let s = angaxis[3].sin();

            let mut rmat = Array2::<f64>::zeros((3, 3).f());

            rmat[[0, 0]] = c + (1.0_f64 - c) * (angaxis[0] * angaxis[0]);
            rmat[[1, 0]] = (1.0_f64 - c) * (angaxis[0] * angaxis[1]) + s * angaxis[2];
            rmat[[2, 0]] = (1.0_f64 - c) * (angaxis[0] * angaxis[2]) - s * angaxis[1];

            rmat[[0, 1]] = (1.0_f64 - c) * (angaxis[0] * angaxis[1]) - s * angaxis[2];
            rmat[[1, 1]] = c + (1.0_f64 - c) * (angaxis[1] * angaxis[1]);
            rmat[[2, 1]] = (1.0_f64 - c) * (angaxis[1] * angaxis[2]) + s * angaxis[0];

            rmat[[0, 2]] = (1.0_f64 - c) * (angaxis[0] * angaxis[2]) + s * angaxis[1];
            rmat[[1, 2]] = (1.0_f64 - c) * (angaxis[1] * angaxis[2]) - s * angaxis[0];
            rmat[[2, 2]] = c + (1.0_f64 - c) * (angaxis[2] * angaxis[2]);

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

    ///Converts the axis-angle representation over to a rotation matrix which has the following properties
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
    }

    ///Copies the orientation structure in self over to the other AngAxis structure 
    ///This operation is done inplace and does not create a new structure
    fn to_ang_axis_inplace(&self, ang_axis: &mut AngAxis){
        let mut ori = ang_axis.ori_view_mut();

        let new_nelem = ori.len_of(Axis(1));
        let nelem = self.ori.len_of(Axis(1));

        assert!(new_nelem == nelem, 
        "The number of elements in the original ori field do no match up with the new field.
        The old field had {} elements, and the new field has {} elements",
        nelem, new_nelem);

        ori.assign(&self.ori);

    }

    ///Converts the axis-angle representation over to a compact angle-axis representation which has the following properties
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

        azip!(mut angaxis_comp (ori.axis_iter_mut(Axis(1))), ref angaxis (self.ori.axis_iter(Axis(1))) in {
            angaxis_comp[0] = angaxis[0] * angaxis[3];
            angaxis_comp[1] = angaxis[1] * angaxis[3];
            angaxis_comp[2] = angaxis[2] * angaxis[3];
        });

    }

    ///Converts the axis-angle representation over to a Rodrigues vector representation which has the following properties
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

        let inv2 = 1.0_f64/2.0_f64;

        azip!(mut rodvec (ori.axis_iter_mut(Axis(1))), ref angaxis (self.ori.axis_iter(Axis(1))) in {
            rodvec[0] = angaxis[0];
            rodvec[1] = angaxis[1];
            rodvec[2] = angaxis[2];
            rodvec[3] = f64::tan(inv2 * angaxis[3]);
        });
    }
    
    ///Converts the axis-angle representation over to a compact Rodrigues vector representation which has the following properties
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

        let inv2 = 1.0_f64/2.0_f64;

        azip!(mut rodvec (ori.axis_iter_mut(Axis(1))), ref angaxis (self.ori.axis_iter(Axis(1))) in {
            let tan2 = f64::tan(inv2 * angaxis[3]);
            rodvec[0] = angaxis[0] * tan2;
            rodvec[1] = angaxis[1] * tan2;
            rodvec[2] = angaxis[2] * tan2;
        });
    }

    ///Converts the axis-angle representation over to a unit quaternion representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    ///This operation is done inplace and does not create a new structure
    fn to_quat_inplace(&self, quat: &mut Quat){
        let mut ori = quat.ori_view_mut();

        let new_nelem = ori.len_of(Axis(1));
        let nelem = self.ori.len_of(Axis(1));

        assert!(new_nelem == nelem, 
        "The number of elements in the original ori field do no match up with the new field.
        The old field had {} elements, and the new field has {} elements",
        nelem, new_nelem);

        let inv2 = 1.0_f64 / 2.0_f64;

        azip!(mut quat (ori.axis_iter_mut(Axis(1))), ref angaxis (self.ori.axis_iter(Axis(1))) in {
            let s = f64::sin(inv2 * angaxis[3]); 

            quat[0] = f64::cos(inv2 * angaxis[3]);
            quat[1] = s * angaxis[0];
            quat[2] = s * angaxis[1];
            quat[3] = s * angaxis[2];
        });
    }

    ///Converts the axis-angle representation over to a homochoric representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    ///This operation is done inplace and does not create a new structure
    fn to_homochoric_inplace(&self, homochoric: &mut Homochoric){
        let mut ori = homochoric.ori.view_mut();

        let new_nelem = ori.len_of(Axis(1));
        let nelem = self.ori.len_of(Axis(1));

        assert!(new_nelem == nelem, 
        "The number of elements in the original ori field do no match up with the new field.
        The old field had {} elements, and the new field has {} elements",
        nelem, new_nelem);

        let inv3  = 1.0_f64 / 3.0_f64;
        let inv34 = 3.0_f64 / 4.0_f64;

        azip!(mut homoch (ori.axis_iter_mut(Axis(1))), ref angaxis (self.ori.axis_iter(Axis(1))) in {
            let pow_term    = inv34 * (angaxis[3] - angaxis[3].sin()); 

            homoch[0] = angaxis[0];
            homoch[1] = angaxis[1];
            homoch[2] = angaxis[2];
            homoch[3] = pow_term.powf(inv3);
        });
    }


}//End of Impl OriConv for AngAxis

///A series of commonly used operations to rotate vector data by a given rotation
impl RotVector for AngAxis{

    ///rot_vector takes in a 2D array view of a series of vectors. It then rotates these vectors using the
    ///given Axis-angle representation. The newly rotated vectors are then returned. This function requires the
    ///number of elements in the Axis-angle representation to be either 1 or nelems.
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
        Axis-angle representation structure, or their must only be one element in Axis-angle representation. The final case is
        that there must only be one element in the vector field. There are
        currently {} elements in vector and {} elements in Axis-angle representation",
        nelems, rnelems);

        let mnelems = cmp::max(rnelems, nelems);
        let mut rvec = Array2::<f64>::zeros((3, mnelems).f());

        //We need to see if we have more than one Axis-angle representation that we're multiplying by
        if rnelems == nelems {
            //The rotations here can be given by the following set of equations as found on Wikipedia:
            //https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula#Statement
            azip!(mut rvec (rvec.axis_iter_mut(Axis(1))), ref vec (vec.axis_iter(Axis(1))), 
            ref ang_axis (self.ori.axis_iter(Axis(1))) in {
                ang_axis_rot_vec(&ang_axis, &vec, rvec);    
            });
        } else if rnelems == 1{
            //We just have one Axis-angle representation so perform pretty much the above to get all of our values
            let ang_axis = self.ori.subview(Axis(1), 0);

            azip!(mut rvec (rvec.axis_iter_mut(Axis(1))), ref vec (vec.axis_iter(Axis(1))) in {  
                ang_axis_rot_vec(&ang_axis, &vec, rvec);   
            });
        } else {
            //We just have one vector so perform pretty much the above to get all of our values
            let vec = vec.subview(Axis(1), 0);

            azip!(mut rvec (rvec.axis_iter_mut(Axis(1))), ref ang_axis (self.ori.axis_iter(Axis(1))) in {  
                ang_axis_rot_vec(&ang_axis, &vec, rvec);  
            });
        }//End if-else
        //Now we just need to return the rvec value
        rvec
    }//End of rot_vector

    ///rot_vector_mut takes in a 2D array view of a series of vectors and a mutable 2D ArrayView of the 
    ///rotated vector. It then rotates these vectors using the given Axis-angle representation. The newly rotated
    /// vectors are assigned to the supplied rotated vector, rvec. This function requires the
    ///number of elements in the Axis-angle representation to be either 1 or nelems.
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
        "The number of elements in the unrotated vector or axis-angle field must be equal to the number of elements
        in the supplied rotated vector field. There are currently {} elements in the unrotated vector or axis-angle
        field and {} elements in the rotated vector field", 
        mnelems, rvnelems);

        assert!( (nelems == rnelems) | (rnelems == 1) | (nelems == 1), 
        "The number of elements in the vector field must be equal to the number of elements in the
        Axis-angle representation structure, or their must only be one element in Axis-angle representation. The final case is
        that there must only be one element in the vector field. There are
        currently {} elements in vector and {} elements in Axis-angle representation",
        nelems, rnelems);

        //We need to see if we have more than one Axis-angle representation that we're multiplying by
        if rnelems == nelems {
            //The rotations here can be given by the following set of equations as found on Wikipedia:
            //https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula#Statement
            azip!(mut rvec (rvec.axis_iter_mut(Axis(1))), ref vec (vec.axis_iter(Axis(1))), 
            ref ang_axis (self.ori.axis_iter(Axis(1))) in {
                ang_axis_rot_vec(&ang_axis, &vec, rvec);   
            });
        } else if rnelems == 1{
            //We just have one Axis-angle representation so perform pretty much the above to get all of our values
            let ang_axis = self.ori.subview(Axis(1), 0);

            azip!(mut rvec (rvec.axis_iter_mut(Axis(1))), ref vec (vec.axis_iter(Axis(1))) in {  
                ang_axis_rot_vec(&ang_axis, &vec, rvec);
            });
        } else{
            //We just have one vector so perform pretty much the above to get all of our values
            let vec = vec.subview(Axis(1), 0);

            azip!(mut rvec (rvec.axis_iter_mut(Axis(1))), ref ang_axis (self.ori.axis_iter(Axis(1))) in {  
                ang_axis_rot_vec(&ang_axis, &vec, rvec);  
            });
        }//End of if-else
    }//End of rot_vector_mut

    ///rot_vector_inplace takes in a mutable 2D array view of a series of vectors. It then rotates these vectors using the
    ///given Axis-angle representation. The newly rotated vectors are assigned to original vector. This function requires the
    ///number of elements in the Axis-angle representation to be either 1 or nelems where vec has nelems in it.
    ///If this condition is not met the function will error out.
    ///vec - the vector to be rotated must have dimensions 3xnelems.
    fn rot_vector_inplace(&self, mut vec: ArrayViewMut2<f64>){

        let nelems = vec.len_of(Axis(1));
        let rnelems = self.ori.len_of(Axis(1));

        let rows  = vec.len_of(Axis(0));
        assert!((rows == 3), "The number of rows must be 3. The number of rows provided is {}", rows); 

        assert!( (nelems == rnelems) | (rnelems == 1), 
        "The number of elements in the vector field must be equal to the number of elements in the
        Axis-angle representation structure, or their must only be one element in Axis-angle representation. There are
        currently {} elements in vector and {} elements in Axis-angle representation",
        nelems, rnelems);

        //We need to see if we have more than one Axis-angle representation that we're multiplying by
        if rnelems == nelems {
            //The rotations here can be given by the following set of equations as found on Wikipedia:
            //https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula#Statement
            azip!(mut vec (vec.axis_iter_mut(Axis(1))), ref ang_axis (self.ori.axis_iter(Axis(1))) in {
                let mut rvec = Array1::<f64>::zeros((3).f());
                ang_axis_rot_vec(&ang_axis, &vec.view(), rvec.view_mut());
                vec.assign({&rvec});    
            });
        } else{
            //We just have one Axis-angle representation so perform pretty much the above to get all of our values
            let ang_axis = self.ori.subview(Axis(1), 0);

            azip!(mut vec (vec.axis_iter_mut(Axis(1))) in {
                let mut rvec = Array1::<f64>::zeros((3).f());
                ang_axis_rot_vec(&ang_axis, &vec.view(), rvec.view_mut());
                vec.assign({&rvec});  
            });
        }//End of if-else
    }//End of rot_vector_inplace
}//Endo of Impl RotVector

///All of the axis-angle vector rotation operations can be described by using the below series of functions.
///This also reduces the amount of repetive code that existed earlier within rot_vector. 
fn ang_axis_rot_vec(ang_axis: &ArrayView1<f64>, vec: &ArrayView1<f64>, mut rvec: ArrayViewMut1<f64>) {
    let (sin_theta, cos_theta) = ang_axis[3].sin_cos();
    let min_cos = 1.0_f64 - cos_theta;
    let dp_mcos = min_cos * (ang_axis[0] * vec[0] + ang_axis[1] * vec[1] + ang_axis[2] * vec[2]);
    let mut cross_prod = Array1::<f64>::zeros((3).f());

    cross_prod[0] = -ang_axis[2] * vec[1] + ang_axis[1] * vec[2];
    cross_prod[1] = ang_axis[2] * vec[0] - ang_axis[0] * vec[1];
    cross_prod[2] = -ang_axis[1] * vec[0] + ang_axis[0] * vec[1];

    rvec[0] = vec[0] * cos_theta + cross_prod[0] * sin_theta + ang_axis[0] * dp_mcos;
    rvec[1] = vec[1] * cos_theta + cross_prod[1] * sin_theta + ang_axis[1] * dp_mcos;
    rvec[2] = vec[2] * cos_theta + cross_prod[2] * sin_theta + ang_axis[2] * dp_mcos;
}
