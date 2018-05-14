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
///A structure that holds an array of Rodrigues vectors
#[derive(Clone, Debug)]
pub struct RodVec{
    ori: Array2<f64>,
}
///A structure that holds an array of compact Rodrigues vectors
#[derive(Clone, Debug)]
pub struct RodVecComp{
    ori: Array2<f64>,
}

impl RodVec{

    ///Creates an array of zeros for the initial Rodrigues vector parameterization when data is not fed into it
    pub fn new(size: usize) -> RodVec{
        assert!(size > 0, "Size inputted: {}, was not greater than 0", size);

        let mut ori = Array2::<f64>::zeros((4, size).f());

        azip!(mut rodvec (ori.axis_iter_mut(Axis(1))) in {rodvec[2] = 1.0_f64});

        RodVec{
            ori,
        }
    }//End of new

    ///Creates a Rodrigues vector parameterization  type with the supplied data as long as the supplied data is in the following format
    ///shape (4, nelems), memory order = fortran/column major.
    ///If it doesn't fit those standards it will fail.
    pub fn new_init(ori: Array2<f64>) -> RodVec{

        let nrow = ori.rows();

        assert!(nrow == 4, "Number of rows of array was: {}, which is not equal to 4", nrow);
        //We need to deal with a borrowing of ori here, so we need to have strides dropped at one point.
        {
            let strides = ori.strides();

            assert!(strides[0] == 1, "The memory stride is not column major (f order)");
        }

        RodVec{
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
}//End of Impl of RodVec

///The orientation conversions of a series of Rodrigues vectors to a number of varying different orientation
///representations commonly used in material orientation processing. 
impl OriConv for RodVec{
    ///Converts the Rodrigues vector representation over to Bunge angles which has the following properties
    ///shape (3, nelems), memory order = fortran/column major.
    fn to_bunge(&self) -> Bunge{
        let rmat = self.to_rmat();
        //When a pure conversion doesn't exist we just use the already existing ones in other orientation
        //representations 
        rmat.to_bunge()   
    }//End of to_bunge


    ///Converts the Rodrigues vector representation over to a rotation matrix which has the following properties
    ///shape (3, 3, nelems), memory order = fortran/column major.
    fn to_rmat(&self) -> RMat{
        let ang_axis = self.to_ang_axis();
        //We could convert this to a pure converesion if we wanted to save on memory usage later on
        ang_axis.to_rmat()
    }//End of to_rmat

    ///Converts the Rodrigues vector representation over to axis-angle representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    fn to_ang_axis(&self) -> AngAxis{

        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((4, nelems).f());

        azip!(mut angaxis (ori.axis_iter_mut(Axis(1))), ref rodvec (self.ori.axis_iter(Axis(1))) in {
            angaxis[0] = rodvec[0];
            angaxis[1] = rodvec[1];
            angaxis[2] = rodvec[2];
            angaxis[3] = 2.0_f64 * rodvec[3].atan();
        });

        AngAxis::new_init(ori)
    }//End of to_ang_axis

    ///Converts the Rodrigues vector representation over to a compact axial vector representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    fn to_ang_axis_comp(&self) -> AngAxisComp{
        let ang_axis = self.to_ang_axis();
        ang_axis.to_ang_axis_comp()
    }//End of to_ang_axis_comp

    ///Returns a clone of the Rodrigues vector.
    fn to_rod_vec(&self) -> RodVec{
        self.clone()
    }//End of to_rod_vec

    ///Converts the Rodrigues vector representation over to a compact Rodrigues which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    fn to_rod_vec_comp(&self) -> RodVecComp{
        //We first convert to a Rodrigues vector representation. Then we scale our normal vector by our the rotation
        //angle which is the fourth component of our axis-angle vector.
        //If we want to be more efficient about this in the future with out as many copies used we can reuse a lot of the code
        //used in the to_ang_axis code. However, we will end up with a lot of similar/repeated code then. We could put that
        //code in a helper function that isn't seen.

        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((3, nelems).f());

        azip!(mut rodvec_comp (ori.axis_iter_mut(Axis(1))), ref rodvec (self.ori.axis_iter(Axis(1))) in {
            rodvec_comp[0] = rodvec[0] * rodvec[3];
            rodvec_comp[1] = rodvec[1] * rodvec[3];
            rodvec_comp[2] = rodvec[2] * rodvec[3];
        });

        RodVecComp::new_init(ori)
    }//End of to_rod_vec_comp

    ///Converts the Rodrigues vector representation over to a unit quaternion which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    fn to_quat(&self) -> Quat{
        //Will replace this with a more direct conversion later on
        let ang_axis = self.to_ang_axis();
        ang_axis.to_quat()
    }//End of to_quat

    ///Converts the Rodrigues vector representation over to a homochoric representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    fn to_homochoric(&self) -> Homochoric{
        let ang_axis = self.to_ang_axis();
        ang_axis.to_homochoric()
    }//End of to_homochoric
    ///This operation is done inplace and does not create a new structure
    fn to_bunge_inplace(&self, bunge: &mut Bunge){
        let rmat = self.to_rmat();
        rmat.to_bunge_inplace(bunge);
    }
    ///This operation is done inplace and does not create a new structure
    fn to_rmat_inplace(&self, rmat: &mut RMat){
        let ang_axis = self.to_ang_axis();
        ang_axis.to_rmat_inplace(rmat);
    }
    ///This operation is done inplace and does not create a new structure
    fn to_ang_axis_inplace(&self, ang_axis: &mut AngAxis){
        let mut ori = ang_axis.ori_view_mut();

        let new_nelem = ori.len_of(Axis(1));
        let nelem = self.ori.len_of(Axis(1));

        assert!(new_nelem == nelem, 
        "The number of elements in the original ori field do no match up with the new field.
        The old field had {} elements, and the new field has {} elements",
        nelem, new_nelem);

        azip!(mut angaxis (ori.axis_iter_mut(Axis(1))), ref rodvec (self.ori.axis_iter(Axis(1))) in {
            angaxis[0] = rodvec[0];
            angaxis[1] = rodvec[1];
            angaxis[2] = rodvec[2];
            angaxis[3] = 2.0_f64 * rodvec[3].atan();
        });
    }
    ///This operation is done inplace and does not create a new structure
    fn to_ang_axis_comp_inplace(&self, ang_axis_comp: &mut AngAxisComp){
        let ang_axis = self.to_ang_axis();
        ang_axis.to_ang_axis_comp_inplace(ang_axis_comp);
    }
    ///This operation is done inplace and does not create a new structure
    fn to_rod_vec_inplace(&self, rod_vec: &mut RodVec){
        let mut ori = rod_vec.ori_view_mut();

        let new_nelem = ori.len_of(Axis(1));
        let nelem = self.ori.len_of(Axis(1));

        assert!(new_nelem == nelem, 
        "The number of elements in the original ori field do no match up with the new field.
        The old field had {} elements, and the new field has {} elements",
        nelem, new_nelem);

        ori.assign(&self.ori);
    }
    ///This operation is done inplace and does not create a new structure
    fn to_rod_vec_comp_inplace(&self, rod_vec_comp: &mut RodVecComp){
        let mut ori = rod_vec_comp.ori_view_mut();

        let new_nelem = ori.len_of(Axis(1));
        let nelem = self.ori.len_of(Axis(1));

        assert!(new_nelem == nelem, 
        "The number of elements in the original ori field do no match up with the new field.
        The old field had {} elements, and the new field has {} elements",
        nelem, new_nelem);

        azip!(mut rodvec_comp (ori.axis_iter_mut(Axis(1))), ref rodvec (self.ori.axis_iter(Axis(1))) in {
            rodvec_comp[0] = rodvec[0] * rodvec[3];
            rodvec_comp[1] = rodvec[1] * rodvec[3];
            rodvec_comp[2] = rodvec[2] * rodvec[3];
        });

    }
    ///This operation is done inplace and does not create a new structure
    fn to_quat_inplace(&self, quat: &mut Quat){
        let ang_axis = self.to_ang_axis();
        ang_axis.to_quat_inplace(quat);
    }
    ///This operation is done inplace and does not create a new structure
    fn to_homochoric_inplace(&self, homochoric: &mut Homochoric){
        let ang_axis = self.to_ang_axis();
        ang_axis.to_homochoric_inplace(homochoric);
    }


}//End of impl Ori_Conv of Rodrigues Vector

///A series of commonly used operations to rotate vector data by a given rotation
impl RotVector for RodVec{

    ///rot_vector takes in a 2D array view of a series of vectors. It then rotates these vectors using the
    ///given Rodrigues vectors. The newly rotated vectors are then returned. This function requires the
    ///number of elements in the Rodrigues vectors to be either 1 or nelems.
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
        Rodrigues vector structure, or their must only be one element in Rodrigues vector. The final case is
        that there must only be one element in the vector field. There are
        currently {} elements in vector and {} elements in Rodrigues vector",
        nelems, rnelems);

        let mnelems = cmp::max(rnelems, nelems);
        let mut rvec = Array2::<f64>::zeros((3, mnelems).f());

        //We need to see if we have more than one rotation matrix that we're multiplying by
        if rnelems == nelems {
            //The rotations here can be given by the following set of equations as found on Wikipedia:
            //https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula#Statement
            azip!(mut rvec (rvec.axis_iter_mut(Axis(1))), ref vec (vec.axis_iter(Axis(1))), 
            ref rod_vec (self.ori.axis_iter(Axis(1))) in {
                rod_vec_rot_vec(&rod_vec, &vec, rvec);     
            });
        } else if rnelems == 1{
            //We just have one Rodrigues vector so perform pretty much the above to get all of our values
            let rod_vec = self.ori.subview(Axis(1), 0);

            azip!(mut rvec (rvec.axis_iter_mut(Axis(1))), ref vec (vec.axis_iter(Axis(1))) in {  
                rod_vec_rot_vec(&rod_vec, &vec, rvec); 
            });
        }else{
            //We just have one Rodrigues vector so perform pretty much the above to get all of our values
            let vec = vec.subview(Axis(1), 0);

            azip!(mut rvec (rvec.axis_iter_mut(Axis(1))), ref rod_vec (self.ori.axis_iter(Axis(1))) in {  
                rod_vec_rot_vec(&rod_vec, &vec, rvec); 
            }); 
        }//End if-else
        //Now we just need to return the rvec value
        rvec
    }//End of rot_vector

    ///rot_vector_mut takes in a 2D array view of a series of vectors and a mutable 2D ArrayView of the 
    ///rotated vector. It then rotates these vectors using the given Rodrigues vector. The newly rotated
    /// vectors are assigned to the supplied rotated vector, rvec. This function requires the
    ///number of elements in the Rodrigues vector to be either 1 or nelems.
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
        "The number of elements in the unrotated vector or Rodrigues vector field must be equal to the number of elements
        in the supplied rotated vector field. There are currently {} elements in the unrotated vector or Rodrigues vector
        field and {} elements in the rotated vector field", 
        mnelems, rvnelems);

        assert!( (nelems == rnelems) | (rnelems == 1) | (nelems == 1), 
        "The number of elements in the vector field must be equal to the number of elements in the
        Rodrigues vector structure, or their must only be one element in Rodrigues vector. The final case is
        that there must only be one element in the vector field. There are
        currently {} elements in vector and {} elements in Rodrigues vector",
        nelems, rnelems);

        //We need to see if we have more than one Rodrigues vector that we're multiplying by
        if rnelems == nelems {
            //The rotations here can be given by the following set of equations as found on Wikipedia:
            //https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula#Statement
            azip!(mut rvec (rvec.axis_iter_mut(Axis(1))), ref vec (vec.axis_iter(Axis(1))), 
            ref rod_vec (self.ori.axis_iter(Axis(1))) in {
                rod_vec_rot_vec(&rod_vec, &vec, rvec);    
            });
        } else if rnelems == 1{
            //We just have one Rodrigues vector so perform pretty much the above to get all of our values
            let rod_vec = self.ori.subview(Axis(1), 0);

            azip!(mut rvec (rvec.axis_iter_mut(Axis(1))), ref vec (vec.axis_iter(Axis(1))) in {  
                rod_vec_rot_vec(&rod_vec, &vec, rvec); 
            });
        }else{
            //We just have one Rodrigues vector so perform pretty much the above to get all of our values
            let vec = vec.subview(Axis(1), 0);

            azip!(mut rvec (rvec.axis_iter_mut(Axis(1))), ref rod_vec (self.ori.axis_iter(Axis(1))) in {  
                rod_vec_rot_vec(&rod_vec, &vec, rvec); 
            }); 
        }//End if-else
    }//End of rot_vector_mut

    ///rot_vector_inplace takes in a mutable 2D array view of a series of vectors. It then rotates these vectors using the
    ///given Rodrigues vectors. The newly rotated vectors are assigned to original vector. This function requires the
    ///number of elements in the Rodrigues vector to be either 1 or nelems where vec has nelems in it.
    ///If this condition is not met the function will error out.
    ///vec - the vector to be rotated must have dimensions 3xnelems.
    fn rot_vector_inplace(&self, mut vec: ArrayViewMut2<f64>){

        let nelems = vec.len_of(Axis(1));
        let rnelems = self.ori.len_of(Axis(1));

        let rows  = vec.len_of(Axis(0));
        assert!((rows == 3), "The number of rows must be 3. The number of rows provided is {}", rows); 

        assert!( (nelems == rnelems) | (rnelems == 1), 
        "The number of elements in the vector field must be equal to the number of elements in the
        Rodrigues vector structure, or their must only be one element in Rodrigues vector. There are
        currently {} elements in vector and {} elements in Rodrigues vector",
        nelems, rnelems);

        //We need to see if we have more than one Rodrigues vector that we're multiplying by
        if rnelems == nelems {
            //The rotations here can be given by the following set of equations as found on Wikipedia:
            //https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula#Statement
            azip!(mut vec (vec.axis_iter_mut(Axis(1))), ref rod_vec (self.ori.axis_iter(Axis(1))) in {
                let mut rvec = Array1::<f64>::zeros((3).f());
                rod_vec_rot_vec(&rod_vec, &vec.view(), rvec.view_mut());
                vec.assign({&rvec});    
            });
        } else{
            //We just have one Rodrigues vector so perform pretty much the above to get all of our values
            let rod_vec = self.ori.subview(Axis(1), 0);

            azip!(mut vec (vec.axis_iter_mut(Axis(1))) in {
                let mut rvec = Array1::<f64>::zeros((3).f());
                rod_vec_rot_vec(&rod_vec, &vec.view(), rvec.view_mut());
                vec.assign({&rvec});  
            });
        }//End if-else
    }//End of rot_vector_inplace
}//Endo of Impl RotVector

///All of the Rodrigues vector vector rotation operations can be described by using the below series of functions.
///This also reduces the amount of repetive code that existed earlier within rot_vector. 
fn rod_vec_rot_vec(rod_vec: &ArrayView1<f64>, vec: &ArrayView1<f64>, mut rvec: ArrayViewMut1<f64>) {
    let (sin_theta, cos_theta) = f64::sin_cos(2.0_f64 * rod_vec[3].atan());
    let min_cos = 1.0_f64 - cos_theta;
    let dp_mcos = min_cos * (rod_vec[0] * vec[0] + rod_vec[1] * vec[1] + rod_vec[2] * vec[2]);
    let mut cross_prod = Array1::<f64>::zeros((3).f());

    cross_prod[0] = -rod_vec[2] * vec[1] + rod_vec[1] * vec[2];
    cross_prod[1] = rod_vec[2] * vec[0] - rod_vec[0] * vec[1];
    cross_prod[2] = -rod_vec[1] * vec[0] + rod_vec[0] * vec[1];

    rvec[0] = vec[0] * cos_theta + cross_prod[0] * sin_theta + rod_vec[0] * dp_mcos;
    rvec[1] = vec[1] * cos_theta + cross_prod[1] * sin_theta + rod_vec[1] * dp_mcos;
    rvec[2] = vec[2] * cos_theta + cross_prod[2] * sin_theta + rod_vec[2] * dp_mcos;
}

impl RodVecComp{

    ///Creates an array of zeros for the initial compact Rodrigues vector parameterization when data is not fed into it
    pub fn new(size: usize) -> RodVecComp{
        assert!(size > 0, "Size inputted: {}, was not greater than 0", size);

        let ori = Array2::<f64>::zeros((3, size).f());

        RodVecComp{
            ori,
        }
    }//End of new

    ///Creates a compact Rodrigues vector parameterization type with the supplied data as long as the supplied data is in the following format
    ///shape (4, nelems), memory order = fortran/column major.
    ///If it doesn't fit those standards it will fail.
    pub fn new_init(ori: Array2<f64>) -> RodVecComp{

        let nrow = ori.rows();

        assert!(nrow == 3, "Number of rows of array was: {}, which is not equal to 4", nrow);
        //We need to deal with a borrowing of ori here, so we need to have strides dropped at one point.
        {
            let strides = ori.strides();

            assert!(strides[0] == 1, "The memory stride is not column major (f order)");
        }

        RodVecComp{
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
}//End of Impl of RodVecComp

///The orientation conversions of a series of Rodrigues vectors to a number of varying different orientation
///representations commonly used in material orientation processing. 
impl OriConv for RodVecComp{
    ///Converts the Rodrigues vector representation over to Bunge angles which has the following properties
    ///shape (3, nelems), memory order = fortran/column major.
    fn to_bunge(&self) -> Bunge{
        let rmat = self.to_rmat();
        //When a pure conversion doesn't exist we just use the already existing ones in other orientation
        //representations 
        rmat.to_bunge()   
    }//End of to_bunge


    ///Converts the Rodrigues vector representation over to a rotation matrix which has the following properties
    ///shape (3, 3, nelems), memory order = fortran/column major.
    fn to_rmat(&self) -> RMat{
        let ang_axis = self.to_ang_axis();
        //We could convert this to a pure converesion if we wanted to save on memory usage later on
        ang_axis.to_rmat()
    }//End of to_rmat

    ///Converts the Rodrigues vector representation over to axis-angle representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    fn to_ang_axis(&self) -> AngAxis{
        let rod_vec = self.to_rod_vec();
        rod_vec.to_ang_axis()
    }//End of to_ang_axis

    ///Converts the Rodrigues vector representation over to a compact axial vector representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    fn to_ang_axis_comp(&self) -> AngAxisComp{
        let ang_axis = self.to_ang_axis();
        ang_axis.to_ang_axis_comp()
    }//End of to_ang_axis_comp

    ///Converts the compact Rodrigues vector representation over to a Rodrigues which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    fn to_rod_vec(&self) -> RodVec{

        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((4, nelems).f());

        let tol = std::f64::EPSILON;

        azip!(mut rodvec (ori.axis_iter_mut(Axis(1))), ref rodvec_comp (self.ori.axis_iter(Axis(1))) in {
            let norm_rodvec = f64::sqrt({
                rodvec_comp[0] * rodvec_comp[0] 
                + rodvec_comp[1] * rodvec_comp[1] 
                + rodvec_comp[2] * rodvec_comp[2]
                });
            //If we follow the same convention that we use with quaternions for cases with no rotation
            //then we set it equal to the following vector with the no rotation ([0, 0, 1], 0)
            if norm_rodvec.abs() < tol {
                rodvec[2] = 1.0_f64;
            }else{
                let inv_norm_rodvec = 1.0_f64 / norm_rodvec;
                rodvec[0] = rodvec_comp[0] * inv_norm_rodvec;
                rodvec[1] = rodvec_comp[1] * inv_norm_rodvec;
                rodvec[2] = rodvec_comp[2] * inv_norm_rodvec;
                rodvec[3] = norm_rodvec;
            }
        });

        RodVec::new_init(ori)
    }//End of to_rod_vec

    ///Returns a clone of the compact Rodrigues vector structure
    fn to_rod_vec_comp(&self) -> RodVecComp{
        self.clone()
    }//End of to_rod_vec_comp

    ///Converts the compact Rodrigues vector representation over to a unit quaternion which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    fn to_quat(&self) -> Quat{
        //Will replace this with a more direct conversion later on
        let ang_axis = self.to_ang_axis();
        ang_axis.to_quat()
    }//End of to_quat

    ///Converts the compact Rodrigues vector representation over to a homochoric representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    fn to_homochoric(&self) -> Homochoric{
        let ang_axis = self.to_ang_axis();
        ang_axis.to_homochoric()
    }//End of to_homochoric
    ///This operation is done inplace and does not create a new structure
    fn to_bunge_inplace(&self, bunge: &mut Bunge){
        let rmat = self.to_rmat();
        rmat.to_bunge_inplace(bunge);
    }
    ///This operation is done inplace and does not create a new structure
    fn to_rmat_inplace(&self, rmat: &mut RMat){
        let ang_axis = self.to_ang_axis();
        ang_axis.to_rmat_inplace(rmat);
    }
    ///This operation is done inplace and does not create a new structure
    fn to_ang_axis_inplace(&self, ang_axis: &mut AngAxis){
        let rod_vec = self.to_rod_vec();
        rod_vec.to_ang_axis_inplace(ang_axis);
    }
    ///This operation is done inplace and does not create a new structure
    fn to_ang_axis_comp_inplace(&self, ang_axis_comp: &mut AngAxisComp){
        let ang_axis = self.to_ang_axis();
        ang_axis.to_ang_axis_comp_inplace(ang_axis_comp);
    }
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

        azip!(mut rodvec (ori.axis_iter_mut(Axis(1))), ref rodvec_comp (self.ori.axis_iter(Axis(1))) in {
            let norm_rodvec = f64::sqrt({
                rodvec_comp[0] * rodvec_comp[0] 
                + rodvec_comp[1] * rodvec_comp[1] 
                + rodvec_comp[2] * rodvec_comp[2]
                });
            //If we follow the same convention that we use with quaternions for cases with no rotation
            //then we set it equal to the following vector with the no rotation ([0, 0, 1], 0)
            if norm_rodvec.abs() < tol {
                rodvec[2] = 1.0_f64;
            }else{
                let inv_norm_rodvec = 1.0_f64 / norm_rodvec;
                rodvec[0] = rodvec_comp[0] * inv_norm_rodvec;
                rodvec[1] = rodvec_comp[1] * inv_norm_rodvec;
                rodvec[2] = rodvec_comp[2] * inv_norm_rodvec;
                rodvec[3] = norm_rodvec;
            }
        });
    }
    ///This operation is done inplace and does not create a new structure
    fn to_rod_vec_comp_inplace(&self, rod_vec_comp: &mut RodVecComp){
        let mut ori = rod_vec_comp.ori_view_mut();

        let new_nelem = ori.len_of(Axis(1));
        let nelem = self.ori.len_of(Axis(1));

        assert!(new_nelem == nelem, 
        "The number of elements in the original ori field do no match up with the new field.
        The old field had {} elements, and the new field has {} elements",
        nelem, new_nelem);

        ori.assign(&self.ori);
    }
    ///This operation is done inplace and does not create a new structure
    fn to_quat_inplace(&self, quat: &mut Quat){
        let ang_axis = self.to_ang_axis();
        ang_axis.to_quat_inplace(quat);
    }
    ///This operation is done inplace and does not create a new structure
    fn to_homochoric_inplace(&self, homochoric: &mut Homochoric){
        let ang_axis = self.to_ang_axis();
        ang_axis.to_homochoric_inplace(homochoric);
    }


}//End of impl Ori_Conv of Compact Rodrigues Vector 