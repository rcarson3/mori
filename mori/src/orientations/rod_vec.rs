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


impl RodVec{

    ///Creates an array of zeros for the initial Rodrigues vector parameterization when data is not fed into it
    pub fn new(size: usize) -> RodVec{
        assert!(size > 0, "Size inputted: {}, was not greater than 0", size);

        let mut ori = Array2::<f64>::zeros((4, size).f());

        azip!((mut rod_vec in ori.axis_iter_mut(Axis(1))) {rod_vec[2] = 1.0_f64});

        RodVec{
            ori,
        }
    }//End of new

    ///Creates a Rodrigues vector parameterization  type with the supplied data as long as the supplied data is in the following format
    ///shape (4, nelems), memory order = fortran/column major.
    ///If it doesn't fit those standards it will fail.
    pub fn new_init(ori: Array2<f64>) -> RodVec{

        let nrow = ori.nrows();

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

    ///Returns a new RodVec that is equal to the equivalent of transposing a rotation matrix.
    ///It turns out this is simply the negative of the normal vector due to the vector being formed
    ///from an axial vector of the rotation matrix --> Rmat\^T = -Rx where Rx is the axial vector.
    pub fn transpose(&self) -> RodVec{
        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((4, nelems).f());

        let f = |mut rod_vec_t: ArrayViewMut1::<f64>, ref rod_vec: ArrayView1::<f64>| {
            rod_vec_t[0] = -1.0_f64 * rod_vec[0];
            rod_vec_t[1] = -1.0_f64 * rod_vec[1];
            rod_vec_t[2] = -1.0_f64 * rod_vec[2];
            rod_vec_t[3] = rod_vec[3];
        };
        
        azip!((rod_vec_t in ori.axis_iter_mut(Axis(1)), rod_vec in self.ori.axis_iter(Axis(1))) {
            f(rod_vec_t, rod_vec);
        });

        RodVec::new_init(ori)
    }

    ///Performs the equivalent of transposing a rotation matrix on the internal orientations.
    ///It turns out this is simply the negative of the normal vector due to the vector being formed
    ///from an axial vector of the rotation matrix --> Rmat\^T = -Rx where Rx is the axial vector.
    pub fn transpose_inplace(&mut self){
        let f = |mut rod_vec_t: ArrayViewMut1::<f64>| {
            rod_vec_t[0] *= -1.0_f64;
            rod_vec_t[1] *= -1.0_f64;
            rod_vec_t[2] *= -1.0_f64;
        };

        azip!((rod_vec_t in self.ori.axis_iter_mut(Axis(1))) {
            f(rod_vec_t);
        });
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

        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array3::<f64>::zeros((3, 3, nelems).f());

        let f = |mut rmat: ArrayViewMut2::<f64>, ref rod_vec: ArrayView1::<f64>| {
            let phi = rod_vec[3].atan() * 2.0_f64;
            let c = phi.cos();
            let s = phi.sin();

            rmat[[0, 0]] = c + (1.0_f64 - c) * (rod_vec[0] * rod_vec[0]);
            rmat[[1, 0]] = (1.0_f64 - c) * (rod_vec[0] * rod_vec[1]) + s * rod_vec[2];
            rmat[[2, 0]] = (1.0_f64 - c) * (rod_vec[0] * rod_vec[2]) - s * rod_vec[1];

            rmat[[0, 1]] = (1.0_f64 - c) * (rod_vec[0] * rod_vec[1]) - s * rod_vec[2];
            rmat[[1, 1]] = c + (1.0_f64 - c) * (rod_vec[1] * rod_vec[1]);
            rmat[[2, 1]] = (1.0_f64 - c) * (rod_vec[1] * rod_vec[2]) + s * rod_vec[0];

            rmat[[0, 2]] = (1.0_f64 - c) * (rod_vec[0] * rod_vec[2]) + s * rod_vec[1];
            rmat[[1, 2]] = (1.0_f64 - c) * (rod_vec[1] * rod_vec[2]) - s * rod_vec[0];
            rmat[[2, 2]] = c + (1.0_f64 - c) * (rod_vec[2] * rod_vec[2]);
        };

        azip!((rmat in ori.axis_iter_mut(Axis(2)), rod_vec in self.ori.axis_iter(Axis(1))) {
            f(rmat, rod_vec);
        });

        RMat::new_init(ori)
    }//End of to_rmat

    ///Converts the Rodrigues vector representation over to axis-angle representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    fn to_ang_axis(&self) -> AngAxis{

        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((4, nelems).f());

        let f = |mut ang_axis: ArrayViewMut1::<f64>, ref rod_vec: ArrayView1::<f64>| {
            ang_axis[0] = rod_vec[0];
            ang_axis[1] = rod_vec[1];
            ang_axis[2] = rod_vec[2];
            ang_axis[3] = 2.0_f64 * rod_vec[3].atan();
        };

        azip!((ang_axis in ori.axis_iter_mut(Axis(1)), rod_vec in self.ori.axis_iter(Axis(1))) {
            f(ang_axis, rod_vec);
        });

        AngAxis::new_init(ori)
    }//End of to_ang_axis

    ///Converts the Rodrigues vector representation over to a compact axial vector representation which has the following properties
    ///shape (3, nelems), memory order = fortran/column major.
    fn to_ang_axis_comp(&self) -> AngAxisComp{
        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((3, nelems).f());

        let f = |mut ang_axis: ArrayViewMut1::<f64>, ref rod_vec: ArrayView1::<f64>| {
            let phi = 2.0_f64 * rod_vec[3].atan();
            ang_axis[0] = rod_vec[0] * phi;
            ang_axis[1] = rod_vec[1] * phi;
            ang_axis[2] = rod_vec[2] * phi;
        };

        azip!((ang_axis in ori.axis_iter_mut(Axis(1)), rod_vec in self.ori.axis_iter(Axis(1))) {
            f(ang_axis, rod_vec);
        });

        AngAxisComp::new_init(ori)
    }//End of to_ang_axis_comp

    ///Returns a clone of the Rodrigues vector.
    fn to_rod_vec(&self) -> RodVec{
        self.clone()
    }//End of to_rod_vec

    ///Converts the Rodrigues vector representation over to a compact Rodrigues which has the following properties
    ///shape (3, nelems), memory order = fortran/column major.
    fn to_rod_vec_comp(&self) -> RodVecComp{

        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((3, nelems).f());

        let f = |mut rod_vec_comp: ArrayViewMut1::<f64>, ref rod_vec: ArrayView1::<f64>| {
            rod_vec_comp[0] = rod_vec[0] * rod_vec[3];
            rod_vec_comp[1] = rod_vec[1] * rod_vec[3];
            rod_vec_comp[2] = rod_vec[2] * rod_vec[3];
        };

        azip!((rod_vec_comp in ori.axis_iter_mut(Axis(1)), rod_vec in self.ori.axis_iter(Axis(1))) {
            f(rod_vec_comp, rod_vec);
        });

        RodVecComp::new_init(ori)
    }//End of to_rod_vec_comp

    ///Converts the Rodrigues vector representation over to a unit quaternion which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    fn to_quat(&self) -> Quat{

        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((4, nelems).f());

        let f = |mut quat: ArrayViewMut1::<f64>, ref rod_vec: ArrayView1::<f64>| {
            let phi = rod_vec[3].atan();
            let s = f64::sin(phi); 

            quat[0] = f64::cos(phi);
            quat[1] = s * rod_vec[0];
            quat[2] = s * rod_vec[1];
            quat[3] = s * rod_vec[2];
        };

        azip!((quat in ori.axis_iter_mut(Axis(1)), rod_vec in self.ori.axis_iter(Axis(1))) {
            f(quat, rod_vec);
        });

        Quat::new_init(ori)
    }//End of to_quat

    ///Converts the Rodrigues vector representation over to a homochoric representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    fn to_homochoric(&self) -> Homochoric{
        let ang_axis = self.to_ang_axis();
        ang_axis.to_homochoric()
    }//End of to_homochoric

    ///Converts the Rodrigues vector representation over to Bunge angles which has the following properties
    ///shape (3, nelems), memory order = fortran/column major.
    ///This operation is done inplace and does not create a new structure
    fn to_bunge_inplace(&self, bunge: &mut Bunge){
        let rmat = self.to_rmat();
        rmat.to_bunge_inplace(bunge);
    }

    ///Converts the Rodrigues vector representation over to a rotation matrix which has the following properties
    ///shape (3, 3, nelems), memory order = fortran/column major.
    ///This operation is done inplace and does not create a new structure
    fn to_rmat_inplace(&self, rmat: &mut RMat){

        let mut ori = rmat.ori_view_mut();

        let new_nelem = ori.len_of(Axis(1));
        let nelem = self.ori.len_of(Axis(1));

        assert!(new_nelem == nelem, 
        "The number of elements in the original ori field do no match up with the new field.
        The old field had {} elements, and the new field has {} elements",
        nelem, new_nelem);

        let f = |mut rmat: ArrayViewMut2::<f64>, ref rod_vec: ArrayView1::<f64>| {
            let phi = rod_vec[3].atan() * 2.0_f64;
            let c = phi.cos();
            let s = phi.sin();

            rmat[[0, 0]] = c + (1.0_f64 - c) * (rod_vec[0] * rod_vec[0]);
            rmat[[1, 0]] = (1.0_f64 - c) * (rod_vec[0] * rod_vec[1]) + s * rod_vec[2];
            rmat[[2, 0]] = (1.0_f64 - c) * (rod_vec[0] * rod_vec[2]) - s * rod_vec[1];

            rmat[[0, 1]] = (1.0_f64 - c) * (rod_vec[0] * rod_vec[1]) - s * rod_vec[2];
            rmat[[1, 1]] = c + (1.0_f64 - c) * (rod_vec[1] * rod_vec[1]);
            rmat[[2, 1]] = (1.0_f64 - c) * (rod_vec[1] * rod_vec[2]) + s * rod_vec[0];

            rmat[[0, 2]] = (1.0_f64 - c) * (rod_vec[0] * rod_vec[2]) + s * rod_vec[1];
            rmat[[1, 2]] = (1.0_f64 - c) * (rod_vec[1] * rod_vec[2]) - s * rod_vec[0];
            rmat[[2, 2]] = c + (1.0_f64 - c) * (rod_vec[2] * rod_vec[2]);
        };

        azip!((rmat in ori.axis_iter_mut(Axis(2)), rod_vec in self.ori.axis_iter(Axis(1))) {
            f(rmat, rod_vec);
        });
    }

    ///Converts the Rodrigues vector representation over to an axial vector representation which has the following properties
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

        let f = |mut ang_axis: ArrayViewMut1::<f64>, ref rod_vec: ArrayView1::<f64>| {
            ang_axis[0] = rod_vec[0];
            ang_axis[1] = rod_vec[1];
            ang_axis[2] = rod_vec[2];
            ang_axis[3] = 2.0_f64 * rod_vec[3].atan();
        };

        azip!((ang_axis in ori.axis_iter_mut(Axis(1)), rod_vec in self.ori.axis_iter(Axis(1))) {
            f(ang_axis, rod_vec);
        });
    }

    ///Converts the Rodrigues vector representation over to a compact axial vector representation which has the following properties
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

        let f = |mut ang_axis: ArrayViewMut1::<f64>, ref rod_vec: ArrayView1::<f64>| {
            let phi = 2.0_f64 * rod_vec[3].atan();
            ang_axis[0] = rod_vec[0] * phi;
            ang_axis[1] = rod_vec[1] * phi;
            ang_axis[2] = rod_vec[2] * phi;
        };

        azip!((ang_axis in ori.axis_iter_mut(Axis(1)), rod_vec in self.ori.axis_iter(Axis(1))) {
            f(ang_axis, rod_vec);
        });
    }

    ///Returns a clone of the Rodrigues vector.
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

    ///Converts the Rodrigues vector representation over to a compact Rodrigues which has the following properties
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

        let f = |mut rod_vec_comp: ArrayViewMut1::<f64>, ref rod_vec: ArrayView1::<f64>| {
            rod_vec_comp[0] = rod_vec[0] * rod_vec[3];
            rod_vec_comp[1] = rod_vec[1] * rod_vec[3];
            rod_vec_comp[2] = rod_vec[2] * rod_vec[3];
        };

        azip!((rod_vec_comp in ori.axis_iter_mut(Axis(1)), rod_vec in self.ori.axis_iter(Axis(1))) {
            f(rod_vec_comp, rod_vec);
        });

    }

    ///Converts the Rodrigues vector representation over to a unit quaternion which has the following properties
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

        let f = |mut quat: ArrayViewMut1::<f64>, ref rod_vec: ArrayView1::<f64>| {
            let phi = rod_vec[3].atan();
            let s = f64::sin(phi); 

            quat[0] = f64::cos(phi);
            quat[1] = s * rod_vec[0];
            quat[2] = s * rod_vec[1];
            quat[3] = s * rod_vec[2];
        };

        azip!((quat in ori.axis_iter_mut(Axis(1)), rod_vec in self.ori.axis_iter(Axis(1))) {
            f(quat, rod_vec);
        });
    }

    ///Converts the Rodrigues vector representation over to a homochoric representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
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

            azip!((rvec in rvec.axis_iter_mut(Axis(1)), ref vec in vec.axis_iter(Axis(1)), 
            ref rod_vec in self.ori.axis_iter(Axis(1))) {
                rod_vec_rot_vec(&rod_vec, &vec, rvec);     
            });
        } else if rnelems == 1{
            //We just have one Rodrigues vector so perform pretty much the above to get all of our values
            let rod_vec = self.ori.index_axis(Axis(1), 0);

            azip!((rvec in rvec.axis_iter_mut(Axis(1)), ref vec in vec.axis_iter(Axis(1))) {  
                rod_vec_rot_vec(&rod_vec, &vec, rvec); 
            });
        }else{
            //We just have one Rodrigues vector so perform pretty much the above to get all of our values
            let vec = vec.index_axis(Axis(1), 0);

            azip!((rvec in rvec.axis_iter_mut(Axis(1)), ref rod_vec in self.ori.axis_iter(Axis(1))) {  
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
            azip!((rvec in rvec.axis_iter_mut(Axis(1)), ref vec in vec.axis_iter(Axis(1)), 
            ref rod_vec in self.ori.axis_iter(Axis(1))) {
                rod_vec_rot_vec(&rod_vec, &vec, rvec);    
            });
        } else if rnelems == 1{
            //We just have one Rodrigues vector so perform pretty much the above to get all of our values
            let rod_vec = self.ori.index_axis(Axis(1), 0);

            azip!((rvec in rvec.axis_iter_mut(Axis(1)), ref vec in vec.axis_iter(Axis(1))) {  
                rod_vec_rot_vec(&rod_vec, &vec, rvec); 
            });
        }else{
            //We just have one Rodrigues vector so perform pretty much the above to get all of our values
            let vec = vec.index_axis(Axis(1), 0);

            azip!((rvec in rvec.axis_iter_mut(Axis(1)), ref rod_vec in self.ori.axis_iter(Axis(1))) {  
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
            azip!((mut vec in vec.axis_iter_mut(Axis(1)), ref rod_vec in self.ori.axis_iter(Axis(1))) {
                let mut rvec = Array1::<f64>::zeros((3).f());
                rod_vec_rot_vec(&rod_vec, &vec.view(), rvec.view_mut());
                vec.assign({&rvec});    
            });
        } else{
            //We just have one Rodrigues vector so perform pretty much the above to get all of our values
            let rod_vec = self.ori.index_axis(Axis(1), 0);

            azip!((mut vec in vec.axis_iter_mut(Axis(1))) {
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
    cross_prod[1] = rod_vec[2] * vec[0] - rod_vec[0] * vec[2];
    cross_prod[2] = -rod_vec[1] * vec[0] + rod_vec[0] * vec[1];

    rvec[0] = vec[0] * cos_theta + cross_prod[0] * sin_theta + rod_vec[0] * dp_mcos;
    rvec[1] = vec[1] * cos_theta + cross_prod[1] * sin_theta + rod_vec[1] * dp_mcos;
    rvec[2] = vec[2] * cos_theta + cross_prod[2] * sin_theta + rod_vec[2] * dp_mcos;
}