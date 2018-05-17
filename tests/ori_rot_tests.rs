extern crate mori;
#[macro_use]
extern crate ndarray;
extern crate csv;

use mori::orientations::*;
use csv::ReaderBuilder;
use ndarray::prelude::*;


#[test]
fn rmat_rot_vec(){
    //First test does the rmat and vec have same number of elems
    let nori = 1551;
    let vec = read_rod_comp_file();

    let rmat = RMat::new(nori);

    let rvec = rmat.rot_vector(vec.view());

    let comp = vec.all_close(&rvec, 1e-14);

    assert!(comp);
    //Second test does the rmat has 1 elem and vec has nelems
    let rmat = RMat::new(1);

    let rvec = rmat.rot_vector(vec.view());

    let comp = vec.all_close(&rvec, 1e-14);

    assert!(comp);
    //Last case does the rmat has nelems and vec has 1 elem
    let rmat = RMat::new(nori);

    let mut vec_comp = Array2::<f64>::zeros((3, nori).f());

    let mut vec = Array2::<f64>::zeros((3,1).f());

    vec[[0,0]] = 1.0_f64;

    azip!(mut vec (vec_comp.axis_iter_mut(Axis(1))) in {vec[0] = 1.0_f64});

    let rvec = rmat.rot_vector(vec.view());

    let comp = vec_comp.all_close(&rvec, 1e-14);

    assert!(comp);
}


#[test]
fn rmat_rot_vec_mut(){
    //First test does the rmat and vec have same number of elems
    let nori = 1551;
    let vec = read_rod_comp_file();
    let mut rvec = Array2::<f64>::zeros((3, nori).f());

    let rmat = RMat::new(nori);

    rmat.rot_vector_mut(vec.view(), rvec.view_mut());

    let comp = vec.all_close(&rvec, 1e-14);

    assert!(comp);
    //Second test does the rmat has 1 elem and vec has nelems
    let rmat = RMat::new(1);
    let mut rvec = Array2::<f64>::zeros((3, nori).f());

    rmat.rot_vector_mut(vec.view(), rvec.view_mut());


    let comp = vec.all_close(&rvec, 1e-14);

    assert!(comp);
    //Last case does the rmat has nelems and vec has 1 elem
    let rmat = RMat::new(nori);

    let mut vec_comp = Array2::<f64>::zeros((3, nori).f());
    let mut rvec = Array2::<f64>::zeros((3, nori).f());

    let mut vec = Array2::<f64>::zeros((3,1).f());

    vec[[0,0]] = 1.0_f64;

    azip!(mut vec (vec_comp.axis_iter_mut(Axis(1))) in {vec[0] = 1.0_f64});

    rmat.rot_vector_mut(vec.view(), rvec.view_mut());

    let comp = vec_comp.all_close(&rvec, 1e-14);

    assert!(comp);
}

#[test]
fn rmat_rot_vec_inplace(){
    //First test does the rmat and vec have same number of elems
    let nori = 1551;
    let mut vec = read_rod_comp_file();
    let vec_comp = vec.clone();

    let rmat = RMat::new(nori);

    rmat.rot_vector_inplace(vec.view_mut());

    let comp = vec_comp.all_close(&vec, 1e-14);

    assert!(comp);
    //Second test does the rmat has 1 elem and vec has nelems
    let rmat = RMat::new(1);
    let mut vec = vec_comp.clone();

    rmat.rot_vector_inplace(vec.view_mut());


    let comp = vec_comp.all_close(&vec, 1e-14);

    assert!(comp);
}

#[test]
fn quat_rot_vec(){
    //First test does the rmat and vec have same number of elems
    let nori = 1551;
    let ori = read_rod_file();

    let rod = RodVec::new_init(ori.clone());

    let rmat = rod.to_rmat();
    let quat = rod.to_quat();

    let vec = Array2::<f64>::ones((3, nori).f());

    let rvec_rmat = rmat.rot_vector(vec.view());
    let rvec_quat = quat.rot_vector(vec.view());

    let comp = rvec_rmat.all_close(&rvec_quat, 1e-14);

    assert!(comp);
    //Second test does the rmat has 1 elem and vec has nelems
    let rmat = RMat::new(1);
    let quat = Quat::new(1);

    let rvec_rmat = rmat.rot_vector(vec.view());
    let rvec_quat = quat.rot_vector(vec.view());

    let comp = rvec_rmat.all_close(&rvec_quat, 1e-14);

    assert!(comp);
}


#[test]
fn quat_rot_vec_mut(){
    //First test does the rmat and vec have same number of elems
    let nori = 1551;
    let ori = read_rod_file();

    let rod = RodVec::new_init(ori.clone());

    let rmat = rod.to_rmat();
    let quat = rod.to_quat();

    let vec = Array2::<f64>::ones((3, nori).f());
    let mut rvec_rmat = vec.clone();
    let mut rvec_quat = vec.clone();

    rmat.rot_vector_mut(vec.view(), rvec_rmat.view_mut());
    quat.rot_vector_mut(vec.view(), rvec_quat.view_mut());

    let comp = rvec_rmat.all_close(&rvec_quat, 1e-14);

    assert!(comp);
    //Second test does the rmat has 1 elem and vec has nelems
    let rmat = RMat::new(1);
    let quat = Quat::new(1);

    let mut rvec_rmat = vec.clone();
    let mut rvec_quat = vec.clone();

    rmat.rot_vector_mut(vec.view(), rvec_rmat.view_mut());
    quat.rot_vector_mut(vec.view(), rvec_quat.view_mut());

    let comp = rvec_rmat.all_close(&rvec_quat, 1e-14);

    assert!(comp);
    //Last case does the rmat has nelems and vec has 1 elem
    let rmat = RMat::new(nori);
    let quat = Quat::new(nori);

    let mut rvec_rmat = Array2::<f64>::zeros((3, nori).f());
    let mut rvec_quat = Array2::<f64>::zeros((3, nori).f());

    let mut vec = Array2::<f64>::zeros((3,1).f());

    vec[[0,0]] = 1.0_f64;

    rmat.rot_vector_mut(vec.view(), rvec_rmat.view_mut());
    quat.rot_vector_mut(vec.view(), rvec_quat.view_mut());

    let comp = rvec_rmat.all_close(&rvec_quat, 1e-14);

    assert!(comp);
}

#[test]
fn quat_rot_vec_inplace(){
    //First test does the rmat and vec have same number of elems
    let nori = 1551;
    let ori = read_rod_file();

    let rod = RodVec::new_init(ori.clone());

    let rmat = rod.to_rmat();
    let quat = rod.to_quat();

    let vec = Array2::<f64>::ones((3, nori).f());
    let mut rvec_rmat = vec.clone();
    let mut rvec_quat = vec.clone();

    rmat.rot_vector_inplace(rvec_rmat.view_mut());
    quat.rot_vector_inplace(rvec_quat.view_mut());

    let comp = rvec_rmat.all_close(&rvec_quat, 1e-14);

    assert!(comp);
    //Second test does the rmat has 1 elem and vec has nelems
    let rmat = RMat::new(1);
    let quat = Quat::new(1);

    let mut rvec_rmat = vec.clone();
    let mut rvec_quat = vec.clone();

    rmat.rot_vector_inplace(rvec_rmat.view_mut());
    quat.rot_vector_inplace(rvec_quat.view_mut());

    let comp = rvec_rmat.all_close(&rvec_quat, 1e-14);

    assert!(comp);
}

#[test]
fn angaxis_rot_vec(){
    //First test does the rmat and vec have same number of elems
    let nori = 1551;
    let ori = read_rod_file();

    let rod = RodVec::new_init(ori.clone());

    let rmat = rod.to_rmat();
    let ang_axis = rod.to_ang_axis();

    let vec = Array2::<f64>::ones((3, nori).f());

    let rvec_rmat = rmat.rot_vector(vec.view());
    let rvec_angaxis = ang_axis.rot_vector(vec.view());

    let comp = rvec_rmat.all_close(&rvec_angaxis, 1e-14);

    assert!(comp);
    //Second test does the rmat has 1 elem and vec has nelems
    let rmat = RMat::new(1);
    let ang_axis = AngAxis::new(1);

    let rvec_rmat = rmat.rot_vector(vec.view());
    let rvec_angaxis = ang_axis.rot_vector(vec.view());

    let comp = rvec_rmat.all_close(&rvec_angaxis, 1e-14);

    assert!(comp);

}


#[test]
fn angaxis_rot_vec_mut(){
    //First test does the rmat and vec have same number of elems
    let nori = 1551;
    let ori = read_rod_file();

    let rod = RodVec::new_init(ori.clone());

    let rmat = rod.to_rmat();
    let ang_axis = rod.to_ang_axis();

    let vec = Array2::<f64>::ones((3, nori).f());
    let mut rvec_rmat = vec.clone();
    let mut rvec_angaxis = vec.clone();

    rmat.rot_vector_mut(vec.view(), rvec_rmat.view_mut());
    ang_axis.rot_vector_mut(vec.view(), rvec_angaxis.view_mut());

    let comp = rvec_rmat.all_close(&rvec_angaxis, 1e-14);

    assert!(comp);
    //Second test does the rmat has 1 elem and vec has nelems
    let rmat = RMat::new(1);
    let ang_axis = AngAxis::new(1);

    let mut rvec_rmat = vec.clone();
    let mut rvec_angaxis = vec.clone();

    rmat.rot_vector_mut(vec.view(), rvec_rmat.view_mut());
    ang_axis.rot_vector_mut(vec.view(), rvec_angaxis.view_mut());

    let comp = rvec_rmat.all_close(&rvec_angaxis, 1e-14);

    assert!(comp);

    //Last case does the rmat has nelems and vec has 1 elem
    let rmat = RMat::new(nori);
    let ang_axis = AngAxis::new(nori);

    let mut rvec_rmat = Array2::<f64>::zeros((3, nori).f());
    let mut rvec_angaxis = Array2::<f64>::zeros((3, nori).f());

    let mut vec = Array2::<f64>::zeros((3,1).f());

    vec[[0,0]] = 1.0_f64;

    rmat.rot_vector_mut(vec.view(), rvec_rmat.view_mut());
    ang_axis.rot_vector_mut(vec.view(), rvec_angaxis.view_mut());

    let comp = rvec_rmat.all_close(&rvec_angaxis, 1e-14);

    assert!(comp);
}

#[test]
fn angaxis_rot_vec_inplace(){
    //First test does the rmat and vec have same number of elems
    let nori = 1551;
    let ori = read_rod_file();

    let rod = RodVec::new_init(ori.clone());

    let rmat = rod.to_rmat();
    let ang_axis = rod.to_ang_axis();

    let vec = Array2::<f64>::ones((3, nori).f());
    let mut rvec_rmat = vec.clone();
    let mut rvec_angaxis = vec.clone();

    rmat.rot_vector_inplace(rvec_rmat.view_mut());
    ang_axis.rot_vector_inplace(rvec_angaxis.view_mut());

    let comp = rvec_rmat.all_close(&rvec_angaxis, 1e-14);

    assert!(comp);
    //Second test does the rmat has 1 elem and vec has nelems
    let rmat = RMat::new(1);
    let ang_axis = AngAxis::new(1);

    let mut rvec_rmat = vec.clone();
    let mut rvec_angaxis = vec.clone();

    rmat.rot_vector_inplace(rvec_rmat.view_mut());
    ang_axis.rot_vector_inplace(rvec_angaxis.view_mut());

    let comp = rvec_rmat.all_close(&rvec_angaxis, 1e-14);

    assert!(comp);
}


#[test]
fn rodvec_rot_vec(){
    //First test does the rmat and vec have same number of elems
    let nori = 1551;
    let ori = read_rod_file();

    let rod = RodVec::new_init(ori.clone());

    let rmat = rod.to_rmat();

    let vec = Array2::<f64>::ones((3, nori).f());

    let rvec_rmat = rmat.rot_vector(vec.view());
    let rvec_rodvec = rod.rot_vector(vec.view());

    let comp = rvec_rmat.all_close(&rvec_rodvec, 1e-14);

    assert!(comp);
    //Second test does the rmat has 1 elem and vec has nelems
    let rmat = RMat::new(1);
    let rod  = RodVec::new(1);

    let vec = Array2::<f64>::ones((3, nori).f());

    let rvec_rmat = rmat.rot_vector(vec.view());
    let rvec_rodvec = rod.rot_vector(vec.view());

    let comp = rvec_rmat.all_close(&rvec_rodvec, 1e-14);

    assert!(comp);
}


#[test]
fn rodvec_rot_vec_mut(){
    //First test does the rmat and vec have same number of elems
    let nori = 1551;
    let ori = read_rod_file();

    let rod = RodVec::new_init(ori.clone());

    let rmat = rod.to_rmat();
    
    let vec = Array2::<f64>::ones((3, nori).f());
    let mut rvec_rmat = vec.clone();
    let mut rvec_rodvec = vec.clone();

    rmat.rot_vector_mut(vec.view(), rvec_rmat.view_mut());
    rod.rot_vector_mut(vec.view(), rvec_rodvec.view_mut());

    let comp = rvec_rmat.all_close(&rvec_rodvec, 1e-14);

    assert!(comp);
    //Second test does the rmat has 1 elem and vec has nelems
    let rmat = RMat::new(1);
    let rod  = RodVec::new(1);
    let mut rvec_rmat = vec.clone();
    let mut rvec_rodvec = vec.clone();

    rmat.rot_vector_mut(vec.view(), rvec_rmat.view_mut());
    rod.rot_vector_mut(vec.view(), rvec_rodvec.view_mut());

    let comp = rvec_rmat.all_close(&rvec_rodvec, 1e-14);

    assert!(comp);

    //Last case does the rmat has nelems and vec has 1 elem
    let rmat = RMat::new(nori);
    let rod_vec = RodVec::new(nori);

    let mut rvec_rmat = Array2::<f64>::zeros((3, nori).f());
    let mut rvec_rodvec = Array2::<f64>::zeros((3, nori).f());

    let mut vec = Array2::<f64>::zeros((3,1).f());

    vec[[0,0]] = 1.0_f64;

    rmat.rot_vector_mut(vec.view(), rvec_rmat.view_mut());
    rod_vec.rot_vector_mut(vec.view(), rvec_rodvec.view_mut());

    let comp = rvec_rmat.all_close(&rvec_rodvec, 1e-14);

    assert!(comp);
}

#[test]
fn par_rodvec_rot_vec_mut(){
    //First test does the rmat and vec have same number of elems
    let nori = 1551;
    let ori = read_rod_file();

    let rod = RodVec::new_init(ori.clone());

    let rmat = rod.to_rmat();
    
    let vec = Array2::<f64>::ones((3, nori).f());
    let mut rvec_rmat = vec.clone();
    let mut rvec_rodvec = vec.clone();

    rmat.rot_vector_mut(vec.view(), rvec_rmat.view_mut());
    rod.par_rot_vector_mut(vec.view(), rvec_rodvec.view_mut());

    let comp = rvec_rmat.all_close(&rvec_rodvec, 1e-14);

    assert!(comp);
    //Second test does the rmat has 1 elem and vec has nelems
    let rmat = RMat::new(1);
    let rod  = RodVec::new(1);
    let mut rvec_rmat = vec.clone();
    let mut rvec_rodvec = vec.clone();

    rmat.rot_vector_mut(vec.view(), rvec_rmat.view_mut());
    rod.par_rot_vector_mut(vec.view(), rvec_rodvec.view_mut());

    let comp = rvec_rmat.all_close(&rvec_rodvec, 1e-14);

    assert!(comp);

    //Last case does the rmat has nelems and vec has 1 elem
    let rmat = RMat::new(nori);
    let rod_vec = RodVec::new(nori);

    let mut rvec_rmat = Array2::<f64>::zeros((3, nori).f());
    let mut rvec_rodvec = Array2::<f64>::zeros((3, nori).f());

    let mut vec = Array2::<f64>::zeros((3,1).f());

    vec[[0,0]] = 1.0_f64;

    rmat.rot_vector_mut(vec.view(), rvec_rmat.view_mut());
    rod_vec.par_rot_vector_mut(vec.view(), rvec_rodvec.view_mut());

    let comp = rvec_rmat.all_close(&rvec_rodvec, 1e-14);

    assert!(comp);
}

#[test]
fn rodvec_rot_vec_inplace(){
    //First test does the rmat and vec have same number of elems
    let nori = 1551;
    let ori = read_rod_file();

    let rod = RodVec::new_init(ori.clone());

    let rmat = rod.to_rmat();

    let vec = Array2::<f64>::ones((3, nori).f());
    let mut rvec_rmat = vec.clone();
    let mut rvec_rodvec = vec.clone();

    rmat.rot_vector_inplace(rvec_rmat.view_mut());
    rod.rot_vector_inplace(rvec_rodvec.view_mut());

    let comp = rvec_rmat.all_close(&rvec_rodvec, 1e-14);

    assert!(comp);
    //Second test does the rmat has 1 elem and vec has nelems
    let rmat = RMat::new(1);
    let rod = RodVec::new(1);

    rmat.rot_vector_inplace(rvec_rmat.view_mut());
    rod.rot_vector_inplace(rvec_rodvec.view_mut());

    let comp = rvec_rmat.all_close(&rvec_rodvec, 1e-14);

    assert!(comp);
}

fn read_rod_file() -> Array2<f64>{
    let nori = 1551;
    let mut rod_ori = Array2::<f64>::zeros((4, nori).f());

    let mut i = 0;

    let mut rdr = ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b' ')
        .from_path("tests/rods.txt").unwrap();

    for result in rdr.records(){
        let record = result.unwrap();
        
        rod_ori[[0, i]] = record[0].parse().unwrap();
        rod_ori[[1, i]] = record[1].parse().unwrap();
        rod_ori[[2, i]] = record[2].parse().unwrap();
        rod_ori[[3, i]] = record[3].parse().unwrap();

        i += 1;
    }

    rod_ori
}


fn read_rod_comp_file() -> Array2<f64>{
    let nori = 1551;
    let mut rod_comp_ori = Array2::<f64>::zeros((3, nori).f());

    let mut rdr = ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b' ')
        .from_path("tests/comp_rods.txt").unwrap();

    let mut i = 0;

    for result in rdr.records(){
        let record = result.unwrap();
        
        rod_comp_ori[[0, i]] = record[0].parse().unwrap();
        rod_comp_ori[[1, i]] = record[1].parse().unwrap();
        rod_comp_ori[[2, i]] = record[2].parse().unwrap();

        i += 1;
    }

    rod_comp_ori
}