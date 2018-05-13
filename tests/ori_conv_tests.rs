extern crate mori;
extern crate ndarray;
extern crate csv;

use mori::orientations::*;
use csv::ReaderBuilder;
use ndarray::prelude::*;

//All of the below tests use a set of Rodrigues Vector data pulled from a Fundamental
//Region of cubic crystals located at the origin. It therefore provides a comprehensive
//set of conversions to try against.
#[test]
fn rod2rodcomp(){

    let rod_ori = read_rod_file();

    let rod_comp_ori = read_rod_comp_file();

    let rod = RodVec::new_init(rod_ori);

    let rod_comp = rod.to_rod_vec_comp();

    let rod_comp_ori2 = rod_comp.ori_view();

    let comp = rod_comp_ori.all_close(&rod_comp_ori2, 1e-14);

    assert!(comp);
}

#[test]
fn rod_edit(){

    let rod_ori = read_rod_file();

    let mut rod = RodVec::new_init(rod_ori.clone());
    //If we want to actually change the value of rod.ori we need to place the changes in a limited scope
    //such that the mutable reference is dropped once we are done changing stuff. If we don't the compiler
    //will complain about the mutable reference not being dropped until later.
    {
        let mut rod_ori2 = rod.ori_view_mut();

        rod_ori2[[0,0]] = 1.0_f64;

        assert_eq!(rod_ori2[[0,0]], 1.0_f64);

    }

    let rod2 = rod.to_rod_vec();

    let rod_ori2 = rod2.ori_view();
    assert_eq!(rod_ori2[[0,0]], 1.0_f64);

}


#[test]
fn rod2angaxis2rod(){

    let rod_ori = read_rod_file();

    let rod = RodVec::new_init(rod_ori.clone());

    let ang_axis = rod.to_ang_axis();
    let rod2 = ang_axis.to_rod_vec();
    let rod_ori2 = rod2.ori_view();

    let comp = rod_ori.all_close(&rod_ori2, 1e-14);

    assert!(comp);
}

#[test]
fn rod2angaxiscomp2rod(){

    let rod_ori = read_rod_file();

    let rod = RodVec::new_init(rod_ori.clone());

    let ang_axis_comp = rod.to_ang_axis_comp();
    let rod2 = ang_axis_comp.to_rod_vec();
    let rod_ori2 = rod2.ori_view();

    let comp = rod_ori.all_close(&rod_ori2, 1e-14);

    assert!(comp);
}


#[test]
fn rod2rmat2rod(){

    let rod_ori = read_rod_file();

    let rod = RodVec::new_init(rod_ori.clone());

    let rmat = rod.to_rmat();
    let rod2 = rmat.to_rod_vec();
    let rod_ori2 = rod2.ori_view();

    let comp = rod_ori.all_close(&rod_ori2, 1e-14);

    assert!(comp);
}

#[test]
fn rod2bunge2rod(){

    let rod_ori = read_rod_file();

    let rod = RodVec::new_init(rod_ori.clone());

    let bunge = rod.to_bunge();
    let rod2 = bunge.to_rod_vec();
    let rod_ori2 = rod2.ori_view();

    let comp = rod_ori.all_close(&rod_ori2, 1e-14);

    assert!(comp);
}


#[test]
fn rod2quat2rod(){

    let rod_ori = read_rod_file();

    let rod = RodVec::new_init(rod_ori.clone());

    let quat = rod.to_quat();
    let rod2 = quat.to_rod_vec();
    let rod_ori2 = rod2.ori_view();

    let comp = rod_ori.all_close(&rod_ori2, 1e-14);

    assert!(comp);
}

#[test]
fn rmat2bunge2rmat(){

    let rod_ori = read_rod_file();

    let rod = RodVec::new_init(rod_ori.clone());

    let rmat = rod.to_rmat();
    let bunge = rmat.to_bunge();
    let rmat2 = bunge.to_rmat();

    let rmat_ori = rmat.ori_view();
    let rmat_ori2 = rmat2.ori_view();

    let comp = rmat_ori.all_close(&rmat_ori2, 1e-14);

    assert!(comp);
}

#[test]
fn quat2bunge2quat(){

    let rod_ori = read_rod_file();

    let rod = RodVec::new_init(rod_ori.clone());

    let quat = rod.to_quat();
    let bunge = quat.to_bunge();
    let quat2 = bunge.to_quat();

    let quat_ori = quat.ori_view();
    let quat_ori2 = quat2.ori_view();

    let comp = quat_ori.all_close(&quat_ori2, 1e-14);

    assert!(comp);
}

#[test]
fn quat2rmat2quat(){

    let rod_ori = read_rod_file();

    let rod = RodVec::new_init(rod_ori.clone());

    let quat = rod.to_quat();
    let rmat = quat.to_rmat();
    let quat2 = rmat.to_quat();

    let quat_ori = quat.ori_view();
    let quat_ori2 = quat2.ori_view();

    let comp = quat_ori.all_close(&quat_ori2, 1e-14);

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