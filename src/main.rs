use std::env;
use std::fs::File;
use std::io::prelude::*;

const NSPEEDS: usize = 9;
const FINALSTATEFILE: &'static str = "final_state.dat";
const AVVELSFILE: &'static str = "av_vels.dat";

/* struct to hold the parameter values */
struct Param {
    nx: usize,            /* no. of cells in x-direction */
    ny: usize,           /* no. of cells in y-direction */
    max_iters: usize,      /* no. of iterations */
    reynolds_dim: i32,  /* dimension for Reynolds number */
    density: f64,       /* density per link */
    accel: f64,         /* density redistribution */
    omega: f64,         /* relaxation parameter */
}

/* struct to hold the 'speed' values */

#[derive(Clone, Copy)]
struct Speed {
    speeds: [f64; NSPEEDS],
}

fn main() {
    println!("Hello rust!");

    let args: Vec<String> = env::args().collect();

    let paramfile = &args[1];
    let obstaclefile = &args[2];

    println!("Searching for {}", paramfile);
    println!("In file {}", obstaclefile);

    let (params, cells, tmp_cells, obstacles, av_vels) = initialise(paramfile, obstaclefile);
    println!("nx: {}", params.nx);

}


// int initialise(const char* paramfile, const char* obstaclefile,
//                t_param* params, t_speed** cells_ptr, t_speed** tmp_cells_ptr,
//                int** obstacles_ptr, double** av_vels_ptr)
// {
fn initialise(paramfile: &str, obstaclefile: &str) -> (Param, Vec<Speed>, Vec<Speed>, Vec<u8>, Vec<f64>) {
    // int    xx, yy;         /* generic array indices */
    // int    blocked;        /* indicates whether a cell is blocked by an obstacle */
    // int    retval;         /* to hold return value for checking */


    let mut f = File::open(paramfile).expect("could not open input parameter file");

    let mut contents = String::new();
    f.read_to_string(&mut contents)
        .expect("something went wrong reading the file");


    println!("{}", contents);
    let mut lines = contents.lines();

    let nx: usize = lines.next().unwrap().parse().unwrap();
    let ny: usize = lines.next().unwrap().parse().unwrap();
    let max_iters: usize = lines.next().unwrap().parse().unwrap();
    let reynolds_dim: i32 = lines.next().unwrap().parse().unwrap();
    let density: f64 = lines.next().unwrap().parse().unwrap();
    let accel: f64 = lines.next().unwrap().parse().unwrap();
    let omega: f64 = lines.next().unwrap().parse().unwrap();

    let params = Param{nx, ny, max_iters, reynolds_dim, density, accel, omega};


    // /* main grid */
    //let mut cells: Vec<Speed> = Vec::with_capacity(params.ny * params.nx);
    //let mut tmp_cells: Vec<Speed> = Vec::with_capacity(params.ny * params.nx);

    let mut cells: Vec<Speed> = vec![Speed{speeds: [0.0; NSPEEDS]}; params.ny * params.nx];
    let mut tmp_cells: Vec<Speed> = vec![Speed{speeds: [0.0; NSPEEDS]}; params.ny * params.nx];


    /* initialise densities */
    let w0 = params.density * 4.0 / 9.0;
    let w1 = params.density      / 9.0;
    let w2 = params.density      / 36.0;
    println! ("Here") ;
    for ii in 0..params.ny {
        for jj in 0..params.nx {
            /* centre */
            cells[ii * params.nx + jj].speeds[0] = w0;
            /* axis directions */
            cells[ii * params.nx + jj].speeds[1] = w1;
            cells[ii * params.nx + jj].speeds[2] = w1;
            cells[ii * params.nx + jj].speeds[3] = w1;
            cells[ii * params.nx + jj].speeds[4] = w1;
            /* diagonals */
            cells[ii * params.nx + jj].speeds[5] = w2;
            cells[ii * params.nx + jj].speeds[6] = w2;
            cells[ii * params.nx + jj].speeds[7] = w2;
            cells[ii * params.nx + jj].speeds[8] = w2;
        }
    }

    /* the map of obstacles */
    let mut obstacles = vec![0 as u8; (params.ny * params.nx)];

    /* open the obstacle data file */
    let mut f = File::open(obstaclefile).expect("could not open obstacle parameter file");

    let mut contents = String::new();
    f.read_to_string(&mut contents)
        .expect("something went wrong reading the obstacle file");

    let mut lines = contents.lines();
    for line in lines {
        let res: Vec<u16> = line.split(" ").map(|s| s.parse().unwrap()).collect();
        obstacles[res[1] as usize * params.nx + res[0] as usize] = res[2] as u8;
        println!("{} {} {}", res[0], res[1], res[2]);
    }

    /*
     ** allocate space to hold a record of the avarage velocities computed
     ** at each timestep
     */
    let mut av_vels = vec![0.0; params.max_iters];

    return (params, cells, tmp_cells, obstacles, av_vels);
}