use std::env;
use std::fs::File;
use std::io::prelude::*;

const NSPEEDS: usize = 9;
const FINALSTATEFILE: &'static str = "final_state.dat";
const AVVELSFILE: &'static str = "av_vels.dat";

/* struct to hold the parameter values */
struct Param {
    nx: i32,            /* no. of cells in x-direction */
    ny: i32,           /* no. of cells in y-direction */
    max_iters: i32,      /* no. of iterations */
    reynolds_dim: i32,  /* dimension for Reynolds number */
    density: f64,       /* density per link */
    accel: f64,         /* density redistribution */
    omega: f64,         /* relaxation parameter */
}

/* struct to hold the 'speed' values */
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

    let params = initialise(paramfile);
    println!("nx: {}", params.nx);
}


// int initialise(const char* paramfile, const char* obstaclefile,
//                t_param* params, t_speed** cells_ptr, t_speed** tmp_cells_ptr,
//                int** obstacles_ptr, double** av_vels_ptr)
// {
fn initialise(paramfile: &str) -> Param {
  // char   message[1024];  /* message buffer */
  // FILE*   fp;            /* file pointer */
  // int    xx, yy;         /* generic array indices */
  // int    blocked;        /* indicates whether a cell is blocked by an obstacle */
  // int    retval;         /* to hold return value for checking */


    let mut f = File::open(paramfile).expect("could not open input parameter file");

    let mut contents = String::new();
    f.read_to_string(&mut contents)
        .expect("something went wrong reading the file");


    println!("{}", contents);
    let mut lines = contents.lines();

    let nx: i32 = lines.next().unwrap().parse().unwrap();
    let ny: i32 = lines.next().unwrap().parse().unwrap();
    let max_iters: i32 = lines.next().unwrap().parse().unwrap();
    let reynolds_dim: i32 = lines.next().unwrap().parse().unwrap();
    let density: f64 = lines.next().unwrap().parse().unwrap();
    let accel: f64 = lines.next().unwrap().parse().unwrap();
    let omega: f64 = lines.next().unwrap().parse().unwrap();

    let params = Param{nx, ny, max_iters, reynolds_dim, density, accel, omega};


    return params;


  // /* read in the parameter values */
  // retval = fscanf(fp, "%d\n", &(params->nx));

  // if (retval != 1) die("could not read param file: nx", __LINE__, __FILE__);

  // retval = fscanf(fp, "%d\n", &(params->ny));

  // if (retval != 1) die("could not read param file: ny", __LINE__, __FILE__);

  // retval = fscanf(fp, "%d\n", &(params->maxIters));

  // if (retval != 1) die("could not read param file: maxIters", __LINE__, __FILE__);

  // retval = fscanf(fp, "%d\n", &(params->reynolds_dim));

  // if (retval != 1) die("could not read param file: reynolds_dim", __LINE__, __FILE__);

  // retval = fscanf(fp, "%lf\n", &(params->density));

  // if (retval != 1) die("could not read param file: density", __LINE__, __FILE__);

  // retval = fscanf(fp, "%lf\n", &(params->accel));

  // if (retval != 1) die("could not read param file: accel", __LINE__, __FILE__);

  // retval = fscanf(fp, "%lf\n", &(params->omega));

  // if (retval != 1) die("could not read param file: omega", __LINE__, __FILE__);

  // /* and close up the file */
  // fclose(fp);

  /*
  ** Allocate memory.
  **
  ** Remember C is pass-by-value, so we need to
  ** pass pointers into the initialise function.
  **
  ** NB we are allocating a 1D array, so that the
  ** memory will be contiguous.  We still want to
  ** index this memory as if it were a (row major
  ** ordered) 2D array, however.  We will perform
  ** some arithmetic using the row and column
  ** coordinates, inside the square brackets, when
  ** we want to access elements of this array.
  **
  ** Note also that we are using a structure to
  ** hold an array of 'speeds'.  We will allocate
  ** a 1D array of these structs.
  */

  // /* main grid */
  // *cells_ptr = (t_speed*)malloc(sizeof(t_speed) * (params->ny * params->nx));

  // if (*cells_ptr == NULL) die("cannot allocate memory for cells", __LINE__, __FILE__);

  // /* 'helper' grid, used as scratch space */
  // *tmp_cells_ptr = (t_speed*)malloc(sizeof(t_speed) * (params->ny * params->nx));

  // if (*tmp_cells_ptr == NULL) die("cannot allocate memory for tmp_cells", __LINE__, __FILE__);

  // /* the map of obstacles */
  // *obstacles_ptr = malloc(sizeof(int) * (params->ny * params->nx));

  // if (*obstacles_ptr == NULL) die("cannot allocate column memory for obstacles", __LINE__, __FILE__);

  // /* initialise densities */
  // double w0 = params->density * 4.0 / 9.0;
  // double w1 = params->density      / 9.0;
  // double w2 = params->density      / 36.0;

  // for (int ii = 0; ii < params->ny; ii++)
  // {
  //   for (int jj = 0; jj < params->nx; jj++)
  //   {
  //     /* centre */
  //     (*cells_ptr)[ii * params->nx + jj].speeds[0] = w0;
  //     /* axis directions */
  //     (*cells_ptr)[ii * params->nx + jj].speeds[1] = w1;
  //     (*cells_ptr)[ii * params->nx + jj].speeds[2] = w1;
  //     (*cells_ptr)[ii * params->nx + jj].speeds[3] = w1;
  //     (*cells_ptr)[ii * params->nx + jj].speeds[4] = w1;
  //     /* diagonals */
  //     (*cells_ptr)[ii * params->nx + jj].speeds[5] = w2;
  //     (*cells_ptr)[ii * params->nx + jj].speeds[6] = w2;
  //     (*cells_ptr)[ii * params->nx + jj].speeds[7] = w2;
  //     (*cells_ptr)[ii * params->nx + jj].speeds[8] = w2;
  //   }
  // }

  // /* first set all cells in obstacle array to zero */
  // for (int ii = 0; ii < params->ny; ii++)
  // {
  //   for (int jj = 0; jj < params->nx; jj++)
  //   {
  //     (*obstacles_ptr)[ii * params->nx + jj] = 0;
  //   }
  // }

  // /* open the obstacle data file */
  // fp = fopen(obstaclefile, "r");

  // if (fp == NULL)
  // {
  //   sprintf(message, "could not open input obstacles file: %s", obstaclefile);
  //   die(message, __LINE__, __FILE__);
  // }

  // /* read-in the blocked cells list */
  // while ((retval = fscanf(fp, "%d %d %d\n", &xx, &yy, &blocked)) != EOF)
  // {
  //   /* some checks */
  //   if (retval != 3) die("expected 3 values per line in obstacle file", __LINE__, __FILE__);

  //   if (xx < 0 || xx > params->nx - 1) die("obstacle x-coord out of range", __LINE__, __FILE__);

  //   if (yy < 0 || yy > params->ny - 1) die("obstacle y-coord out of range", __LINE__, __FILE__);

  //   if (blocked != 1) die("obstacle blocked value should be 1", __LINE__, __FILE__);

  //   /* assign to array */
  //   (*obstacles_ptr)[yy * params->nx + xx] = blocked;
  // }

  // /* and close the file */
  // fclose(fp);

  // /*
  // ** allocate space to hold a record of the avarage velocities computed
  // ** at each timestep
  // */
  // *av_vels_ptr = (double*)malloc(sizeof(double) * params->maxIters);
}