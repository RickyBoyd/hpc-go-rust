package main

import (
	"fmt"
	"os"
    "math"
    "syscall"
)

const NSPEEDS int = 9
const FINALSTATEFILE string =  "final_state.dat"
const AVVELSFILE string = "av_vels.dat"

/* struct to hold the parameter values */
type t_param struct {
  nx int            /* no. of cells in x-direction */
  ny int            /* no. of cells in y-direction */
  maxIters int      /* no. of iterations */
  reynolds_dim int  /* dimension for Reynolds number */
  density float64   /* density per link */
  accel float64     /* density redistribution */
  omega float64     /* relaxation parameter */
}

/* struct to hold the 'speed' values */
type t_speed struct {
  speeds [NSPEEDS]float64;
}

func main() {

    var paramfile string         /* name of the input parameter file */
    var obstaclefile string      /* name of a the input obstacle file */
    var params t_param           /* struct to hold parameter values */
    var cells []t_speed          /* grid containing fluid densities */
    var tmp_cells []t_speed      /* scratch space */
    var obstacles []int32        /* grid indicating which cells are blocked */
    var av_vels []float64        /* a record of the av. velocity computed for each timestep*/ 
    var timstr syscall.Timeval           /* structure to hold elapsed time */
    var ru syscall.Rusage                /* structure to hold CPU time--system and user */
    var usrtim float64           /* floating point number to record elapsed user CPU time */
    var systim float64           /* floating point number to record elapsed system CPU time */

    /* parse the command line */
    args := os.Args
    if len(args) != 3 {
        usage(args[0]);
    } else {
        paramfile = args[1];
        obstaclefile = args[2];
    }
    fmt.Println("Paramfile is", paramfile)
    fmt.Println("Obstacle file is", obstaclefile)
    /* initialise our data structures and load values from file */
    params = initialise(paramfile, obstaclefile, &cells, &tmp_cells, &obstacles, &av_vels);

    /* iterate for maxIters timesteps */
    syscall.Gettimeofday(&timstr);
    tic := float64(timstr.Sec) + (float64(timstr.Usec) / 1000000.0);

    for tt := 0; tt < params.maxIters; tt++ {
        timestep(params, cells, tmp_cells, obstacles);
        av_vels[tt] = av_velocity(params, cells, obstacles);
        //fmt.Printf("==timestep: %d==\n", tt);
        //fmt.Printf("av velocity: %.12E\n", av_vels[tt]);
        //fmt.Printf("tot density: %.12E\n", total_density(params, cells));
    }
    syscall.Gettimeofday(&timstr);
    toc := float64(timstr.Sec) + (float64(timstr.Usec) / 1000000.0);
    syscall.Getrusage(syscall.RUSAGE_SELF, &ru);
    timstr = ru.Utime;
    usrtim = float64(timstr.Sec) + (float64(timstr.Usec) / 1000000.0);
    timstr = ru.Stime;
    systim = float64(timstr.Sec) + (float64(timstr.Usec) / 1000000.0);
    fmt.Println("==done==");
    fmt.Printf("Reynolds number:\t\t%.12E\n", calc_reynolds(params, cells, obstacles));
    fmt.Printf("Elapsed time:\t\t\t%.6f (s)\n", toc - tic);
    fmt.Printf("Elapsed user CPU time:\t\t%.6f (s)\n", usrtim);
    fmt.Printf("Elapsed system CPU time:\t%.6f (s)\n", systim);
    write_values(params, cells, obstacles, av_vels);
}

func timestep(params t_param, cells []t_speed, tmp_cells []t_speed, obstacles []int32) {
    accelerate_flow(params, cells, obstacles);
    propagate(params, cells, tmp_cells);
    rebound(params, cells, tmp_cells, obstacles);
    collision(params, cells, tmp_cells, obstacles);
}

func accelerate_flow(params t_param, cells []t_speed, obstacles []int32) {
    /* compute weighting factors */
    w1 := params.density * params.accel / 9.0;
    w2 := params.density * params.accel / 36.0;

    /* modify the 2nd row of the grid */
    ii := params.ny - 2;

    for jj := 0; jj < params.nx; jj++ {
        /* if the cell is not occupied and
        ** we don't send a negative density */
        if (obstacles[ii * params.nx + jj] == 0 &&
            (cells[ii * params.nx + jj].speeds[3] - w1) > 0.0 && 
            (cells[ii * params.nx + jj].speeds[6] - w2) > 0.0 && 
            (cells[ii * params.nx + jj].speeds[7] - w2) > 0.0) {
            
            /* increase 'east-side' densities */
            cells[ii * params.nx + jj].speeds[1] += w1;
            cells[ii * params.nx + jj].speeds[5] += w2;
            cells[ii * params.nx + jj].speeds[8] += w2;
            /* decrease 'west-side' densities */
            cells[ii * params.nx + jj].speeds[3] -= w1;
            cells[ii * params.nx + jj].speeds[6] -= w2;
            cells[ii * params.nx + jj].speeds[7] -= w2;
        }
    }

}

func propagate(params t_param, cells []t_speed, tmp_cells []t_speed) {
    /* loop over _all_ cells */
    for ii := 0; ii < params.ny; ii++ {
        for jj := 0; jj < params.nx; jj++ {
            /* determine indices of axis-direction neighbours
            ** respecting periodic boundary conditions (wrap around) */
            y_n := (ii + 1) % params.ny;
            x_e := (jj + 1) % params.nx;
            var y_s int
            if y_s = ii - 1; ii == 0 {
                y_s = ii+params.ny-1
            }
            var x_w int
            if x_w = jj - 1; jj == 0 {
                x_w = jj + params.nx - 1
            }
            /* propagate densities to neighbouring cells, following
            ** appropriate directions of travel and writing into
            ** scratch space grid */
            tmp_cells[ii * params.nx + jj].speeds[0]  = cells[ii * params.nx + jj].speeds[0]; /* central cell, no movement */
            tmp_cells[ii * params.nx + x_e].speeds[1] = cells[ii * params.nx + jj].speeds[1]; /* east */
            tmp_cells[y_n * params.nx + jj].speeds[2]  = cells[ii * params.nx + jj].speeds[2]; /* north */
            tmp_cells[ii * params.nx + x_w].speeds[3] = cells[ii * params.nx + jj].speeds[3]; /* west */
            tmp_cells[y_s * params.nx + jj].speeds[4]  = cells[ii * params.nx + jj].speeds[4]; /* south */
            tmp_cells[y_n * params.nx + x_e].speeds[5] = cells[ii * params.nx + jj].speeds[5]; /* north-east */
            tmp_cells[y_n * params.nx + x_w].speeds[6] = cells[ii * params.nx + jj].speeds[6]; /* north-west */
            tmp_cells[y_s * params.nx + x_w].speeds[7] = cells[ii * params.nx + jj].speeds[7]; /* south-west */
            tmp_cells[y_s * params.nx + x_e].speeds[8] = cells[ii * params.nx + jj].speeds[8]; /* south-east */
        }
    }
}

func rebound(params t_param, cells []t_speed, tmp_cells []t_speed, obstacles []int32) {
    /* loop over the cells in the grid */
    for ii := 0; ii < params.ny; ii++ {
        for jj := 0; jj < params.nx; jj++ {
            /* if the cell contains an obstacle */
            if obstacles[ii * params.nx + jj] != 0 {
                /* called after propagate, so taking values from scratch space
                ** mirroring, and writing into main grid */
                cells[ii * params.nx + jj].speeds[1] = tmp_cells[ii * params.nx + jj].speeds[3];
                cells[ii * params.nx + jj].speeds[2] = tmp_cells[ii * params.nx + jj].speeds[4];
                cells[ii * params.nx + jj].speeds[3] = tmp_cells[ii * params.nx + jj].speeds[1];
                cells[ii * params.nx + jj].speeds[4] = tmp_cells[ii * params.nx + jj].speeds[2];
                cells[ii * params.nx + jj].speeds[5] = tmp_cells[ii * params.nx + jj].speeds[7];
                cells[ii * params.nx + jj].speeds[6] = tmp_cells[ii * params.nx + jj].speeds[8];
                cells[ii * params.nx + jj].speeds[7] = tmp_cells[ii * params.nx + jj].speeds[5];
                cells[ii * params.nx + jj].speeds[8] = tmp_cells[ii * params.nx + jj].speeds[6];
            }
        }
    }
}

func collision(params t_param, cells []t_speed, tmp_cells []t_speed, obstacles []int32) {
    c_sq := 1.0 / 3.0; /* square of speed of sound */
    w0   := 4.0 / 9.0;  /* weighting factor */
    w1   := 1.0 / 9.0;  /* weighting factor */
    w2   := 1.0 / 36.0; /* weighting factor */

    /* loop over the cells in the grid
    ** NB the collision step is called after
    ** the propagate step and so values of interest
    ** are in the scratch-space grid */
    for ii := 0; ii < params.ny; ii++ {
        for jj := 0; jj < params.nx; jj++ {
            /* don't consider occupied cells */
            if obstacles[ii * params.nx + jj] == 0 {
                /* compute local density total */
                local_density := 0.0;

                for kk := 0; kk < NSPEEDS; kk++ {
                    local_density += tmp_cells[ii * params.nx + jj].speeds[kk];
                }

                /* compute x velocity component */
                u_x := (tmp_cells[ii * params.nx + jj].speeds[1] + 
                        tmp_cells[ii * params.nx + jj].speeds[5] + 
                        tmp_cells[ii * params.nx + jj].speeds[8] - 
                        (tmp_cells[ii * params.nx + jj].speeds[3] + 
                        tmp_cells[ii * params.nx + jj].speeds[6] + 
                        tmp_cells[ii * params.nx + jj].speeds[7])) / local_density;
                /* compute y velocity component */
                u_y := (tmp_cells[ii * params.nx + jj].speeds[2] + 
                        tmp_cells[ii * params.nx + jj].speeds[5] + 
                        tmp_cells[ii * params.nx + jj].speeds[6] - 
                       (tmp_cells[ii * params.nx + jj].speeds[4] + 
                        tmp_cells[ii * params.nx + jj].speeds[7] + 
                        tmp_cells[ii * params.nx + jj].speeds[8])) / local_density;

                /* velocity squared */
                u_sq := u_x * u_x + u_y * u_y

                /* directional velocity components */
                var u [NSPEEDS]float64
                u[1] =   u_x        /* east */
                u[2] =         u_y  /* north */
                u[3] = - u_x;        /* west */
                u[4] =       - u_y  /* south */
                u[5] =   u_x + u_y  /* north-east */
                u[6] = - u_x + u_y  /* north-west */
                u[7] = - u_x - u_y  /* south-west */
                u[8] =   u_x - u_y  /* south-east */

                /* equilibrium densities */
                var d_equ [NSPEEDS]float64
                /* zero velocity density: weight w0 */
                d_equ[0] = w0 * local_density * (1.0 - u_sq / (2.0 * c_sq))
                /* axis speeds: weight w1 */
                d_equ[1] = w1 * local_density * (1.0 + u[1] / c_sq  + 
                                (u[1] * u[1]) / (2.0 * c_sq * c_sq) - 
                                        u_sq / (2.0 * c_sq));
                d_equ[2] = w1 * local_density * (1.0 + u[2] / c_sq + 
                                (u[2] * u[2]) / (2.0 * c_sq * c_sq) -
                                         u_sq / (2.0 * c_sq));
                d_equ[3] = w1 * local_density * (1.0 + u[3] / c_sq + 
                                (u[3] * u[3]) / (2.0 * c_sq * c_sq) - 
                                        u_sq / (2.0 * c_sq));
                d_equ[4] = w1 * local_density * (1.0 + u[4] / c_sq + 
                                (u[4] * u[4]) / (2.0 * c_sq * c_sq) - 
                                        u_sq / (2.0 * c_sq));
                /* diagonal speeds: weight w2 */
                d_equ[5] = w2 * local_density * (1.0 + u[5] / c_sq + 
                                (u[5] * u[5]) / (2.0 * c_sq * c_sq) - 
                                        u_sq / (2.0 * c_sq));
                d_equ[6] = w2 * local_density * (1.0 + u[6] / c_sq + 
                                (u[6] * u[6]) / (2.0 * c_sq * c_sq) - 
                                        u_sq / (2.0 * c_sq));
                d_equ[7] = w2 * local_density * (1.0 + u[7] / c_sq + 
                                (u[7] * u[7]) / (2.0 * c_sq * c_sq) - 
                                        u_sq / (2.0 * c_sq));
                d_equ[8] = w2 * local_density * (1.0 + u[8] / c_sq + 
                                (u[8] * u[8]) / (2.0 * c_sq * c_sq) - 
                                        u_sq / (2.0 * c_sq));

                /* relaxation step */
                for kk := 0; kk < NSPEEDS; kk++ {
                    cells[ii * params.nx + jj].speeds[kk] = tmp_cells[ii * params.nx + jj].speeds[kk] + 
                                                params.omega * (d_equ[kk] - tmp_cells[ii * params.nx + jj].speeds[kk]);
                }
            }
        }
    }
}


func calc_reynolds( params t_param, cells []t_speed, obstacles []int32 ) float64 {
    viscosity := 1.0 / 6.0 * (2.0 / params.omega - 1.0)

    return av_velocity(params, cells, obstacles) * float64(params.reynolds_dim) / viscosity;
}

func av_velocity( params t_param, cells []t_speed, obstacles []int32 ) float64 {
    tot_cells := 0;  /* no. of cells used in calculation */
    /* initialise */
    tot_u := 0.0;

    /* loop over all non-blocked cells */
    for ii := 0; ii < params.ny; ii++ {
        for jj := 0; jj < params.nx; jj++ {
            /* ignore occupied cells */
            if obstacles[ii * params.nx + jj] == 0 {
                /* local density total */
                local_density := 0.0;

                for kk := 0; kk < NSPEEDS; kk++ {
                    local_density += cells[ii * params.nx + jj].speeds[kk];
                }

                /* x-component of velocity */
                u_x := (cells[ii * params.nx + jj].speeds[1] + cells[ii * params.nx + jj].speeds[5] + 
                        cells[ii * params.nx + jj].speeds[8] - (cells[ii * params.nx + jj].speeds[3] + 
                        cells[ii * params.nx + jj].speeds[6] + cells[ii * params.nx + jj].speeds[7])) / local_density;
                /* compute y velocity component */
                u_y := (cells[ii * params.nx + jj].speeds[2] + cells[ii * params.nx + jj].speeds[5] + 
                        cells[ii * params.nx + jj].speeds[6] - (cells[ii * params.nx + jj].speeds[4] + 
                        cells[ii * params.nx + jj].speeds[7] + cells[ii * params.nx + jj].speeds[8])) / local_density;
                /* accumulate the norm of x- and y- velocity components */
                tot_u += math.Sqrt((u_x * u_x) + (u_y * u_y));
                /* increase counter of inspected cells */
                tot_cells += 1
            }
        }
    }
    return tot_u / float64(tot_cells)
}

func total_density(params t_param, cells []t_speed) float64 {
    total := 0.0;  /* accumulator */

    for ii := 0; ii < params.ny; ii++ {
        for jj := 0; jj < params.nx; jj++ {
            for kk := 0; kk < NSPEEDS; kk++ {
                total += cells[ii * params.nx + jj].speeds[kk];
            }
        }
    }
    return total
}

func initialise(paramfile string, obstaclefile string, cells_ptr *[]t_speed, tmp_cells_ptr *[]t_speed, obstacles_ptr *[]int32, av_vels_ptr *[]float64) t_param {
    var message string     /* message buffer */
    var xx, yy int         /* generic array indices */
    var blocked int32      /* indicates whether a cell is blocked by an obstacle */
    var retval int
    
    params := t_param{}

    fp, err := os.Open(paramfile)
    if err != nil {
        fmt.Sprintf(message, "could not open input parameter file: %s", paramfile);
        panic(message);
    }
    defer fp.Close()

    /* read in the parameter values */
    _, err = fmt.Fscanf(fp, "%d\n", &(params.nx))
    if err != nil {
        panic("could not read param file: nx")
    } 

    _, err = fmt.Fscanf(fp, "%d\n", &(params.ny))
    if err != nil {
        panic("could not read param file: ny")
    }

    _, err = fmt.Fscanf(fp, "%d\n", &(params.maxIters))
    if err != nil {
        panic("could not read param file: maxIters")
    }

    _, err = fmt.Fscanf(fp, "%d\n", &(params.reynolds_dim))
    if err != nil {
        panic("could not read param file: reynolds_dim")
    }

    _, err = fmt.Fscanf(fp, "%f\n", &(params.density))
    if err != nil {
        panic("could not read param file: density")
    }

    _, err = fmt.Fscanf(fp, "%f\n", &(params.accel))
    if err != nil {
        panic("could not read param file: accel")
    }

    _, err = fmt.Fscanf(fp, "%f\n", &(params.omega))
    if err != nil {
        panic("could not read param file: omega")
    }

    /* main grid */
    (*cells_ptr) = make([]t_speed, params.ny * params.nx, params.ny * params.nx)

    if (*cells_ptr) == nil {
        panic("cannot allocate memory for cells")
    }

    /* 'helper' grid, used as scratch space */
    (*tmp_cells_ptr) = make([]t_speed, params.ny * params.nx, params.ny * params.nx)

    if (*tmp_cells_ptr) == nil {
        panic("cannot allocate memory for tmp_cells");
    }

    /* the map of obstacles */
    (*obstacles_ptr) = make([]int32, params.nx*params.ny, params.nx * params.ny)

    if (*obstacles_ptr) == nil{
        panic("cannot allocate column memory for obstacles");
    } 

    /* initialise densities */
    w0 := params.density * 4.0 / 9.0;
    w1 := params.density      / 9.0;
    w2 := params.density      / 36.0;

    for ii := 0; ii < params.ny; ii++ {
        for jj := 0; jj < params.nx; jj++ {
            /* centre */
            (*cells_ptr)[ii * params.nx + jj].speeds[0] = w0;
            /* axis directions */
            (*cells_ptr)[ii * params.nx + jj].speeds[1] = w1;
            (*cells_ptr)[ii * params.nx + jj].speeds[2] = w1;
            (*cells_ptr)[ii * params.nx + jj].speeds[3] = w1;
            (*cells_ptr)[ii * params.nx + jj].speeds[4] = w1;
            /* diagonals */
            (*cells_ptr)[ii * params.nx + jj].speeds[5] = w2;
            (*cells_ptr)[ii * params.nx + jj].speeds[6] = w2;
            (*cells_ptr)[ii * params.nx + jj].speeds[7] = w2;
            (*cells_ptr)[ii * params.nx + jj].speeds[8] = w2;
        }
    }

    /* first set all cells in obstacle array to zero */
    for ii := 0; ii < params.ny; ii++ {
        for jj := 0; jj < params.nx; jj++ {
            (*obstacles_ptr)[ii * params.nx + jj] = 0;
        }
    }

    /* open the obstacle data file */
    fp_obstacles, err := os.Open(obstaclefile)

    if (err != nil) {
        fmt.Sprintf(message, "could not open input obstacles file: %s", obstaclefile)
        panic(message)
    }

    defer fp_obstacles.Close()

    /* read-in the blocked cells list */
    
    for err == nil {
        /* some checks */
        retval, err = fmt.Fscanf(fp_obstacles, "%d %d %d\n", &xx, &yy, &blocked)
        if err != nil { break; }
        if retval != 3 {
            panic("expected 3 values per line in obstacle file")
        }

        if xx < 0 || xx > params.nx - 1 {
            panic("obstacle x-coord out of range")
        }

        if yy < 0 || yy > params.ny - 1 {
            panic("obstacle y-coord out of range")
        }

        if blocked != 1 {
            panic("obstacle blocked value should be 1")
        }

        /* assign to array */
        (*obstacles_ptr)[yy * params.nx + xx] = blocked;
    }


    /*
    ** allocate space to hold a record of the avarage velocities computed
    ** at each timestep
    */
    (*av_vels_ptr) = make([]float64, params.maxIters, params.maxIters)

    return params
}

func write_values(params t_param, cells []t_speed, obstacles []int32, av_vels []float64) {
    const c_sq float64 = 1.0 / 3.0    /* sq. of speed of sound */
    var local_density float64         /* per grid cell sum of densities */
    var pressure float64              /* fluid pressure in grid cell */
    var u_x float64                   /* x-component of velocity in grid cell */
    var u_y float64                   /* y-component of velocity in grid cell */
    var u float64                     /* norm--root of summed squares--of u_x and u_y */

    fp, err := os.OpenFile(FINALSTATEFILE, os.O_WRONLY|os.O_CREATE, 0666)


    if err != nil {
        panic("could not open file output file")
    }
    defer fp.Close()

    for ii := 0; ii < params.ny; ii++ {
        for jj := 0; jj < params.nx; jj++ {
            /* an occupied cell */
            if obstacles[ii * params.nx + jj] != 0 {
                u_x, u_y, u = 0.0, 0.0, 0.0
                pressure = params.density * c_sq;
            } else {
                local_density = 0.0

                for kk := 0; kk < NSPEEDS; kk++ {
                    local_density += cells[ii * params.nx + jj].speeds[kk]
                }

                /* compute x velocity component */
                u_x = (cells[ii * params.nx + jj].speeds[1] + cells[ii * params.nx + jj].speeds[5] + 
                       cells[ii * params.nx + jj].speeds[8] - (cells[ii * params.nx + jj].speeds[3] + 
                       cells[ii * params.nx + jj].speeds[6] + cells[ii * params.nx + jj].speeds[7])) / local_density;
                /* compute y velocity component */
                u_y = (cells[ii * params.nx + jj].speeds[2] + cells[ii * params.nx + jj].speeds[5] + 
                       cells[ii * params.nx + jj].speeds[6] - (cells[ii * params.nx + jj].speeds[4] + 
                       cells[ii * params.nx + jj].speeds[7] + cells[ii * params.nx + jj].speeds[8])) / local_density;
                /* compute norm of velocity */
                u = math.Sqrt((u_x * u_x) + (u_y * u_y))
                /* compute pressure */
                pressure = local_density * c_sq
            }

            /* write to file */
            fmt.Fprintf(fp, "%d %d %.12E %.12E %.12E %.12E %d\n", jj, ii, u_x, u_y, u, pressure, obstacles[ii * params.nx + jj]);
        }
    }


    fp_vels, err_vels := os.OpenFile(AVVELSFILE, os.O_WRONLY|os.O_CREATE, 0666)

    if err_vels != nil {
        panic("could not open file output file")
    }

    defer fp_vels.Close()

    for ii := 0; ii < params.maxIters; ii++ {
        fmt.Fprintf(fp_vels, "%d:\t%.12E\n", ii, av_vels[ii])
    }
}

func usage(exe string) {
    fmt.Printf("Usage: %s <paramfile> <obstaclefile>\n", exe)
    os.Exit(1)
}