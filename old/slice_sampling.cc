#include <climits>
#include <numeric>
#include <string.h>

#include <gsl/gsl_math.h>
#include "slice_sampling.h"

#include "tools.h"
#include "hilbert.h"
#include "sampling.h"
#include "error_estimate.h"
#include "defs.h"

SliceSampling::SliceSampling(size_t _range_delta) : 
    range_delta(_range_delta),
    initialized(false)
{
}


SliceSampling::~SliceSampling()
{
}


void SliceSampling::Initialize()
{
    this->total_bits = NDIM * NBITS_PER_DIM;

    mpz_init(this->U);
    mpz_init(this->N);
    mpz_init(this->B);
    mpz_init(this->zero);
    mpz_init(this->xprime);
    mpf_init(this->uniform);

    mpz_set_ui(this->B, 0);
    mpz_set_ui(this->zero, 0);
    mpz_set_d(this->B, pow(2,this->total_bits));

    gmp_randinit_default(this->rand_state);

    //size_t const bits = sizeof(double) * 8;
    this->initialized = true;
}

#define GRID_FACTOR (double)(((uint64_t)(1)<<(NBITS_PER_DIM+1)) - 1)

//linearly scale ndim [0,1] (double) cube coords into (uint64_t)
//[0,max_64bit_int] grid_coords
void expand_to_grid_coords(const double *real_coord,
                           uint64_t *grid_coord)
{
    for (int i = 0; i != NDIM; ++i)
        grid_coord[i] = (uint64_t)(real_coord[i] * GRID_FACTOR);
    
}

void contract_from_grid_coords(const uint64_t *grid_coord,
                               double *real_coord)
{
    real_coord[NDIM] = 1;
    for (int i = 0; i != NDIM; ++i)
        real_coord[NDIM] -= real_coord[i] = (double)grid_coord[i] / GRID_FACTOR;
}


//randomly step in the integer grid with a step size up to
//<current_range> bits
void SliceSampling::grid_step(int current_range, 
                              const uint64_t *xg,
                              uint64_t *xgprime)
{

    mpz_t x;
    mpz_init(x);

    uint64_t xgu[NDIM], xgu_prime[NDIM];

    Hilbert_to_int(xg, x);

    mpz_urandomb(this->N, this->rand_state, current_range);
    mpz_sub(x, x, this->U);
    if (mpz_cmp(x, this->zero) < 0)
    {
        //xprime < 0: wrap it
        mpz_add(x, x, this->B);
        int cmp = mpz_cmp(x, this->zero);
        assert(cmp > 0);
    }
    int_to_Hilbert(x, xgu);

    mpz_xor(x, x, this->N);
    int_to_Hilbert(x, xgu_prime);

//     mpz_add(*xnew, *xnew, this->U);
//     if (mpz_cmp(*xnew, this->B) > 0)
//     {
//         //xprime > this->B:  wrap it
//         mpz_sub(*xnew, *xnew, this->B);
//         int cmp = mpz_cmp(*xnew, this->B);
//         assert(cmp < 0);
//     }

    for (size_t d = 0; d != NDIM; ++d)
        xgprime[d] = xg[d] + xgu_prime[d] - xgu[d];

}


int in_simplex(const double *x)
{
    int i;
    for (i = 0; i != NUM_NUCS; ++i)
        if (x[i] < 0) return 0;

    double sum = 0;
    for (i = 0; i != NUM_NUCS; ++i)
        sum += x[i];

    return gsl_fcmp(sum, 1.0, 1e-10) == 0;
}

/* sample from the search interval until a point is within the slice
   (integrand value is above y).  at each point not within the slice,
   shrink the interval so as to retain the initial point x within it.
   interval starts as a neighborhood of x with initial_range bits of
   width.  updates xcoord, xcoord_grid and yprime with latest values
*/
int SliceSampling::step_in(ErrorEstimate *posterior,
                           const uint64_t *xg,
                           double y,
                           int initial_range,
                           uint64_t *xgp)
{

    int current_range = initial_range + this->range_delta;
    mpz_urandomb(this->U, this->rand_state, this->total_bits);
    double xrp[NUM_NUCS], yp;

    do
    {
        current_range -= this->range_delta;
        this->grid_step(current_range, xg, xgp);
        contract_from_grid_coords(xgp, xrp);
        yp = in_simplex(xrp) ? posterior->log_pdf_trunc(xrp) : -DBL_MAX;
    }
    while (yp < y && (! std::equal(xg, xg + NDIM, xgp)));
    
    return current_range;
}



/* expand the search interval until a test point is found that lies
   outside (below) the slice defined by level y.  start in the
   neighborhood of x with a neighborhood of intitial_range
   bits. returns the range found.
*/
int SliceSampling::step_out(ErrorEstimate *posterior, 
                            const uint64_t *xg,
                            double y,
                            int initial_range)
{

    int current_range = initial_range - this->range_delta;
    mpz_urandomb(this->U, this->rand_state, this->total_bits);
    uint64_t xgp[NDIM];
    double xrp[NUM_NUCS], yp;

    do
    {
        current_range += this->range_delta;

        this->grid_step(current_range, xg, xgp);
        contract_from_grid_coords(xgp, xrp);
        yp = in_simplex(xrp) ? posterior->log_pdf_trunc(xrp) : -DBL_MAX;
    }
    while (current_range != (int)this->total_bits && yp >= y);

    return current_range;
}


/* selects an auxiliary coordinate between U(0, f(x)).  The auxiliary
   coordinate is also sampled in log space, and the integrand given is
   assumed to be log(integrand-of-interest)
   !!! untested code
 */
double SliceSampling::choose_auxiliary_coord(ErrorEstimate *posterior,
                                             const double *x)
{

    mpf_urandomb(this->uniform, this->rand_state, 64);
    double y = posterior->log_pdf_trunc(x) + log(mpf_get_d(this->uniform));
    return y;
}


/* xprime and yprime are the proposed new point for slice sampling they
   are only accepted if yprime < y.  once accepted, the slice sampling
   continues, using (xprime, yprime) as the new (x, y)

   In order to conform to sampling::print_quantiles, sample_points_flat
   shall be populated with 4-tuples of normalized points, rather than 3.
   The 4th is just auto-filled in
*/
void SliceSampling::sample(ErrorEstimate *posterior,
                           const double *starting_x,
                           int initial_range,
                           size_t every_nth,
                           double *sample_points_flat,
                           size_t num_samples)
{

    if (! this->initialized)
    {
        fprintf(stderr, "Must first initialize SliceSampling\n");
        exit(1);
    }

    if (static_cast<size_t>(initial_range) > this->total_bits)
    {
        fprintf(stderr, "initial_range (%i) must be <= total bits (%Zu)\n",
                initial_range, this->total_bits);
        exit(1);
    }

    // grid coordinate coorresponding to xrp
    uint64_t xg[NDIM];

    // grid coordinate representing 'previous' xg, in the markov process
    uint64_t xgp[NDIM];

    // coordinate in real space, that the posterior is actually evaluated in
    double xrp[NUM_NUCS];

    //initialize
    int current_range = initial_range;
    memcpy(xrp, starting_x, sizeof(xrp));

    expand_to_grid_coords(xrp, xg);
    double y = this->choose_auxiliary_coord(posterior, xrp);

    double *point = sample_points_flat;
    size_t si, se = num_samples * every_nth;
    for (si = 0; si != se; ++si)
    {
        current_range = initial_range;
        current_range = this->step_out(posterior, xg, y, current_range);
        current_range = this->step_in(posterior, xg, y, current_range, xgp);
        contract_from_grid_coords(xgp, xrp);

        if (si % every_nth == 0)
        {
            memcpy(point, xrp, sizeof(xrp));
            point += NUM_NUCS;
        }

        //update markov chain
        memcpy(xg, xgp, sizeof(xg));

        //choose y coordinate from new x coordinate
        y = this->choose_auxiliary_coord(posterior, xrp);
        
    }
}


/*
void expand_to_grid_coords(double const* cube_coord, size_t const ndim, 
                           uint64_t * grid_coord, size_t const nbits_per_dim)
{
    //rescale the number from [0,1] to [0, max_integer]
    double factor = static_cast<double>((static_cast<uint64_t>(1)<<(nbits_per_dim+1)) - 1);

    for (size_t i = 0; i != ndim; ++i)
        grid_coord[i] = static_cast<uint64_t>(cube_coord[i] * factor);
}


//linearly scale ndim (uint64_t) [0,max_64bit_int] grid_coords into
//[0,1] (double) cube coords
void contract_from_grid_coords(const uint64_t *grid_coord, size_t const ndim,
                               double * cube_coord, size_t const nbits_per_dim)
{
    double factor = static_cast<double>((static_cast<uint64_t>(1)<<(nbits_per_dim+1)) - 1);

    for (size_t i = 0; i != ndim; ++i)
        cube_coord[i] = static_cast<double>(grid_coord[i]) / factor;
}
*/



/* initialize xcoord, xcoord_grid (the projection of xcoord onto an
   integer grid), and xprime (xcoord_grid mapped onto the hilbert
   curve) from starting_x (plain unit_hypercube x coordinates)
 */
/*
void SliceSampling::initialize_starting_point(double const* starting_x, 
                                              size_t const ndim)
{

    mpf_urandomb(this->uniform, this->rand_state, 61);

    std::copy(starting_x, starting_x + ndim, this->xcoord);
    
    expand_to_grid_coords(this->xcoord, this->xcoord_grid);

    Hilbert_to_int(this->xcoord_grid, ndim, NBITS_PER_DIM, this->xprime);
    
}    
*/


/* sample a bits interval of current_range in the neighborhood of x,
   storing the result in xprime
 */
/*
void SliceSampling::sample_bits_interval(int const current_range, 
                                         mpz_t const& x,
                                         uint64_t * xg,
                                         mpz_t * xnew,
                                         uint64_t * xgnew)
{

    mpz_urandomb(this->N, this->rand_state, current_range);
        
    mpz_sub(*xnew, x, this->U);
        
    if (mpz_cmp(*xnew, this->zero) < 0)
    {
        //xprime < 0: wrap it
        mpz_add(*xnew, *xnew, this->B);
        int cmp = mpz_cmp(*xnew, this->zero);
        assert(cmp > 0);
    }
    int_to_Hilbert(*xnew, NBITS_PER_DIM, xgnew, NDIM);

    mpz_xor(*xnew, *xnew, this->N);
    mpz_add(*xnew, *xnew, this->U);
    if (mpz_cmp(*xnew, this->B) > 0)
    {
        //xprime > this->B:  wrap it
        mpz_sub(*xnew, *xnew, this->B);
        int cmp = mpz_cmp(*xnew, this->B);
        assert(cmp < 0);
    }

}
*/
