#include <iostream>

// ViennaCL includes
#include "viennacl/compressed_matrix.hpp"
#include "viennacl/io/matrix_market.hpp"
#include "viennacl/linalg/amg.hpp"
#include "viennacl/linalg/bicgstab.hpp"
#include "viennacl/linalg/cg.hpp"
#include "viennacl/linalg/gmres.hpp"
#include "viennacl/linalg/jacobi_precond.hpp"
#include "viennacl/linalg/prod.hpp"
#include "viennacl/scalar.hpp"
#include "viennacl/vector.hpp"

extern "C" {
#include "matrix_csr.h"
}

#ifdef DEBUG
#define DEBUG_PRINT(s) s
#else
#define DEBUG_PRINT(s)
#endif

typedef double ScalarType;

struct viennacl_solver_ctx {
    int use_amg;
    idx_t N;
    double rtol = 1E-5;
    viennacl::compressed_matrix<ScalarType> A;
    viennacl::vector<ScalarType> b;
    viennacl::vector<ScalarType> x_guess;
    viennacl::linalg::cg_solver<viennacl::vector<ScalarType>> *cg = nullptr;
    viennacl::linalg::amg_precond<viennacl::compressed_matrix<ScalarType>> *amg = nullptr;

    viennacl_solver_ctx(struct matrix_csr A_csr, int use_amg) {
        N = A_csr.num_rows;
        this->use_amg = use_amg;
        this->rtol = 1E-5;

        // set up vectors
        this->b = viennacl::vector<ScalarType>(N);
        this->x_guess = viennacl::vector<ScalarType>(N);
        this->x_guess.clear(); // set all entries to zero
        DEBUG_PRINT(printf("Set up vectors\n"));

        A.set(
            A_csr.row_ptr,
            A_csr.col,
            A_csr.val,
            A_csr.num_rows,
            A_csr.num_cols,
            A_csr.num_nonzeros
        );

        // Initialise AMG solver
        {
            DEBUG_PRINT(printf("Initialising AMG preconditioner\n"));

            viennacl::linalg::amg_tag my_amg_tag;
            //my_amg_tag.set_coarsening_method(viennacl::linalg::AMG_COARSENING_METHOD_ONEPASS);
            my_amg_tag.set_coarsening_method(viennacl::linalg::AMG_COARSENING_METHOD_MIS2_AGGREGATION);
            //my_amg_tag.set_interpolation_method(viennacl::linalg::AMG_INTERPOLATION_METHOD_DIRECT);
            my_amg_tag.set_interpolation_method(
                    viennacl::linalg::AMG_INTERPOLATION_METHOD_SMOOTHED_AGGREGATION);
            my_amg_tag.set_strong_connection_threshold(0.25);
            my_amg_tag.set_jacobi_weight(0.67);
            my_amg_tag.set_presmooth_steps(1);
            my_amg_tag.set_postsmooth_steps(1);

            DEBUG_PRINT(printf("Creating AMG precond\n"));
            this->amg = new viennacl::linalg::amg_precond<viennacl::compressed_matrix<ScalarType> >(A, my_amg_tag);

            if (use_amg) {
                printf("Setting up AMG precond ... ");
                fflush(stdout);
                amg->setup();
                printf("done\n");
            }
            DEBUG_PRINT(printf("Done with AMG precond\n"));
        }
    }

    ~viennacl_solver_ctx() {
        delete this->cg;
        if (this->amg != nullptr) {
            delete this->amg;
        }
    }
};
typedef struct viennacl_solver_ctx viennacl_solver_ctx_t;

template<typename VectorT, typename NumericT, typename MatrixT>
bool monitor_amg(VectorT const & current_approx, NumericT residual_estimate, void *user_data)
{
    // Extract residual:
    struct viennacl_solver_ctx const *ctx = reinterpret_cast<struct viennacl_solver_ctx *>(user_data);

    // Form residual r = b - A*x, taking an initial guess into account: r = b - A * (current_approx + x_initial)
    VectorT x = current_approx + ctx->x_guess;
    VectorT residual = ctx->b - viennacl::linalg::prod(ctx->A, x);

    double true_norm = viennacl::linalg::norm_2(residual) / viennacl::linalg::norm_2(ctx->b);
    //std::cout << "Residual estimate vs. true residual: " << residual_estimate << " vs. " << true_norm << std::endl;

    bool should_terminate = true_norm < ctx->rtol;
    return should_terminate;
}


extern "C"
{
int solver_viennacl_init(void **viennacl_solver_ctx_ptr,
    struct matrix_csr A_csr, int use_amg)
{
    assert(A_csr.num_rows == A_csr.num_cols);
    DEBUG_PRINT(printf("num_rows: %llu  num_rows: %llu  num_nonzeros: %llu\n",
        A_csr.num_rows,
        A_csr.num_cols,
        A_csr.num_nonzeros
    ));
#if 0
    if (A_csr.num_rows < 10) {
        for (int i = 0; i <= A_csr.num_rows; i++) {
            printf("row_ptr[%d] = %d\n", i, A_csr.row_ptr[i]);
        }
        for (int i = 0; i < A_csr.num_nonzeros; i++) {
            printf("col[%d] = %6d; ", i, A_csr.col[i]);
            printf("val[%d] = %6g\n", i, A_csr.val[i]);
        }
    }
#endif

    viennacl_solver_ctx_t *ctx = new viennacl_solver_ctx_t(A_csr, use_amg);

    *viennacl_solver_ctx_ptr = ctx;
    DEBUG_PRINT(printf("Initialising succeeded\n"));
    return EXIT_SUCCESS;
}

int solver_viennacl_solve(viennacl_solver_ctx_t *ctx, ScalarType *b, ScalarType *x, idx_t N)
{
    // Check that vector lengths matches matrix dimensions
    assert(ctx->A.size1() == N);

    // copy b array into b vector
    copy(b, b+N, ctx->b.begin());
    DEBUG_PRINT(printf("Copied b array into b vector\n"));

    // Set up CG tag
    double rtol = ctx->rtol;
    int max_iters = -1;

    /* when using AMG, we want to use the monitor to determine when to terminate,
    so we set rtol=0 and max_iters = -1 so that the solver does not terminate
    unless the monitor signals it to do so */
    if (ctx->use_amg) {
        rtol = 0;
    }

    viennacl::linalg::cg_tag my_cg_tag(rtol, max_iters);

    // configure termination criterion for CG without AMG
    /* the overhead of evaluating an extra sparse-matrix product in the monitor
    is substantial, so we instead set the absolute tolerance */
    if (!ctx->use_amg) {
        double atol = ctx->rtol * viennacl::linalg::norm_2(ctx->b);
        my_cg_tag.abs_tolerance(atol);
    }

    viennacl::linalg::cg_solver<viennacl::vector<ScalarType>> cg(my_cg_tag);

    // set initial guess to the previous solution
    cg.set_initial_guess(ctx->x_guess);

    if (!ctx->use_amg) {
        DEBUG_PRINT(printf("Solving without AMG\n"));
        // without preconditioner
        ctx->x_guess = cg(ctx->A, ctx->b);
        DEBUG_PRINT(printf("Solved without AMG\n"));
    } else {
        DEBUG_PRINT(printf("Solving with AMG\n"));
        assert(ctx->amg != nullptr);
        // with AMG preconditioner

        // set monitor to decide when to terminate
        cg.set_monitor(
            monitor_amg<viennacl::vector<ScalarType>, ScalarType, viennacl::compressed_matrix<ScalarType>>,
            ctx
        );
        ctx->x_guess = cg(ctx->A, ctx->b, *ctx->amg);
        DEBUG_PRINT(printf("Solved with AMG\n"));
    }

    // copy x vector back to array
    copy(ctx->x_guess.begin(), ctx->x_guess.end(), x);
    DEBUG_PRINT(printf("Copied x vector into x array\n"));

    return EXIT_SUCCESS;
}

int solver_viennacl_free(viennacl_solver_ctx_t *ctx)
{
    DEBUG_PRINT(printf("In solver_viennacl_free\n"));
    delete ctx;
    return EXIT_SUCCESS;
}
} // extern "C"


