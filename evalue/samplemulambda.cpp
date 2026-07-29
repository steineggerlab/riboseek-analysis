// samplemulambda <inputfile> <outputfile>
// inputfile is a text file with one float (score) per line
// outputfile is the output file to write the mu and lambda values to
#include <random>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <algorithm>
#include <fstream>

/* Function: Lawless416()
 * Date:     SRE, Thu Nov 13 11:48:50 1997 [St. Louis]
 *
 * Purpose:  Equation 4.1.6 from [Lawless82], pg. 143, and
 *           its first derivative with respect to lambda,
 *           for finding the ML fit to EVD lambda parameter.
 *           This equation gives a result of zero for the maximum
 *           likelihood lambda.
 *
 *           Can either deal with a histogram or an array.
 *
 *           Warning: beware overflow/underflow issues! not bulletproof.
 *
 * Args:     x      - array of sample values (or x-axis of a histogram)
 *           y      - NULL (or y-axis of a histogram)
 *           n      - number of samples (or number of histogram bins)
 *           lambda - a lambda to test
 *           ret_f  - RETURN: 4.1.6 evaluated at lambda
 *           ret_df - RETURN: first derivative of 4.1.6 evaluated at lambda
 *
 * Return:   (void)
 */
void
Lawless416(float *x, int *y, int n, float lambda, float *ret_f, float *ret_df)
{

    double esum;			/* \sum e^(-lambda xi)      */
    double xesum;			/* \sum xi e^(-lambda xi)   */
    double xxesum;		/* \sum xi^2 e^(-lambda xi) */
    double xsum;			/* \sum xi                  */
    double mult;			/* histogram count multiplier */
    double total;			/* total samples            */
    int i;


    esum = xesum = xsum  = xxesum = total = 0.;
    for (i = 0; i < n; i++)
    {
        mult = (y == NULL) ? 1. : (double) y[i];
        xsum   += mult * x[i];
        xesum  += mult * x[i] * exp(-1. * lambda * x[i]);
        xxesum += mult * x[i] * x[i] * exp(-1. * lambda * x[i]);
        esum   += mult * exp(-1. * lambda * x[i]);
        total  += mult;
    }
    *ret_f  = 1./lambda - xsum / total + xesum / esum;
    *ret_df = ((xesum / esum) * (xesum / esum))
              - (xxesum / esum)
              - (1. / (lambda * lambda));

    return;
}

/* Function: EVDMaxLikelyFit()
 * Date:     SRE, Fri Nov 14 07:56:29 1997 [St. Louis]
 *
 * Purpose:  Given a list or a histogram of EVD-distributed samples,
 *           find maximum likelihood parameters lambda and
 *           mu.
 *
 * Algorithm: Uses approach described in [Lawless82]. Solves
 *           for lambda using Newton/Raphson iterations;
 *           then substitutes lambda into Lawless' equation 4.1.5
 *           to get mu.
 *
 *           Newton/Raphson algorithm developed from description in
 *           Numerical Recipes in C [Press88].
 *
 * Args:     x          - list of EVD distributed samples or x-axis of histogram
 *           c          - NULL, or y-axis of histogram
 *           n          - number of samples, or number of histogram bins
 *           ret_mu     : RETURN: ML estimate of mu
 *           ret_lambda : RETURN: ML estimate of lambda
 *
 * Return:   1 on success; 0 on any failure
 */
 int
 EVDMaxLikelyFit(float *x, int *c, int n, float *ret_mu, float *ret_lambda)
 {
     float  lambda, mu;
     float  fx;			/* f(x)  */
     float  dfx;			/* f'(x) */
     double esum;                  /* \sum e^(-lambda xi) */
     double mult;
     double total;
     float  tol = 1e-5;
     int    i;
 
     /* 1. Find an initial guess at lambda: linear regression here?
      */
     lambda = 0.2;
 
     /* 2. Use Newton/Raphson to solve Lawless 4.1.6 and find ML lambda
      */
     for (i = 0; i < 100; i++)
     {
         Lawless416(x, c, n, lambda, &fx, &dfx);
         if (fabs(fx) < tol) break;             /* success */
         lambda = lambda - fx / dfx;	     /* Newton/Raphson is simple */
         if (lambda <= 0.) lambda = 0.001;      /* but be a little careful  */
     }
 
     /* 2.5: If we did 100 iterations but didn't converge, Newton/Raphson failed.
      *      Resort to a bisection search. Worse convergence speed
      *      but guaranteed to converge (unlike Newton/Raphson).
      *      We assume (!?) that fx is a monotonically decreasing function of x;
      *      i.e. fx > 0 if we are left of the root, fx < 0 if we
      *      are right of the root.
      */
     if (i == 100)
     {
         float left, right, mid;
         printf(("EVDMaxLikelyFit(): Newton/Raphson failed; switchover to bisection"));
 
         /* First we need to bracket the root */
         lambda = right = left = 0.2;
         Lawless416(x, c, n, lambda, &fx, &dfx);
         if (fx < 0.)
         {			/* fix right; search left. */
             do {
                 left -= 0.1;
                 if (left < 0.) {
                     printf(("EVDMaxLikelyFit(): failed to bracket root"));
                     return 0;
                 }
                 Lawless416(x, c, n, left, &fx, &dfx);
             } while (fx < 0.);
         }
         else
         {			/* fix left; search right. */
             do {
                 right += 0.1;
                 Lawless416(x, c, n, right, &fx, &dfx);
                 if (right > 100.) {
                     printf(("EVDMaxLikelyFit(): failed to bracket root"));
                     return 0;
                 }
             } while (fx > 0.);
         }
         /* now we bisection search in left/right interval */
         for (i = 0; i < 100; i++)
         {
             mid = (left + right) / 2.;
             Lawless416(x, c, n, mid, &fx, &dfx);
             if (fabs(fx) < tol) break;             /* success */
             if (fx > 0.)	left = mid;
             else          right = mid;
         }
         if (i == 100) {
             printf(("EVDMaxLikelyFit(): even the bisection search failed"));
             return 0;
         }
         lambda = mid;
     }
 
     /* 3. Substitute into Lawless 4.1.5 to find mu
      */
     esum = 0.;
     total = 0.;
     for (i = 0; i < n; i++)
     {
         mult   = (c == NULL) ? 1. : (double) c[i];
         esum  += mult * exp(-1 * lambda * x[i]);
         total += mult;
     }
     mu = -1. * log(esum / total) / lambda;
 
     *ret_lambda = lambda;
     *ret_mu     = mu;
     return 1;
 }

// Main function
int main(int argv, char **argc) {
    // take in the first argument (filename)
    if (argv < 2) {
        printf("Usage: %s <filename>\n", argc[0]);
        return 1;
    }
    const char *filename = argc[1];
    const char *filename_out = argc[2];
    // open the file
    std::ifstream file(filename);
    if (!file) {
        printf("Error opening file: %s\n", filename);
        return 1;
    }
    // read the file into a vector of floats
    std::vector<float> data;
    std::string line;
    while (std::getline(file, line)) {
        data.push_back(std::stof(line));
    }
    file.close();
    // Sort the data numerically
    std::sort(data.begin(), data.end());
    // convert the vector to an array
    float *data_array = new float[data.size()];
    for (size_t i = 0; i < data.size(); i++) {
        data_array[i] = data[i];
    }
    // call EVDMaxLikelyFit
    float mu = 0.0;
    float lambda = 0.0;
    if (!EVDMaxLikelyFit(data_array, NULL, data.size(), &mu, &lambda)) {
        printf("EVDMaxLikelyFit failed\n");
        delete[] data_array;
        return 1;
    } else {
        printf("mu: %f, lambda: %f\n", mu, lambda);
        std::ofstream file_out;
        file_out.open(filename_out, std::ios::app);
        if (!file_out) {
            printf("Error opening file: %s\n", filename_out);
            delete[] data_array;
            return 1;
        }
        file_out << mu << '\t' << lambda << std::endl;
        file_out.close();
    }
    return 0;
}