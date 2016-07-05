#include "generateIncrements.h"

namespace janek {
    boost::multi_array<double, 2> generateIncrements(unsigned long N, double &D, double &tau, double &alpha, Random *random) {
        double ran_norm = random->normal();
        double ran_uni = random->uniform();
        //boost::multi_array<double, 2> increments;
        int dimensions = 3;
        auto Nextended = N * 2;
        typedef boost::multi_array<double, 2> array_type;
        typedef array_type::index index;
        array_type increments(boost::extents[dimensions][N]);
        fftw_plan planforward, planbackward;
        std::complex<double> *innew = new std::complex<double>[Nextended];
        std::complex<double> *in = new std::complex<double>[Nextended];
        std::complex<double> *out = new std::complex<double>[Nextended];
        std::complex<double> *outnew = new std::complex<double>[Nextended];
        planforward = fftw_plan_dft_1d((int) Nextended, reinterpret_cast<fftw_complex *>(in),
                                       reinterpret_cast<fftw_complex *>(out), FFTW_FORWARD, FFTW_ESTIMATE);
        planbackward = fftw_plan_dft_1d((int) Nextended, reinterpret_cast<fftw_complex *>(innew),
                                        reinterpret_cast<fftw_complex *>(outnew), FFTW_BACKWARD, FFTW_ESTIMATE);

        std::complex<double> icomplex(0, 1);
        double factor = sqrt(D * std::pow((double) N, alpha) / (N));

        for (int ii = 0; ii < Nextended; ++ii) {
            if (ii <= N) {
                in[ii] = std::pow(tau, alpha) * (1 - std::pow((double) ii / N, alpha)) / 2;
            }

            if (ii > N) {
                in[ii] = in[2 * N - (ii)];
            }
        }
        fftw_execute(planforward);
        innew[0] = std::complex<double>(0, 0);

        for (index idim = 0; idim < dimensions; ++idim) {
            //std::cout << idim;
            for (auto i = 0; i < Nextended; ++i) {
                if (i > 0 and i < N) {
                    innew[i] = sqrt(out[i]) * ran_norm * std::exp(2 * M_PI * ran_uni * icomplex);
                }

                if (i >= N) {
                    innew[i] = std::conj(innew[2 * N - i]);
                }
            }
            fftw_execute(planbackward);

            for (index i = 0; i < N; ++i) {
                increments[idim][i] = (-outnew[i].real() + outnew[i + 1].real()) * factor;
            }
        }
        fftw_destroy_plan(planforward);
        fftw_destroy_plan(planbackward);

        return increments;

    }
}
