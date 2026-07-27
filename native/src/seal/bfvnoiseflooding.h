#pragma once

#include "seal/util/polyarithsmallmod.h"
#include "seal/util/timer.h"
#include <cstdint>
#include <seal/seal.h>

using namespace std;
using namespace seal;
using namespace seal::util;

namespace seal
{
    class NoiseFlooding : private Evaluator
    {
    private:
        struct CDTEntry
        {
            int64_t k;
            __uint128_t cdf_val;
        };

        struct CDTTable
        {
            double s;
            double c_tilde;
            int64_t B;
            vector<CDTEntry> table;
        };

    public:
        NoiseFlooding(
            const SEALContext &context, const Encryptor &encryptor, uint64_t a_max, double sigma, double tau0,
            double tau1, size_t b = 0)
            : Evaluator(context), context_(context), encryptor_(encryptor),
              first_parms_(context.first_context_data()->parms()), a_max_(a_max),
              t_(first_parms_.plain_modulus().value()), sigma_(sigma), tau0_(tau0), tau1_(tau1), b_(b)
        {
            // Verify parameters.
            if (first_parms_.scheme() != scheme_type::bfv)
            {
                throw invalid_argument("unsupported scheme");
            }

            // Determine the number of CDT tables based on the maximum offset of the center.
            // The upper bound M for rejection sampling is approx. exp(12 * pi * |c_max| / s).
            // Assuming a typical FHE setting (e.g., s = sigma / t ~ 0.05), setting
            // |c_max| = 1e-4 yields M ~ 1.078, guaranteeing an acceptance rate of ~92.7%.
            // Therefore, if |c_max| <= 1e-4, we use a single table (b = 0) to save memory
            // while maintaining high efficiency. Otherwise (e.g., Batch mode), we generate
            // quantized CDT tables (2^b + 1 tables) with b = 8 to prevent the rejection
            // rate from exploding.
            if (b_ == 0)
            {
                double c_max_ = static_cast<double>(a_max_) / static_cast<double>(t_);
                b_ = (abs(c_max_) > 1e-4) ? 8 : 0;
            }

            make_cdt_tables();
        }

        /**
        Samples a discrete Gaussian error over the coset (coeff + t * Z) and returns the noisy coefficient.
        The function internally centers the input coefficient, selects the appropriate quantized CDT table,
        and applies rejection sampling to perfectly correct the residual offset.
        */
        int64_t sample_coset_discrete_gaussian(uint64_t coeff) const
        {
            // Center the plaintext coefficient to the range [-t/2, t/2).
            // Since the centered value can be negative, we must use int64_t.
            int64_t a = static_cast<int64_t>(coeff);
            if (a >= static_cast<int64_t>(t_ >> 1))
            {
                a -= static_cast<int64_t>(t_);
            }

            double s = sigma_ / static_cast<double>(t_);
            double c = -static_cast<double>(a) / static_cast<double>(t_);
            double c_tilde = 0.0;

            // Find the closest quantized table if multiple tables are used (b_ > 0).
            size_t table_idx = 0;
            if (b_ > 0)
            {
                double scaled_c = c * static_cast<double>(1ULL << b_);

                // idx and max_idx represent relative offsets and can be negative (e.g., [-128, 128]),
                // so they must be kept as int64_t to prevent underflow during clamp.
                int64_t idx = static_cast<int64_t>(round(scaled_c));
                int64_t max_idx = 1LL << (b_ - 1);
                idx = clamp(idx, -max_idx, max_idx);
                c_tilde = static_cast<double>(idx) / static_cast<double>(1ULL << b_);

                // Shift the index to a strictly non-negative range [0, 2^b] for safe array access.
                table_idx = static_cast<size_t>(idx + max_idx);
            }

            const CDTTable &cdt = cdts_[table_idx];

            // Calculate the residual offset (delta) between the actual center and the quantized center.
            // Determine the maximum possible bound for delta.
            double delta = c - c_tilde;
            double delta_max = (b_ > 0) ? (1.0 / static_cast<double>(1ULL << (b_ + 1)))
                                        : (static_cast<double>(a_max_) / static_cast<double>(t_));

            // Calculate the upper bound M for the rejection sampling acceptance probability:
            // M = exp(pi * (12 * s * |delta_max| + delta_max^2) / s^2)
            double M_exp = M_PI * (12.0 * s * delta_max + delta_max * delta_max) / (s * s);
            double M = exp(M_exp);

            auto prng = first_parms_.random_generator()->create();

            while (true)
            {
                // Generate a 128-bit uniform random integer u.
                uint64_t rand_h, rand_l;
                prng->generate(sizeof(uint64_t), reinterpret_cast<seal_byte *>(&rand_h));
                prng->generate(sizeof(uint64_t), reinterpret_cast<seal_byte *>(&rand_l));
                __uint128_t u = ((__uint128_t)rand_h << 64) | rand_l;

                // Use binary search to find the smallest k0 such that u < CDF[k0].
                auto it =
                    upper_bound(cdt.table.begin(), cdt.table.end(), u, [](__uint128_t val, const CDTEntry &entry) {
                        return val < entry.cdf_val;
                    });

                // k0 represents the sampled error and can be negative, so it must be int64_t.
                int64_t k0 = (it != cdt.table.end()) ? it->k : cdt.B;

                // Compute the acceptance probability based on the residual offset delta.
                double x = M_PI * (2.0 * (k0 - c_tilde) * delta - delta * delta) / (s * s);
                double accept_prob = exp(x) / M;

                // Generate a uniform random double v in [0, 1) and test for acceptance.
                uint64_t rand_v;
                prng->generate(sizeof(uint64_t), reinterpret_cast<seal_byte *>(&rand_v));
                double v = static_cast<double>(rand_v >> 11) * (1.0 / (1ULL << 53));

                if (v < accept_prob)
                {
                    // The noisy coefficient is a + t_ * k0. Both 'a' and 'k0' are signed,
                    // so the resulting int64_t safely preserves the negative value if present.
                    return a + static_cast<int64_t>(t_) * k0;
                }

                // TODO: delete
                counting();
            }
        }

        /**
        Multiplies a Ciphertext with a Plaintext using the noise flooding technique.
        The result is stored in the destination parameter.
        */
        void multiply_plain_with_noise_flooding(
            const Ciphertext &encrypted, const Plaintext &plain, Ciphertext &destination,
            MemoryPoolHandle pool = MemoryManager::GetPool()) const
        {
            destination = encrypted;
            multiply_plain_with_noise_flooding_inplace(destination, plain, pool);
        }

        /**
        Multiplies a Ciphertext with a Plaintext using the noise flooding technique in-place.
        This function generates a temporary zero-encryption mask and samples a discrete Gaussian error
        over the coset for the plaintext multiplication to achieve circuit privacy.
        */
        void multiply_plain_with_noise_flooding_inplace(
            Ciphertext &encrypted, const Plaintext &plain, MemoryPoolHandle pool = MemoryManager::GetPool()) const
        {
            Ciphertext enc_mask;
            encryptor_.encrypt_zero(enc_mask, { tau0_, tau1_ }, pool);
            multiply_coset_plain(encrypted, plain, pool);
            add_inplace(encrypted, enc_mask);
        }

        /**
        Multiplies a Ciphertext with a Plaintext using the noise flooding technique with a provided random plaintext r.
        The result is stored in the destination parameter.
        */
        void multiply_plain_with_noise_flooding(
            const Ciphertext &encrypted, const Plaintext &plain, const Plaintext &r, Ciphertext &destination,
            MemoryPoolHandle pool = MemoryManager::GetPool()) const
        {
            Ciphertext enc_mask;
            destination = encrypted;
            multiply_plain_with_noise_flooding_inplace(destination, plain, r, pool);
        }

        /**
        Multiplies a Ciphertext with a Plaintext using the noise flooding technique in-place,
        using a provided random plaintext r for symmetric encryption of the mask.
        */
        void multiply_plain_with_noise_flooding_inplace(
            Ciphertext &encrypted, const Plaintext &plain, const Plaintext &r,
            MemoryPoolHandle pool = MemoryManager::GetPool()) const
        {
            Ciphertext enc_mask;
            encryptor_.encrypt(r, enc_mask, { tau0_, tau1_ }, pool);
            multiply_coset_plain(encrypted, plain, pool);
            add_inplace(encrypted, enc_mask);
        }

        /**
        Core multiplication logic for noise flooding.
        Samples a discrete Gaussian error for each coefficient of the plaintext over its corresponding coset,
        transforms the noisy plaintext into NTT form, and performs dyadic multiplication with the ciphertext.
        */
        void multiply_coset_plain(
            Ciphertext &encrypted, const Plaintext &plain, MemoryPoolHandle pool = MemoryManager::GetPool()) const
        {
            // Extract encryption parameters.
            auto &context_data = *context_.get_context_data(encrypted.parms_id());
            auto &parms = context_data.parms();
            auto &coeff_modulus = parms.coeff_modulus();
            size_t coeff_count = parms.poly_modulus_degree();
            size_t coeff_modulus_size = coeff_modulus.size();

            if (parms.scheme() != scheme_type::bfv)
            {
                throw invalid_argument("unsupported scheme");
            }

            uint64_t plain_upper_half_threshold = context_data.plain_upper_half_threshold();
            auto plain_upper_half_increment = context_data.plain_upper_half_increment();
            auto ntt_tables = iter(context_data.small_ntt_tables());

            size_t encrypted_size = encrypted.size();
            size_t plain_coeff_count = plain.coeff_count();
            size_t plain_nonzero_coeff_count = plain.nonzero_coeff_count();

            // Size check
            if (!product_fits_in(encrypted_size, coeff_count, coeff_modulus_size))
            {
                throw logic_error("invalid parameters");
            }

            /*
            Optimizations for constant / monomial multiplication can lead to the
            presence of a timing side-channel in use-cases where the plaintext data
            should also be kept private.
            */
            if (plain_nonzero_coeff_count == 1)
            {
                // Multiplying by a monomial?
                size_t mono_exponent = plain.significant_coeff_count() - 1;
                int64_t coset_sample = sample_coset_discrete_gaussian(plain[mono_exponent]);

                SEAL_ALLOCATE_GET_COEFF_ITER(coset_rns, coeff_modulus_size, pool);

                for (size_t i = 0; i < coeff_modulus_size; i++)
                {
                    auto &q_i = coeff_modulus[i];

                    if (coset_sample >= 0)
                    {
                        coset_rns[i] = q_i.reduce(static_cast<uint64_t>(coset_sample));
                    }
                    else
                    {
                        coset_rns[i] = negate_uint_mod(q_i.reduce(static_cast<uint64_t>(-coset_sample)), q_i);
                    }
                }

                negacyclic_multiply_poly_mono_coeffmod(
                    encrypted, encrypted_size, coset_rns, mono_exponent, coeff_modulus, encrypted, pool);

                return;
            }

            // Generic case: any plaintext polynomial
            // Allocate temporary space for an entire RNS polynomial
            auto temp(allocate_zero_poly(coeff_count, coeff_modulus_size, pool));
            RNSIter temp_iter(temp.get(), coeff_count);
            auto coset_samples = allocate<int64_t>(coeff_count, pool);

            for (size_t j = 0; j < coeff_count; j++)
            {
                uint64_t a_j = (j < plain_coeff_count) ? plain[j] : 0;
                coset_samples[j] = sample_coset_discrete_gaussian(static_cast<int64_t>(a_j));
            }

            SEAL_ITERATE(iter(temp_iter, coeff_modulus), coeff_modulus_size, [&](auto I) {
                auto target_poly = get<0>(I);
                auto &q_i = get<1>(I);

                for (size_t j = 0; j < coeff_count; j++)
                {
                    int64_t coset_sample = coset_samples[j];

                    if (coset_sample >= 0)
                    {
                        target_poly[j] = q_i.reduce(static_cast<uint64_t>(coset_sample));
                    }
                    else
                    {
                        target_poly[j] = negate_uint_mod(q_i.reduce(static_cast<uint64_t>(-coset_sample)), q_i);
                    }
                }
            });

            // Need to multiply each component in encrypted with temp; first step is to
            // transform to NTT form
            ntt_negacyclic_harvey(temp_iter, coeff_modulus_size, ntt_tables);

            SEAL_ITERATE(iter(encrypted), encrypted_size, [&](auto I) {
                SEAL_ITERATE(iter(I, temp_iter, coeff_modulus, ntt_tables), coeff_modulus_size, [&](auto J) {
                    // Lazy reduction
                    ntt_negacyclic_harvey_lazy(get<0>(J), get<3>(J));
                    dyadic_product_coeffmod(get<0>(J), get<1>(J), coeff_count, get<2>(J), get<0>(J));
                    inverse_ntt_negacyclic_harvey(get<0>(J), get<3>(J));
                });
            });
        }

    private:
        /**
        Generates Cumulative Distribution Tables (CDTs) for the discrete Gaussian distribution.
        Depending on the parameter b_, it either generates a single table (b_ = 0) centered at 0.0,
        or multiple quantized tables (b_ > 0) to support arbitrary centers in [-0.5, 0.5].
        */
        void make_cdt_tables()
        {
            int64_t num_tables = (b_ == 0) ? 1 : ((1 << b_) + 1);
            cdts_.resize(num_tables);

            if (b_ == 0)
            {
                construct_cdt_table(cdts_[0], 0.0);
            }
            else
            {
                // Ex) For b_ = 8, the offset is 128.
                // We create tables centered at c_tilde = i / 256 for i in [-128, 128].
                int64_t offset = 1 << (b_ - 1);
                for (int64_t i = -offset; i <= offset; ++i)
                {
                    // Array index mapping: i [-128, 128] maps to index [0, 256].
                    double c_tilde = static_cast<double>(i) / static_cast<double>(1ULL << b_);
                    construct_cdt_table(cdts_[static_cast<size_t>(i + offset)], c_tilde);
                }
            }
        }

        /**
        Constructs a single CDT table for a specific center c_tilde.
        The table covers the range [-B, B] where B is bounded by ceil(6 * s) + 1.
        */
        void construct_cdt_table(CDTTable &cdt, double c_tilde) const
        {
            cdt.s = sigma_ / static_cast<double>(t_);
            cdt.c_tilde = c_tilde;

            // Add +1 to account for the center shift c_tilde.
            cdt.B = static_cast<int64_t>(ceil(6.0 * cdt.s)) + 1;

            // k \in {-B, ..., B}. Total entries = 2B + 1.
            size_t num_entries = static_cast<size_t>(2 * cdt.B + 1);
            vector<long double> pdf(num_entries);
            long double total_sum = 0.0L;

            // Calculate the unnormalized Probability Density Function (PDF):
            // p(k) = exp(-pi * (k - c_tilde)^2 / s^2)
            for (int64_t k = -cdt.B; k <= cdt.B; ++k)
            {
                size_t idx = static_cast<size_t>(k + cdt.B);
                long double exponent = -M_PI * pow(static_cast<long double>(k) - c_tilde, 2.0) / (cdt.s * cdt.s);
                pdf[idx] = exp(exponent);
                total_sum += pdf[idx];
            }

            // Calculate the Cumulative Distribution Function (CDF) scaled to 2^128.
            long double cumulative = 0.0L;
            long double scale_64 = pow(2.0L, 64);
            long double scale_128 = pow(2.0L, 128);

            cdt.table.resize(num_entries);

            for (int64_t k = -cdt.B; k <= cdt.B; ++k)
            {
                size_t idx = static_cast<size_t>(k + cdt.B);

                // Force the last cumulative value to 1.0 to prevent floating-point inaccuracies.
                if (k == cdt.B)
                {
                    cumulative = 1.0L;
                }
                else
                {
                    cumulative += pdf[idx] / total_sum;
                }

                long double scaled_val = cumulative * scale_128;
                __uint128_t uint128_val = 0;

                // Handle potential precision overflow at the boundaries.
                if (scaled_val >= scale_128)
                {
                    uint128_val = ~(__uint128_t)0;
                }
                else
                {
                    long double high = floor(scaled_val / scale_64);
                    long double low = scaled_val - high * scale_64;
                    uint128_val = ((__uint128_t)(static_cast<uint64_t>(high)) << 64) | static_cast<uint64_t>(low);
                }

                cdt.table[idx] = { k, uint128_val };
            }
        }

        const SEALContext &context_;

        const Encryptor &encryptor_;

        const EncryptionParameters &first_parms_;

        const uint64_t a_max_;

        const uint64_t t_;

        const double sigma_;

        const double tau0_;

        const double tau1_;

        size_t b_;

        vector<CDTTable> cdts_;
    };
} // namespace seal