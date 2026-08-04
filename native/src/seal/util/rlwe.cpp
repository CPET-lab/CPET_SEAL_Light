// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT license.

#include "seal/ciphertext.h"
#include "seal/randomgen.h"
#include "seal/randomtostd.h"
#include "seal/util/clipnormal.h"
#include "seal/util/common.h"
#include "seal/util/globals.h"
#include "seal/util/ntt.h"
#include "seal/util/polyarithsmallmod.h"
#include "seal/util/polycore.h"
#include "seal/util/rlwe.h"

using namespace std;

namespace seal
{
    namespace util
    {
        void sample_poly_ternary(
            shared_ptr<UniformRandomGenerator> prng, const EncryptionParameters &parms, uint64_t *destination)
        {
            auto coeff_modulus = parms.coeff_modulus();
            size_t coeff_modulus_size = coeff_modulus.size();
            size_t coeff_count = parms.poly_modulus_degree();

            RandomToStandardAdapter engine(prng);
            uniform_int_distribution<uint64_t> dist(0, 2);

            SEAL_ITERATE(iter(destination), coeff_count, [&](auto &I) {
                uint64_t rand = dist(engine);
                uint64_t flag = static_cast<uint64_t>(-static_cast<int64_t>(rand == 0));
                SEAL_ITERATE(
                    iter(StrideIter<uint64_t *>(&I, coeff_count), coeff_modulus), coeff_modulus_size,
                    [&](auto J) { *get<0>(J) = rand + (flag & get<1>(J).value()) - 1; });
            });
        }

        void sample_poly_ternary_hwt(
            shared_ptr<UniformRandomGenerator> prng, const EncryptionParameters &parms, uint64_t *destination,
            size_t hwt)
        {
            auto coeff_modulus = parms.coeff_modulus();
            size_t coeff_modulus_size = coeff_modulus.size();
            size_t coeff_count = parms.poly_modulus_degree();

            RandomToStandardAdapter engine(prng);
            uniform_int_distribution<uint64_t> dist(0, 1);

            vector<size_t> indices(coeff_count);
            iota(indices.begin(), indices.end(), 0);
            shuffle(indices.begin(), indices.end(), mt19937(random_device{}()));

            vector<bool> sparse(coeff_count, false);
            for (size_t i = 0; i < hwt; ++i)
            {
                sparse[indices[i]] = true;
            }

            size_t idx = 0;
            SEAL_ITERATE(iter(destination), coeff_count, [&](auto &I) {
                uint64_t rand = sparse[idx++] ? (dist(engine) ? 2ULL : 0ULL) : 1ULL;
                uint64_t flag = static_cast<uint64_t>(-static_cast<int64_t>(rand == 0));

                SEAL_ITERATE(
                    iter(StrideIter<uint64_t *>(&I, coeff_count), coeff_modulus), coeff_modulus_size,
                    [&](auto J) { *get<0>(J) = rand + (flag & get<1>(J).value()) - 1; });
            });
        }

        void sample_poly_normal(
            shared_ptr<UniformRandomGenerator> prng, const EncryptionParameters &parms, uint64_t *destination,
            double noise_standard_deviation, uint64_t inverse_scale_factor)
        {
            auto coeff_modulus = parms.coeff_modulus();
            size_t coeff_modulus_size = coeff_modulus.size();
            size_t coeff_count = parms.poly_modulus_degree();
            double noise_standard_deviation_tilde =
                noise_standard_deviation * static_cast<double>(inverse_scale_factor);
            double noise_max_deviation =
                noise_standard_deviation_tilde * global_variables::noise_distribution_width_multiplier;

            if (are_close(noise_max_deviation, 0.0))
            {
                set_zero_poly(coeff_count, coeff_modulus_size, destination);
                return;
            }

            // IEEE 754 double-precision float mantissa limit (2^53 = 9007199254740992.0).
            constexpr double DOUBLE_PRECISION_LIMIT = 9007199254740992.0;

            if (noise_standard_deviation_tilde < DOUBLE_PRECISION_LIMIT)
            {
                // Standard mode: Use default SEAL ClippedNormalDistribution for small standard deviations.
                RandomToStandardAdapter engine(prng);
                ClippedNormalDistribution dist(0, noise_standard_deviation_tilde, noise_max_deviation);

                SEAL_ITERATE(iter(destination), coeff_count, [&](auto &I) {
                    int64_t noise = std::llround(dist(engine) / static_cast<double>(inverse_scale_factor));
                    uint64_t flag = static_cast<uint64_t>(-static_cast<int64_t>(noise < 0));
                    SEAL_ITERATE(
                        iter(StrideIter<uint64_t *>(&I, coeff_count), coeff_modulus), coeff_modulus_size,
                        [&](auto J) { *get<0>(J) = static_cast<uint64_t>(noise) + (flag & get<1>(J).value()); });
                });
            }
            else
            {
                // High-precision mode: Use 2-Gaussian convolution with exact fixed-point integer scaling.
                // This mode prevents precision loss for large standard deviations (>= 2^53).

                // Generate a pair of independent standard normal random variables N(0, 1)
                // using the Box-Muller transform over a 53-bit uniform random bitstream.
                auto unit_normal_pair = [&]() -> pair<double, double> {
                    uint64_t rand1, rand2;
                    prng->generate(sizeof(uint64_t), reinterpret_cast<seal_byte *>(&rand1));
                    prng->generate(sizeof(uint64_t), reinterpret_cast<seal_byte *>(&rand2));

                    // Map uniform 53-bit integers to double in (0, 1] for u1 and [0, 1) for u2.
                    double u1 = (static_cast<double>(rand1 >> 11) + 1.0) * (1.0 / (1ULL << 53));
                    double u2 = static_cast<double>(rand2 >> 11) * (1.0 / (1ULL << 53));

                    double r = sqrt(-2.0 * log(u1));
                    double theta = 2.0 * M_PI * u2;
                    return { r * cos(theta), r * sin(theta) };
                };

                // Convert a double floating-point value y to an exact __int128_t
                // fixed-point integer scaled by 2^(log_scale + F).
                constexpr int F = 50;
                auto exact_scaled = [F](double y, int log_scale) -> __int128_t {
                    int e;
                    double fr = frexp(y, &e);
                    int64_t m = static_cast<int64_t>(fr * 9007199254740992.0); // Exact 53-bit mantissa

                    int sh = e - 53 + log_scale + F;
                    __int128_t val = static_cast<__int128_t>(m);

                    return (sh >= 0) ? (val << sh) : 0;
                };

                __int128_t sigma_128 = static_cast<__int128_t>(noise_standard_deviation_tilde);

                SEAL_ITERATE(iter(destination), coeff_count, [&](auto &I) {
                    double y0, y1;

                    // Bulk tail cut: Resample outliers where |y0| > 10.5.
                    do
                    {
                        auto p = unit_normal_pair();
                        y0 = p.first;
                        y1 = p.second;
                    } while (abs(y0) > 10.5);

                    // Compute primary Gaussian A = sigma * y0 * 2^F.
                    __int128_t A = sigma_128 * exact_scaled(y0, 0);

                    // Compute auxiliary Gaussian B = 2^50 * y1 * 2^F to cover the lower bits.
                    // A 50-bit shift guarantees continuous precision for sigma up to 2^103.
                    __int128_t B = exact_scaled(y1, 50);

                    // Combine primary and auxiliary Gaussians.
                    __int128_t S = A + B;

                    // Compute total divisor combining scaling factor 2^F and inverse_scale_factor (t).
                    __int128_t divisor = static_cast<__int128_t>(inverse_scale_factor) << F;
                    __int128_t half_divisor = divisor >> 1;

                    // Apply exact round-to-nearest division to eliminate C++ integer truncation bias.
                    __int128_t rounded_noise =
                        (S >= 0) ? ((S + half_divisor) / divisor) : ((S - half_divisor) / divisor);

                    // Perform modulo reduction for each RNS prime q_i.
                    SEAL_ITERATE(
                        iter(StrideIter<uint64_t *>(&I, coeff_count), coeff_modulus), coeff_modulus_size, [&](auto J) {
                            uint64_t q_i = get<1>(J).value();

                            // 128-bit modulo reduction with negative result correction.
                            int64_t reduced = static_cast<int64_t>(rounded_noise % static_cast<__int128_t>(q_i));
                            if (reduced < 0)
                            {
                                reduced += q_i;
                            }

                            *get<0>(J) = static_cast<uint64_t>(reduced);
                        });
                });
            }
        }

        void sample_poly_cbd(
            shared_ptr<UniformRandomGenerator> prng, const EncryptionParameters &parms, uint64_t *destination,
            double noise_standard_deviation)
        {
            auto coeff_modulus = parms.coeff_modulus();
            size_t coeff_modulus_size = coeff_modulus.size();
            size_t coeff_count = parms.poly_modulus_degree();

            if (are_close(global_variables::noise_max_deviation, 0.0))
            {
                set_zero_poly(coeff_count, coeff_modulus_size, destination);
                return;
            }

            if (!are_close(noise_standard_deviation, 3.2))
            {
                throw logic_error(
                    "centered binomial distribution only supports standard deviation 3.2; use rounded "
                    "Gaussian instead");
            }

            auto cbd = [&]() {
                unsigned char x[6];
                prng->generate(6, reinterpret_cast<seal_byte *>(x));
                x[2] &= 0x1F;
                x[5] &= 0x1F;
                return hamming_weight(x[0]) + hamming_weight(x[1]) + hamming_weight(x[2]) - hamming_weight(x[3]) -
                       hamming_weight(x[4]) - hamming_weight(x[5]);
            };

            SEAL_ITERATE(iter(destination), coeff_count, [&](auto &I) {
                int32_t noise = cbd();
                uint64_t flag = static_cast<uint64_t>(-static_cast<int64_t>(noise < 0));
                SEAL_ITERATE(
                    iter(StrideIter<uint64_t *>(&I, coeff_count), coeff_modulus), coeff_modulus_size,
                    [&](auto J) { *get<0>(J) = static_cast<uint64_t>(noise) + (flag & get<1>(J).value()); });
            });
        }

        void sample_poly_uniform(
            shared_ptr<UniformRandomGenerator> prng, const EncryptionParameters &parms, uint64_t *destination)
        {
            // Extract encryption parameters
            auto coeff_modulus = parms.coeff_modulus();
            size_t coeff_modulus_size = coeff_modulus.size();
            size_t coeff_count = parms.poly_modulus_degree();
            size_t dest_byte_count = mul_safe(coeff_modulus_size, coeff_count, sizeof(uint64_t));

            constexpr uint64_t max_random = static_cast<uint64_t>(0xFFFFFFFFFFFFFFFFULL);

            // Fill the destination buffer with fresh randomness
            prng->generate(dest_byte_count, reinterpret_cast<seal_byte *>(destination));

            for (size_t j = 0; j < coeff_modulus_size; j++)
            {
                auto &modulus = coeff_modulus[j];
                uint64_t max_multiple = max_random - barrett_reduce_64(max_random, modulus) - 1;
                transform(destination, destination + coeff_count, destination, [&](uint64_t rand) {
                    // This ensures uniform distribution
                    while (rand >= max_multiple)
                    {
                        prng->generate(sizeof(uint64_t), reinterpret_cast<seal_byte *>(&rand));
                    }
                    return barrett_reduce_64(rand, modulus);
                });
                destination += coeff_count;
            }
        }

        void sample_poly_uniform_seal_3_4(
            shared_ptr<UniformRandomGenerator> prng, const EncryptionParameters &parms, uint64_t *destination)
        {
            // Extract encryption parameters
            auto coeff_modulus = parms.coeff_modulus();
            size_t coeff_modulus_size = coeff_modulus.size();
            size_t coeff_count = parms.poly_modulus_degree();

            RandomToStandardAdapter engine(prng);

            constexpr uint64_t max_random = static_cast<uint64_t>(0x7FFFFFFFFFFFFFFFULL);
            for (size_t j = 0; j < coeff_modulus_size; j++)
            {
                auto &modulus = coeff_modulus[j];
                uint64_t max_multiple = max_random - barrett_reduce_64(max_random, modulus) - 1;
                for (size_t i = 0; i < coeff_count; i++)
                {
                    // This ensures uniform distribution
                    uint64_t rand;
                    do
                    {
                        rand = (static_cast<uint64_t>(engine()) << 31) | (static_cast<uint64_t>(engine()) >> 1);
                    } while (rand >= max_multiple);
                    destination[i + j * coeff_count] = barrett_reduce_64(rand, modulus);
                }
            }
        }

        void sample_poly_uniform_seal_3_5(
            shared_ptr<UniformRandomGenerator> prng, const EncryptionParameters &parms, uint64_t *destination)
        {
            // Extract encryption parameters
            auto coeff_modulus = parms.coeff_modulus();
            size_t coeff_modulus_size = coeff_modulus.size();
            size_t coeff_count = parms.poly_modulus_degree();

            RandomToStandardAdapter engine(prng);

            constexpr uint64_t max_random = static_cast<uint64_t>(0xFFFFFFFFFFFFFFFFULL);
            for (size_t j = 0; j < coeff_modulus_size; j++)
            {
                auto &modulus = coeff_modulus[j];
                uint64_t max_multiple = max_random - barrett_reduce_64(max_random, modulus) - 1;
                for (size_t i = 0; i < coeff_count; i++)
                {
                    // This ensures uniform distribution
                    uint64_t rand;
                    do
                    {
                        rand = (static_cast<uint64_t>(engine()) << 32) | static_cast<uint64_t>(engine());
                    } while (rand >= max_multiple);
                    destination[i + j * coeff_count] = barrett_reduce_64(rand, modulus);
                }
            }
        }

        void encrypt_zero_asymmetric(
            const PublicKey &public_key, const SEALContext &context, parms_id_type parms_id, bool is_ntt_form,
            Ciphertext &destination, vector<double> noise_standard_deviations, vector<uint64_t> inverse_scale_factors)
        {
#ifdef SEAL_DEBUG
            if (!is_valid_for(public_key, context))
            {
                throw invalid_argument("public key is not valid for the encryption parameters");
            }
#endif
            // Verify parameters: size must be public_key.size() + 1 (typically 2 + 1 = 3).
            // Index 0: tau_0 for e'_0
            // Index 1: tau_1 for e'_1
            // Index 2: tau_2 for e'_2 (acting as 'u')
            if (noise_standard_deviations.size() != (public_key.data().size() + 1))
            {
                throw invalid_argument("noise_standard_deviations size must match public key size + 1 (typically 3)");
            }
            if (inverse_scale_factors.size() != (public_key.data().size() + 1))
            {
                throw invalid_argument("inverse_scale_factors size must match public key size + 1 (typically 3)");
            }

            // We use a fresh memory pool with `clear_on_destruction' enabled
            MemoryPoolHandle pool = MemoryManager::GetPool(mm_prof_opt::mm_force_new, true);

            auto &context_data = *context.get_context_data(parms_id);
            auto &parms = context_data.parms();
            auto &coeff_modulus = parms.coeff_modulus();
            auto &plain_modulus = parms.plain_modulus();
            size_t coeff_modulus_size = coeff_modulus.size();
            size_t coeff_count = parms.poly_modulus_degree();
            auto ntt_tables = context_data.small_ntt_tables();
            size_t encrypted_size = public_key.data().size();
            scheme_type type = parms.scheme();

            // Make destination have right size and parms_id
            // Ciphertext (c_0,c_1, ...)
            destination.resize(context, parms_id, encrypted_size);
            destination.is_ntt_form() = is_ntt_form;
            destination.scale() = 1.0;
            destination.correction_factor() = 1;

            // c[j] = public_key[j] * u + e[j] in BFV/CKKS = public_key[j] * u + p * e[j] in BGV
            // where e[j] <-- chi, u <-- R_3

            // Create a PRNG; u and the noise/error share the same PRNG
            auto prng = parms.random_generator()->create();

            // Generate u (acting as e'_2) using noise_standard_deviations[2]
            auto u(allocate_poly(coeff_count, coeff_modulus_size, pool));
            if (are_close(noise_standard_deviations[2], 3.2) && inverse_scale_factors[2] == 1)
            {
                sample_poly_cbd(prng, parms, u.get(), noise_standard_deviations[2]);
            }
            else
            {
                sample_poly_normal(prng, parms, u.get(), noise_standard_deviations[2], inverse_scale_factors[2]);
            }
            // sample_poly_ternary(prng, parms, u.get());

            // c[j] = u * public_key[j]
            for (size_t i = 0; i < coeff_modulus_size; i++)
            {
                ntt_negacyclic_harvey(u.get() + i * coeff_count, ntt_tables[i]);
                for (size_t j = 0; j < encrypted_size; j++)
                {
                    dyadic_product_coeffmod(
                        u.get() + i * coeff_count, public_key.data().data(j) + i * coeff_count, coeff_count,
                        coeff_modulus[i], destination.data(j) + i * coeff_count);

                    // Addition with e_0, e_1 is in non-NTT form
                    if (!is_ntt_form)
                    {
                        inverse_ntt_negacyclic_harvey(destination.data(j) + i * coeff_count, ntt_tables[i]);
                    }
                }
            }

            // Generate e_j <-- chi
            // c[j] = public_key[j] * u + e[j] in BFV/CKKS, = public_key[j] * u + p * e[j] in BGV,
            for (size_t j = 0; j < encrypted_size; j++)
            {
                if (are_close(noise_standard_deviations[j], 3.2) && inverse_scale_factors[j] == 1)
                {
                    sample_poly_cbd(prng, parms, u.get(), noise_standard_deviations[j]);
                }
                else
                {
                    sample_poly_normal(prng, parms, u.get(), noise_standard_deviations[j], inverse_scale_factors[j]);
                }
                // SEAL_NOISE_SAMPLER(prng, parms, u.get(), noise_standard_deviation);
                RNSIter gaussian_iter(u.get(), coeff_count);

                // In BGV, p * e is used
                if (type == scheme_type::bgv)
                {
                    if (is_ntt_form)
                    {
                        ntt_negacyclic_harvey_lazy(gaussian_iter, coeff_modulus_size, ntt_tables);
                    }
                    multiply_poly_scalar_coeffmod(
                        gaussian_iter, coeff_modulus_size, plain_modulus.value(), coeff_modulus, gaussian_iter);
                }
                else
                {
                    if (is_ntt_form)
                    {
                        ntt_negacyclic_harvey(gaussian_iter, coeff_modulus_size, ntt_tables);
                    }
                }
                RNSIter dst_iter(destination.data(j), coeff_count);
                add_poly_coeffmod(gaussian_iter, dst_iter, coeff_modulus_size, coeff_modulus, dst_iter);
            }
        }

        void encrypt_zero_symmetric(
            const SecretKey &secret_key, const SEALContext &context, parms_id_type parms_id, bool is_ntt_form,
            bool save_seed, Ciphertext &destination, vector<double> noise_standard_deviations,
            vector<uint64_t> inverse_scale_factors)
        {
#ifdef SEAL_DEBUG
            if (!is_valid_for(secret_key, context))
            {
                throw invalid_argument("secret key is not valid for the encryption parameters");
            }
#endif
            // Verify parameters.
            if (noise_standard_deviations.size() != 1)
            {
                throw invalid_argument("noise_standard_deviations size must be exactly 1 for symmetric encryption");
            }
            if (inverse_scale_factors.size() != 1)
            {
                throw invalid_argument("inverse_scale_factors size must be exactly 1 for symmetric encryption");
            }

            // We use a fresh memory pool with `clear_on_destruction' enabled.
            MemoryPoolHandle pool = MemoryManager::GetPool(mm_prof_opt::mm_force_new, true);

            auto &context_data = *context.get_context_data(parms_id);
            auto &parms = context_data.parms();
            auto &coeff_modulus = parms.coeff_modulus();
            auto &plain_modulus = parms.plain_modulus();
            size_t coeff_modulus_size = coeff_modulus.size();
            size_t coeff_count = parms.poly_modulus_degree();
            auto ntt_tables = context_data.small_ntt_tables();
            size_t encrypted_size = 2;
            scheme_type type = parms.scheme();

            // If a polynomial is too small to store UniformRandomGeneratorInfo,
            // it is best to just disable save_seed. Note that the size needed is
            // the size of UniformRandomGeneratorInfo plus one (uint64_t) because
            // of an indicator word that indicates a seeded ciphertext.
            size_t poly_uint64_count = mul_safe(coeff_count, coeff_modulus_size);
            size_t prng_info_byte_count =
                static_cast<size_t>(UniformRandomGeneratorInfo::SaveSize(compr_mode_type::none));
            size_t prng_info_uint64_count =
                divide_round_up(prng_info_byte_count, static_cast<size_t>(bytes_per_uint64));
            if (save_seed && poly_uint64_count < prng_info_uint64_count + 1)
            {
                save_seed = false;
            }

            destination.resize(context, parms_id, encrypted_size);
            destination.is_ntt_form() = is_ntt_form;
            destination.scale() = 1.0;
            destination.correction_factor() = 1;

            // Create an instance of a random number generator. We use this for sampling
            // a seed for a second PRNG used for sampling u (the seed can be public
            // information. This PRNG is also used for sampling the noise/error below.
            auto bootstrap_prng = parms.random_generator()->create();

            // Sample a public seed for generating uniform randomness
            prng_seed_type public_prng_seed;
            bootstrap_prng->generate(prng_seed_byte_count, reinterpret_cast<seal_byte *>(public_prng_seed.data()));

            // Set up a new default PRNG for expanding u from the seed sampled above
            auto ciphertext_prng = UniformRandomGeneratorFactory::DefaultFactory()->create(public_prng_seed);

            // Generate ciphertext: (c[0], c[1]) = ([-(as+ e)]_q, a) in BFV/CKKS
            // Generate ciphertext: (c[0], c[1]) = ([-(as+pe)]_q, a) in BGV
            uint64_t *c0 = destination.data();
            uint64_t *c1 = destination.data(1);

            // Sample a uniformly at random
            if (is_ntt_form || !save_seed)
            {
                // Sample the NTT form directly
                sample_poly_uniform(ciphertext_prng, parms, c1);
            }
            else if (save_seed)
            {
                // Sample non-NTT form and store the seed
                sample_poly_uniform(ciphertext_prng, parms, c1);
                for (size_t i = 0; i < coeff_modulus_size; i++)
                {
                    // Transform the c1 into NTT representation
                    ntt_negacyclic_harvey(c1 + i * coeff_count, ntt_tables[i]);
                }
            }

            // Sample e <-- chi
            auto noise(allocate_poly(coeff_count, coeff_modulus_size, pool));
            if (are_close(noise_standard_deviations[0], 3.2) && inverse_scale_factors[0] == 1)
            {
                sample_poly_cbd(bootstrap_prng, parms, noise.get(), noise_standard_deviations[0]);
            }
            else
            {
                sample_poly_normal(
                    bootstrap_prng, parms, noise.get(), noise_standard_deviations[0], inverse_scale_factors[0]);
            }
            // SEAL_NOISE_SAMPLER(bootstrap_prng, parms, noise.get(), noise_standard_deviation);

            // Calculate -(as+ e) (mod q) and store in c[0] in BFV/CKKS
            // Calculate -(as+pe) (mod q) and store in c[0] in BGV
            for (size_t i = 0; i < coeff_modulus_size; i++)
            {
                dyadic_product_coeffmod(
                    secret_key.data().data() + i * coeff_count, c1 + i * coeff_count, coeff_count, coeff_modulus[i],
                    c0 + i * coeff_count);
                if (is_ntt_form)
                {
                    // Transform the noise e into NTT representation
                    ntt_negacyclic_harvey(noise.get() + i * coeff_count, ntt_tables[i]);
                }
                else
                {
                    inverse_ntt_negacyclic_harvey(c0 + i * coeff_count, ntt_tables[i]);
                }

                if (type == scheme_type::bgv)
                {
                    // noise = pe instead of e in BGV
                    multiply_poly_scalar_coeffmod(
                        noise.get() + i * coeff_count, coeff_count, plain_modulus.value(), coeff_modulus[i],
                        noise.get() + i * coeff_count);
                }

                // c0 = as + noise
                add_poly_coeffmod(
                    noise.get() + i * coeff_count, c0 + i * coeff_count, coeff_count, coeff_modulus[i],
                    c0 + i * coeff_count);
                // (as + noise, a) -> (-(as + noise), a),
                negate_poly_coeffmod(c0 + i * coeff_count, coeff_count, coeff_modulus[i], c0 + i * coeff_count);
            }

            if (!is_ntt_form && !save_seed)
            {
                for (size_t i = 0; i < coeff_modulus_size; i++)
                {
                    // Transform the c1 into non-NTT representation
                    inverse_ntt_negacyclic_harvey(c1 + i * coeff_count, ntt_tables[i]);
                }
            }

            if (save_seed)
            {
                UniformRandomGeneratorInfo prng_info = ciphertext_prng->info();

                // Write prng_info to destination.data(1) after an indicator word
                c1[0] = static_cast<uint64_t>(0xFFFFFFFFFFFFFFFFULL);
                prng_info.save(reinterpret_cast<seal_byte *>(c1 + 1), prng_info_byte_count, compr_mode_type::none);
            }
        }
    } // namespace util
} // namespace seal
