// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT license.

#pragma once

#include "seal/ciphertext.h"
#include "seal/context.h"
#include "seal/encryptionparams.h"
#include "seal/publickey.h"
#include "seal/randomgen.h"
#include "seal/secretkey.h"
#include <cstdint>
#include <vector>

using namespace std;

namespace seal
{
    namespace util
    {
        /**
        Generate a uniform ternary polynomial and store in RNS representation.

        @param[in] prng A uniform random generator
        @param[in] parms EncryptionParameters used to parameterize an RNS polynomial
        @param[out] destination Allocated space to store a random polynomial
        */
        void sample_poly_ternary(
            std::shared_ptr<UniformRandomGenerator> prng, const EncryptionParameters &parms,
            std::uint64_t *destination);

        // Added by Dice15
        /**
        Generate a sparse ternary polynomial with a specified Hamming weight
        and store it in RNS representation.

        @param[in] prng A uniform random generator
        @param[in] parms EncryptionParameters used to parameterize an RNS polynomial
        @param[out] destination Allocated space to store the random polynomial
        @param[in] hwt The Hamming weight (number of non-zero coefficients) to use
        */
        void sample_poly_ternary_hwt(
            std::shared_ptr<UniformRandomGenerator> prng, const EncryptionParameters &parms, std::uint64_t *destination,
            size_t hwt);

        /**
        Generate a polynomial from a normal distribution and store in RNS representation.
        Samples a noise polynomial from a discrete Gaussian distribution and adds
        the modulo q_i reduced noise to the destination polynomial in RNS form.

        If noise_standard_deviation < 2^53, the function uses the default SEAL
        ClippedNormalDistribution for maximum performance.
        If noise_standard_deviation >= 2^53, it applies exact integer scaling with a
        2-Gaussian convolution (y0 * sigma + y1 * 2^40) to prevent precision loss in
        lower bits caused by floating-point mantissa limitations (IEEE 754).

        @param[in] prng A uniform random generator
        @param[in] parms EncryptionParameters used to parameterize an RNS polynomial
        @param[out] destination Allocated space to store a random polynomial
        @param[in] noise_standard_deviation Standard deviation of the Gaussian distribution.
        */
        void sample_poly_normal(
            std::shared_ptr<UniformRandomGenerator> prng, const EncryptionParameters &parms, std::uint64_t *destination,
            double noise_standard_deviation = global_variables::noise_standard_deviation);

        /**
        Generate a polynomial from a centered binomial distribution and store in RNS representation.

        @param[in] prng A uniform random generator.
        @param[in] parms EncryptionParameters used to parameterize an RNS polynomial
        @param[out] destination Allocated space to store a random polynomial
        */
        void sample_poly_cbd(
            std::shared_ptr<UniformRandomGenerator> prng, const EncryptionParameters &parms, std::uint64_t *destination,
            double noise_standard_deviation = global_variables::noise_standard_deviation);

        /**
        Generate a uniformly random polynomial and store in RNS representation.

        @param[in] prng A uniform random generator
        @param[in] parms EncryptionParameters used to parameterize an RNS polynomial
        @param[out] destination Allocated space to store a random polynomial
        */
        void sample_poly_uniform(
            std::shared_ptr<UniformRandomGenerator> prng, const EncryptionParameters &parms,
            std::uint64_t *destination);

        /**
        Generate a uniformly random polynomial and store in RNS representation.
        This implementation corresponds to Microsoft SEAL 3.4 and earlier.

        @param[in] prng A uniform random generator
        @param[in] parms EncryptionParameters used to parameterize an RNS polynomial
        @param[out] destination Allocated space to store a random polynomial
        */
        void sample_poly_uniform_seal_3_4(
            std::shared_ptr<UniformRandomGenerator> prng, const EncryptionParameters &parms,
            std::uint64_t *destination);

        /**
        Generate a uniformly random polynomial and store in RNS representation.
        This implementation corresponds to Microsoft SEAL 3.5 and earlier.

        @param[in] prng A uniform random generator
        @param[in] parms EncryptionParameters used to parameterize an RNS polynomial
        @param[out] destination Allocated space to store a random polynomial
        */
        void sample_poly_uniform_seal_3_5(
            std::shared_ptr<UniformRandomGenerator> prng, const EncryptionParameters &parms,
            std::uint64_t *destination);

        /**
        Modified by Dice15.
        Create an encryption of zero with a public key and store in a ciphertext.

        @param[in] public_key The public key used for encryption
        @param[in] context The SEALContext containing a chain of ContextData
        @param[in] parms_id Indicates the level of encryption
        @param[in] is_ntt_form If true, store ciphertext in NTT form
        @param[out] destination The output ciphertext - an encryption of zero
        */
        void encrypt_zero_asymmetric(
            const PublicKey &public_key, const SEALContext &context, parms_id_type parms_id, bool is_ntt_form,
            Ciphertext &destination,
            vector<double> noise_standard_deviations =
                vector<double>(2, util::global_variables::noise_standard_deviation));

        /**
        Modified by Dice15.
        Create an encryption of zero with a secret key and store in a ciphertext.

        @param[out] destination The output ciphertext - an encryption of zero
        @param[in] secret_key The secret key used for encryption
        @param[in] context The SEALContext containing a chain of ContextData
        @param[in] parms_id Indicates the level of encryption
        @param[in] is_ntt_form If true, store ciphertext in NTT form
        @param[in] save_seed If true, the second component of ciphertext is
        replaced with the random seed used to sample this component
        */
        void encrypt_zero_symmetric(
            const SecretKey &secret_key, const SEALContext &context, parms_id_type parms_id, bool is_ntt_form,
            bool save_seed, Ciphertext &destination,
            vector<double> noise_standard_deviations =
                vector<double>(1, util::global_variables::noise_standard_deviation));
    } // namespace util
} // namespace seal
