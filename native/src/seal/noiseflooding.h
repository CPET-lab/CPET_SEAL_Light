#pragma once

#include <seal/seal.h>
#include <random>

using namespace std;
using namespace seal;
using namespace seal::util;

namespace seal
{
    inline void make_encryption_of_zero(
        const SEALContext& context, const CoeffEncoder& encoder, const Encryptor& encryptor, const Evaluator& evaluator, 
        int log_noise_bound, int stat_sec, Ciphertext& destination)
    {
        // Extract parameters.
        auto context_data_ptr = context.first_context_data();
        const auto& coeff_modulus = context_data_ptr->parms().coeff_modulus();
        size_t coeff_modulus_size = coeff_modulus.size();
        uint64_t n = context_data_ptr->parms().poly_modulus_degree();
        uint64_t t = context_data_ptr->parms().plain_modulus().value();


        // Encrypt zero.
        Plaintext encoded_zero;
        Ciphertext encrypted_zero;
        encoder.encode(vector<uint64_t>(n), encoded_zero);
        encryptor.encrypt(encoded_zero, encrypted_zero);


        // Generate noise vector.   
        random_device rd;
        mt19937 gen(rd());
        vector<vector<pair<uint64_t, uint64_t>>> noises(n);
        vector<pair<bool, bool>> signs(n);

        int log_B = (stat_sec >> 1) + log_noise_bound - 1;
        while (log_B > 0)
        {
            uniform_int_distribution<uint64_t> noise_dist(0, 1ULL << (log_B > 60 ? 60 : log_B));
            bernoulli_distribution sign_dist(0.5);
            for (size_t i = 0; i < n; i++)
            {
                if (noises[i].empty())
                {
                    signs[i] = { sign_dist(gen) , sign_dist(gen) };
                }

                noises[i].emplace_back();
                noises[i].back() = { noise_dist(gen) ,  noise_dist(gen) };
            }

            log_B -= 60;
        }


        // Add large noise to encrypted zero.
        evaluator.transform_from_ntt_inplace(encrypted_zero);
        Ciphertext::ct_coeff_type* ptr = encrypted_zero.data();

        for (size_t j = 0; j < coeff_modulus_size; j++)
        {
            const auto& q_j = coeff_modulus[j];

            for (size_t i = 0; i < n; i++, ptr++)
            {
                uint64_t E0 = 0ULL;
                uint64_t E1 = 0ULL;
                uint64_t acc_multiplier = 1ULL;

                for (const auto& [e0, e1] : noises[i])
                {
                    E0 = add_uint_mod(E0, multiply_uint_mod(e0, acc_multiplier, q_j), q_j);
                    E1 = add_uint_mod(E1, multiply_uint_mod(e1, acc_multiplier, q_j), q_j);
                    acc_multiplier = multiply_uint_mod(acc_multiplier, q_j.value(), q_j);
                }

                if (signs[i].first && E0)
                {
                    E0 = q_j.value() - E0;
                }

                if (signs[i].second && E1)
                {
                    E1 = q_j.value() - E1;
                }

                uint64_t E = add_uint_mod(E0, E1, q_j);
                E = multiply_uint_mod(E, t, q_j);
                *ptr = add_uint_mod(*ptr, E, q_j);
            }
        }

        evaluator.transform_to_ntt(encrypted_zero, destination);
    }


    inline void noise_flooding_inplace(const Evaluator& evaluator, const Ciphertext& encrypted_zero, Ciphertext& encrypted)
    {
        // Verify parameters.
        if (!encrypted.is_ntt_form() || !encrypted_zero.is_ntt_form())
        {
            throw logic_error("Ciphertext must be in NTT form.");
        }


        // Add the noise to the input ciphertext.
        evaluator.add_inplace(encrypted, encrypted_zero);
    }
}