// Added by Dice15.

#include "seal/coeffencoder.h"
#include "seal/valcheck.h"
#include "seal/util/common.h"
#include <algorithm>
#include <limits>
#include <random>
#include <stdexcept>

using namespace std;
using namespace seal::util;

namespace seal
{
    CoeffEncoder::CoeffEncoder(const SEALContext& context) : context_(context)
    {
        // Verify parameters
        if (!context_.parameters_set())
        {
            throw invalid_argument("encryption parameters are not set correctly");
        }

        auto& context_data = *context_.first_context_data();
        if (context_data.parms().scheme() != scheme_type::bfv && context_data.parms().scheme() != scheme_type::bgv)
        {
            throw invalid_argument("unsupported scheme");
        }
        if (!context_data.qualifiers().using_batching)
        {
            throw invalid_argument("encryption parameters are not valid for batching");
        }

        // Set the slot count
        slots_ = context_data.parms().poly_modulus_degree();
    }

    void CoeffEncoder::encode(
        const vector<uint64_t>& values_matrix, Plaintext& destination) const
    {
        auto& context_data = *context_.first_context_data();

        // Validate input parameters
        size_t values_matrix_size = values_matrix.size();
        if (values_matrix_size > slots_)
        {
            throw invalid_argument("values_matrix size is too large");
        }
#ifdef SEAL_DEBUG
        uint64_t modulus = context_data.parms().plain_modulus().value();
        for (auto v : values_matrix)
        {
            // Validate the i-th input
            if (v >= modulus)
            {
                throw invalid_argument("input value is larger than plain_modulus");
            }
        }
#endif
        // Set destination to full size
        destination.resize(slots_);
        destination.parms_id() = parms_id_zero;

        // Write values to destination coefficients directly without mapping.
         // Convolution-based multiplication mode without NTT transformation.
        for (size_t i = 0; i < values_matrix_size; i++)
        {
            *(destination.data() + i) = values_matrix[i];
        }
        for (size_t i = values_matrix_size; i < slots_; i++)
        {
            *(destination.data() + i) = 0;
        }
    }

    void CoeffEncoder::encode(
        const vector<int64_t>& values_matrix, Plaintext& destination) const
    {
        auto& context_data = *context_.first_context_data();
        uint64_t modulus = context_data.parms().plain_modulus().value();

        // Validate input parameters
        size_t values_matrix_size = values_matrix.size();
        if (values_matrix_size > slots_)
        {
            throw invalid_argument("values_matrix size is too large");
        }
#ifdef SEAL_DEBUG
        uint64_t plain_modulus_div_two = modulus >> 1;
        for (auto v : values_matrix)
        {
            // Validate the i-th input
            if (unsigned_gt(llabs(v), plain_modulus_div_two))
            {
                throw invalid_argument("input value is larger than plain_modulus");
            }
        }
#endif
        // Set destination to full size
        destination.resize(slots_);
        destination.parms_id() = parms_id_zero;

        // Write values to destination coefficients directly without mapping.
        //  Convolution-based multiplication mode without NTT transformation.
        for (size_t i = 0; i < values_matrix_size; i++)
        {
            *(destination.data() + i) = (values_matrix[i] < 0) ? (modulus + static_cast<uint64_t>(values_matrix[i]))
                : static_cast<uint64_t>(values_matrix[i]);
        }
        for (size_t i = values_matrix_size; i < slots_; i++)
        {
            *(destination.data() + i) = 0;
        }
    }
#ifdef SEAL_USE_MSGSL
    void CoeffEncoder::encode(
        gsl::span<const uint64_t> values_matrix, Plaintext& destination) const
    {
        auto& context_data = *context_.first_context_data();

        // Validate input parameters
        size_t values_matrix_size = static_cast<size_t>(values_matrix.size());
        if (values_matrix_size > slots_)
        {
            throw invalid_argument("values_matrix size is too large");
        }
#ifdef SEAL_DEBUG
        uint64_t modulus = context_data.parms().plain_modulus().value();
        for (auto v : values_matrix)
        {
            // Validate the i-th input
            if (v >= modulus)
            {
                throw invalid_argument("input value is larger than plain_modulus");
            }
        }
#endif
        // Set destination to full size
        destination.resize(slots_);
        destination.parms_id() = parms_id_zero;

        // Write values to destination coefficients directly without mapping.
        // Convolution-based multiplication mode without NTT transformation.
        for (size_t i = 0; i < values_matrix_size; i++)
        {
            *(destination.data() + i) = values_matrix[i];
        }
        for (size_t i = values_matrix_size; i < slots_; i++)
        {
            *(destination.data() + i) = 0;
        }
    }

    void CoeffEncoder::encode(
        gsl::span<const int64_t> values_matrix, Plaintext& destination) const
    {
        auto& context_data = *context_.first_context_data();
        uint64_t modulus = context_data.parms().plain_modulus().value();

        // Validate input parameters
        size_t values_matrix_size = static_cast<size_t>(values_matrix.size());
        if (values_matrix_size > slots_)
        {
            throw invalid_argument("values_matrix size is too large");
        }
#ifdef SEAL_DEBUG
        uint64_t plain_modulus_div_two = modulus >> 1;
        for (auto v : values_matrix)
        {
            // Validate the i-th input
            if (unsigned_gt(llabs(v), plain_modulus_div_two))
            {
                throw invalid_argument("input value is larger than plain_modulus");
            }
        }
#endif
        // Set destination to full size
        destination.resize(slots_);
        destination.parms_id() = parms_id_zero;

        // Write values to destination coefficients directly without mapping.
        // convolution`: Convolution-based multiplication mode without NTT transformation.
        for (size_t i = 0; i < values_matrix_size; i++)
        {
            *(destination.data() + i) = (values_matrix[i] < 0) ? (modulus + static_cast<uint64_t>(values_matrix[i]))
                : static_cast<uint64_t>(values_matrix[i]);
        }
        for (size_t i = values_matrix_size; i < slots_; i++)
        {
            *(destination.data() + i) = 0;
        }
    }
#endif
    void CoeffEncoder::decode(
        const Plaintext& plain, vector<uint64_t>& destination, MemoryPoolHandle pool) const
    {
        if (!is_valid_for(plain, context_))
        {
            throw invalid_argument("plain is not valid for encryption parameters");
        }
        if (plain.is_ntt_form())
        {
            throw invalid_argument("plain cannot be in NTT form");
        }
        if (!pool)
        {
            throw invalid_argument("pool is uninitialized");
        }

        auto& context_data = *context_.first_context_data();

        // Set destination size
        destination.resize(slots_);

        // Never include the leading zero coefficient (if present)
        size_t plain_coeff_count = min(plain.coeff_count(), slots_);

        auto temp_dest(allocate_uint(slots_, pool));

        // Make a copy of poly
        set_uint(plain.data(), plain_coeff_count, temp_dest.get());
        set_zero_uint(slots_ - plain_coeff_count, temp_dest.get() + plain_coeff_count);

        // Copy coefficients directly without transformation or reordering.
        // Directly copy the coefficients without transformation or reordering.
        for (size_t i = 0; i < slots_; i++)
        {
            destination[i] = temp_dest[i];
        }
    }

    void CoeffEncoder::decode(
        const Plaintext& plain, vector<int64_t>& destination, MemoryPoolHandle pool) const
    {
        if (!is_valid_for(plain, context_))
        {
            throw invalid_argument("plain is not valid for encryption parameters");
        }
        if (plain.is_ntt_form())
        {
            throw invalid_argument("plain cannot be in NTT form");
        }
        if (!pool)
        {
            throw invalid_argument("pool is uninitialized");
        }

        auto& context_data = *context_.first_context_data();
        uint64_t modulus = context_data.parms().plain_modulus().value();

        // Set destination size
        destination.resize(slots_);

        // Never include the leading zero coefficient (if present)
        size_t plain_coeff_count = min(plain.coeff_count(), slots_);

        auto temp_dest(allocate_uint(slots_, pool));

        // Make a copy of poly
        set_uint(plain.data(), plain_coeff_count, temp_dest.get());
        set_zero_uint(slots_ - plain_coeff_count, temp_dest.get() + plain_coeff_count);

        // Copy coefficients directly without transformation or reordering.
        // Directly copy the coefficients without transformation or reordering.
        uint64_t plain_modulus_div_two = modulus >> 1;
        for (size_t i = 0; i < slots_; i++)
        {
            uint64_t curr_value = temp_dest[i];
            destination[i] = (curr_value > plain_modulus_div_two)
                ? (static_cast<int64_t>(curr_value) - static_cast<int64_t>(modulus))
                : static_cast<int64_t>(curr_value);
        }
    }
#ifdef SEAL_USE_MSGSL
    void CoeffEncoder::decode(
        const Plaintext& plain, gsl::span<uint64_t> destination, MemoryPoolHandle pool) const
    {
        if (!is_valid_for(plain, context_))
        {
            throw invalid_argument("plain is not valid for encryption parameters");
        }
        if (plain.is_ntt_form())
        {
            throw invalid_argument("plain cannot be in NTT form");
        }
        if (!pool)
        {
            throw invalid_argument("pool is uninitialized");
        }

        auto& context_data = *context_.first_context_data();

        if (unsigned_gt(destination.size(), numeric_limits<int>::max()) || unsigned_neq(destination.size(), slots_))
        {
            throw invalid_argument("destination has incorrect size");
        }

        // Never include the leading zero coefficient (if present)
        size_t plain_coeff_count = min(plain.coeff_count(), slots_);

        auto temp_dest(allocate_uint(slots_, pool));

        // Make a copy of poly
        set_uint(plain.data(), plain_coeff_count, temp_dest.get());
        set_zero_uint(slots_ - plain_coeff_count, temp_dest.get() + plain_coeff_count);

        // Copy coefficients directly without transformation or reordering.
        // Directly copy the coefficients without transformation or reordering.
        for (size_t i = 0; i < slots_; i++)
        {
            destination[i] = temp_dest[i];
        }
    }

    void CoeffEncoder::decode(
        const Plaintext& plain, gsl::span<int64_t> destination, MemoryPoolHandle pool) const
    {
        if (!is_valid_for(plain, context_))
        {
            throw invalid_argument("plain is not valid for encryption parameters");
        }
        if (plain.is_ntt_form())
        {
            throw invalid_argument("plain cannot be in NTT form");
        }
        if (!pool)
        {
            throw invalid_argument("pool is uninitialized");
        }

        auto& context_data = *context_.first_context_data();
        uint64_t modulus = context_data.parms().plain_modulus().value();

        if (unsigned_gt(destination.size(), numeric_limits<int>::max()) || unsigned_neq(destination.size(), slots_))
        {
            throw invalid_argument("destination has incorrect size");
        }

        // Never include the leading zero coefficient (if present)
        size_t plain_coeff_count = min(plain.coeff_count(), slots_);

        auto temp_dest(allocate_uint(slots_, pool));

        // Make a copy of poly
        set_uint(plain.data(), plain_coeff_count, temp_dest.get());
        set_zero_uint(slots_ - plain_coeff_count, temp_dest.get() + plain_coeff_count);

        // Copy coefficients directly without transformation or reordering.
        // Directly copy the coefficients without transformation or reordering.
        uint64_t plain_modulus_div_two = modulus >> 1;
        for (size_t i = 0; i < slots_; i++)
        {
            uint64_t curr_value = temp_dest[i];
            destination[i] = (curr_value > plain_modulus_div_two)
                ? (static_cast<int64_t>(curr_value) - static_cast<int64_t>(modulus))
                : static_cast<int64_t>(curr_value);
        }
    }
#endif
} // namespace seal